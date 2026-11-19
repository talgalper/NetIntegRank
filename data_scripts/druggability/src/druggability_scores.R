### converts fpocket _info.txt outputs into one druggability table ###
# usage: Rscript druggability_scores.R <alphafold|swissmodel> <results_dir> <threads> <struct_score.csv>
#                                      <structures_dir> <min_residues> <min_conf>
#
# For each structure: the best pocket's scores and the number of druggable pockets, over all pockets and over
# confident pockets only (pocket confidence = mean confidence of the residues lining it), plus its residue count.
# Column names start with AF_ or SM_ (all except uniprot_id) so the two sources merge easily. SWISS-MODEL confidence
# (QMEANDisCo, 0-1) is multiplied by 100 so both sources use one 0-100 scale.
#
# Outputs, in results_dir:
#   fpocket_druggability_full.csv  structures with at least min_residues residues
#   small_structures.csv           the shorter structures, same columns, left out of the main table
#   fpocket_pockets.csv            every pocket of every structure, with its druggability and confidence

suppressPackageStartupMessages({
  library(readr)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: Rscript druggability_scores.R <alphafold|swissmodel> <results_dir> <threads> <struct_score.csv> ",
       "<structures_dir> <min_residues> <min_conf>")
}
source_type <- args[1]
results_dir <- args[2]
threads     <- as.integer(args[3])
struct_score_file <- args[4]   # written by the confidence filter (filter_structures.sh)
structures_dir <- args[5]      # the structures fpocket was run on
min_residues <- as.integer(args[6])   # structures shorter than this go to small_structures.csv
min_conf <- as.numeric(args[7])       # a pocket is confident if its lining residues' mean confidence is at least this
scores_dir  <- file.path(results_dir, "scores")
if (!source_type %in% c("alphafold", "swissmodel")) stop("source must be alphafold or swissmodel")
prefix <- if (source_type == "alphafold") "AF_" else "SM_"
conf_scale <- if (source_type == "alphafold") 1 else 100   # QMEANDisCo 0-1 -> 0-100


## parse one fpocket _info.txt into a descriptors x pockets table
fpocket_format <- function(txt_file){
  # read in .txt file
  df <- read.delim(txt_file, sep = ",", header = F)

  # Create a list to store the new columns
  new_columns <- list()

  # Loop through every 20th row in column 1 and move it to a new column
  for (i in seq(20, nrow(df), 20)) {
    new_col <- df[(i-19):i, 1]
    new_columns[[length(new_columns)+1]] <- new_col
  }

  # Bind the new columns together
  new_df <- do.call(cbind, new_columns)
  new_df <- as.data.frame(new_df)

  # rename columns
  new_df <- new_df[-1,]
  new_df <- as.data.frame(new_df) # make sure its a data frame
  colnames(new_df) <- paste0("pocket_", 1:ncol(new_df))

  # rename rows
  new_df <- cbind(new_df, do.call("rbind", strsplit(as.character(new_df$`pocket_1`), ":")))
  new_df <- new_df[, -ncol(new_df)]
  row_names <- new_df[, ncol(new_df)]
  new_df <- new_df[, -ncol(new_df)]
  new_df <- as.data.frame(new_df) # make sure its a data frame
  colnames(new_df) <- paste0("pocket_", 1:ncol(new_df))
  rownames(new_df) <- row_names

  # loop through every cell in the data frame and remove text before ":"
  for (i in 1:nrow(new_df)) {
    for (j in 1:ncol(new_df)) {
      new_df[i, j] <- gsub(".*:", "", new_df[i, j])
      new_df[i, j] <- gsub("\t", "", new_df[i, j])
    }
  }

  # convert all cells from character to numeric
  for (col in names(new_df)) {
    new_df[[col]] <- as.numeric(new_df[[col]])
  }

  return(new_df)
}


## per-residue confidence of a structure: B-factor of each residue's CA atom, named by chain + residue number
## (pLDDT for AlphaFold; QMEANDisCo local x 100 for SWISS-MODEL)
residue_conf <- function(pdb_file) {
  lines <- readLines(pdb_file)
  ca <- lines[startsWith(lines, "ATOM") & substr(lines, 13, 16) == " CA "]
  setNames(round(as.numeric(substr(ca, 61, 66)) * conf_scale, 6),    # rounded so 0.29 x 100 is exactly 29
           paste(substr(ca, 22, 22), trimws(substr(ca, 23, 27))))
}


## summarise one structure: best pocket's scores and number of druggable pockets,
## then the same counting only confident pockets
summarise_file <- function(file) {
  path <- file.path(scores_dir, file)
  if (file.size(path) == 0) return(NULL)   # empty fpocket output (no pockets), skipped

  name <- sub("_info\\.txt$", "", file)    # the structure file name without .pdb
  if (source_type == "alphafold") {
    split_name <- unlist(strsplit(file, split = "-"))
    structure_id <- paste0(split_name[2], "-", split_name[3])   # AF-P12345-F1-model_v6 -> P12345-F1
    uniprot_id <- sub("^[^-]+-([^-]+)-.*", "\\1", file)
  } else {
    structure_id <- name
    uniprot_id <- name
  }

  ## all pockets
  data <- fpocket_format(path)
  pocket_max <- max(data[1,])
  druggability_max <- max(data[2,])
  num_drug_pockets <- sum(data[2,] >= 0.4)
  best <- which.max(data["\tDruggability Score", ])

  ## pocket confidence: mean confidence of the residues lining each pocket. fpocket lists the lining atoms of
  ## pocket N in structural/<name>_out/pockets/pocketN_atm.pdb; their residues are looked up in the input structure.
  conf <- residue_conf(file.path(structures_dir, paste0(name, ".pdb")))
  pocket_dir <- file.path(results_dir, "structural", paste0(name, "_out"), "pockets")
  n_pockets <- ncol(data)
  n_files <- length(list.files(pocket_dir, pattern = "^pocket[0-9]+_atm\\.pdb$"))
  if (n_files != n_pockets) stop(name, ": ", n_pockets, " pockets in _info.txt but ", n_files, " pocket files")
  pocket_conf <- numeric(n_pockets)
  n_lining <- integer(n_pockets)
  for (k in seq_len(n_pockets)) {
    lines <- readLines(file.path(pocket_dir, paste0("pocket", k, "_atm.pdb")))
    ## check the pocket file is the same pocket as column k of _info.txt (drug score given to 4 and 3 decimals)
    header_drug <- as.numeric(sub(".*:", "", grep("Drug Score", lines, value = TRUE)[1]))
    if (abs(header_drug - data[2, k]) > 0.0006) stop(name, ": pocket ", k, " drug score differs between files")
    atoms <- lines[startsWith(lines, "ATOM")]
    residues <- unique(paste(substr(atoms, 22, 22), trimws(substr(atoms, 23, 27))))
    if (any(!residues %in% names(conf))) stop(name, ": pocket ", k, " has residues not found in the structure")
    # 15 significant digits, as PocketMiner's window averages (saved and read back), so a mean of exactly 70 counts
    pocket_conf[k] <- signif(mean(conf[residues]), 15)
    n_lining[k] <- length(residues)
  }
  confident <- pocket_conf >= min_conf

  row <- data.frame(uniprot_id = uniprot_id,
                    structure_id = structure_id,
                    n_residues = length(conf),
                    pocket = pocket_max,
                    druggability = druggability_max,
                    num_drug_pockets = num_drug_pockets,
                    hydrophobicity_score = data["\tHydrophobicity score", best],
                    volume_score = data["\tVolume score", best],
                    alpha_max_dist = data["\tCent. of mass - Alpha Sphere max dist", best],
                    conf_best_pocket = pocket_conf[best],
                    druggability_conf = if (any(confident)) max(unlist(data[2, confident])) else NA_real_,
                    num_drug_pockets_conf = sum(data[2, ] >= 0.4 & confident))
  pockets <- data.frame(uniprot_id = uniprot_id,
                        structure_id = structure_id,
                        pocket_number = seq_len(n_pockets),
                        pocket_score = unname(unlist(data[1, ])),
                        druggability = unname(unlist(data[2, ])),
                        pocket_conf = pocket_conf,
                        n_lining_residues = n_lining)
  list(file = file, row = row, pockets = pockets)
}


## summarise all structures in parallel. Any structure that fails stops the script, so a table is never
## made with structures silently missing.
files <- list.files(scores_dir, pattern = "_info\\.txt$")
res <- mclapply(files, function(f) {
  tryCatch(summarise_file(f), error = function(e) paste("Could not summarise", f, ":", conditionMessage(e)))
}, mc.cores = threads)
failed <- vapply(res, is.character, logical(1))
if (any(failed)) stop(sum(failed), " fpocket outputs could not be summarised:\n", paste(head(unlist(res[failed]), 5), collapse = "\n"))
res <- res[!vapply(res, is.null, logical(1))]   # structures where fpocket found no pockets
results <- do.call(rbind, lapply(res, function(x) x$row))
pockets <- do.call(rbind, lapply(res, function(x) x$pockets))
cat(sprintf("%d of %d fpocket outputs summarised (%d with no pockets), %d pockets\n",
            nrow(results), length(files), length(files) - nrow(results), nrow(pockets)))


## add structure confidence from the filter step (SWISS-MODEL x 100), joined by structure file name
struct_conf <- read.csv(struct_score_file)
struct_conf$struct_score <- struct_conf$struct_score * conf_scale
results$filename <- paste0(sub("_info\\.txt$", "", vapply(res, function(x) x$file, character(1))), ".pdb")
results_master <- merge(results, struct_conf[c("filename", "struct_score")], by = "filename")
results_master <- results_master[, c("uniprot_id", "structure_id", "n_residues", "struct_score", "pocket", "druggability",
                                     "num_drug_pockets", "hydrophobicity_score", "volume_score", "alpha_max_dist",
                                     "conf_best_pocket", "druggability_conf", "num_drug_pockets_conf")]

# Order the table by highest druggability scores, and give every column but uniprot_id the source prefix
results_master <- results_master[order(-results_master$druggability), ]
rownames(results_master) <- NULL
names(results_master)[-1] <- paste0(prefix, names(results_master)[-1])
names(pockets)[-1] <- paste0(prefix, names(pockets)[-1])

## main table: structures of at least min_residues residues; the shorter ones go to small_structures.csv
too_short <- results_master[[paste0(prefix, "n_residues")]] < min_residues
write_csv(results_master[!too_short, ], file.path(results_dir, "fpocket_druggability_full.csv"))
write_csv(results_master[too_short, ], file.path(results_dir, "small_structures.csv"))
write_csv(pockets, file.path(results_dir, "fpocket_pockets.csv"))
cat(sprintf("%d structures in fpocket_druggability_full.csv, %d under %d residues in small_structures.csv\n",
            sum(!too_short), sum(too_short), min_residues))
