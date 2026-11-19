### scores structure confidence and builds the filtered structure set that fpocket and PocketMiner score ###
# usage: Rscript struct_conf.R <alphafold|swissmodel> <structures_dir> <filtered_dir> <min_conf> <threads>
#
# Each structure's confidence is the mean B-factor of its atoms (pLDDT for AlphaFold, QMEAN local for SWISS-MODEL).
# SWISS-MODEL models are trimmed (protein ATOM records only, no insertion codes, first chain only) and written out
# without header records; AlphaFold models go in as they are, as links to the raw files.
# Raw input files are never moved or changed.

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
source_type <- args[1]
pdb_dir     <- args[2]
filtered_dir <- args[3]
min_conf    <- as.numeric(args[4])
threads     <- as.integer(args[5])

# fresh pdb/ dir each run (it is rebuilt from the raw files, so nothing is lost)
input_dir <- file.path(filtered_dir, "pdb")
unlink(input_dir, recursive = TRUE)
dir.create(input_dir, recursive = TRUE)

pdb_files <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE)


## score one structure; if it passes, add it to pdb/
score_structure <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)
  pdb_df <- as.data.frame(pdb$atom)
  filename <- basename(pdb_file)

  # SWISS-MODEL trimming: remove non standard residues and inserts, keep first chain
  if (source_type == "swissmodel") {
    pdb_df <- pdb_df[pdb_df$type == "ATOM", ]
    pdb_df <- pdb_df[is.na(pdb_df$insert), ]
    if (length(unique(pdb_df$chain)) > 1) {
      first_chain <- pdb_df$chain[1]
      pdb_df <- pdb_df[pdb_df$chain == first_chain, ]
    }
    uniprot_id <- sub("\\.pdb$", "", filename)
    file_id <- pdb_file
  } else {
    uniprot_id <- sub("^[^-]+-([^-]+)-.*", "\\1", filename)
    file_id <- sub("\\.pdb$", "", filename)
  }

  # check all atoms in each residue share one confidence value
  resno_b_df <- subset(pdb_df, select = c("resno", "b"))
  unique_resno <- unique(resno_b_df$resno)
  b_values <- sapply(unique_resno, function(resno) unique(resno_b_df[resno_b_df$resno == resno, "b"]))
  atom_residue <- all(sapply(b_values, function(b) length(b) == 1))

  # structure confidence = mean B-factor (pLDDT for AlphaFold, QMEAN local for SWISS-MODEL)
  struct_score <- mean(pdb_df$b)

  # passing structures: headerless trimmed copy (SWISS-MODEL) or link to the raw file (AlphaFold)
  if (struct_score >= min_conf) {
    if (source_type == "swissmodel") {
      pdb$atom <- pdb_df
      xyz <- as.numeric(c(t(pdb_df[, c("x", "y", "z")])))
      dim(xyz) <- c(1, length(xyz))
      pdb$xyz <- xyz
      pdb$calpha <- pdb_df$elety == "CA"
      write.pdb(pdb, file = file.path(input_dir, paste0(uniprot_id, ".pdb")))
    } else {
      file.symlink(normalizePath(pdb_file), file.path(input_dir, filename))
    }
  }

  data.frame(filename = filename, file_id = file_id, uniprot_id = uniprot_id,
             atom_residue = atom_residue, struct_score = struct_score)
}


## score all structures in parallel; unreadable files are listed rather than stopping the run
res <- mclapply(pdb_files, function(f) {
  tryCatch(score_structure(f), error = function(e) {
    message("Could not score ", f, ": ", conditionMessage(e))
    NULL
  })
}, mc.cores = threads)

failed <- pdb_files[vapply(res, is.null, logical(1))]
scores <- do.call(rbind, res)
result_df  <- scores[scores$struct_score >= min_conf, ]
removed_df <- scores[scores$struct_score <  min_conf, ]


## warn about mixed per-atom confidence, then write tables
if (any(result_df$atom_residue == FALSE)) {
  message("Warning!: Some atoms within a single residue have varying confidence scores.
          Check struct_score.csv for these structures (look for the FALSE values).")
}

write.csv(result_df,  file.path(filtered_dir, "struct_score.csv"), row.names = F)
write.csv(removed_df, file.path(filtered_dir, "low_conf_struct.csv"), row.names = F)
writeLines(failed,    file.path(filtered_dir, "struct_failed.txt"))

cat(sprintf("%d structures passed (struct_score >= %s), %d below, %d unreadable\n",
            nrow(result_df), min_conf, nrow(removed_df), length(failed)))
