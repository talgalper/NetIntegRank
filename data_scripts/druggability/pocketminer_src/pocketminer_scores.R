#!/usr/bin/env Rscript

### PocketMiner scoring: turns the per-residue predictions into one row per structure ###
### (called by pocketminer_master.sh)
#
# usage: Rscript pocketminer_scores.R <results_dir> <threads> <min_residues> <window> <hit_threshold>
#                                     <alphafold|swissmodel> <structures_dir> <min_conf>
#
# Outputs, in results_dir (columns start with AF_ or SM_, except uniprot_id, so the two sources merge easily):
#   pocketminer_results.csv  structures with at least min_residues residues
#   small_structures.csv     the shorter structures, scored the same way and left out of the main table
#   pocket_results/          per residue: raw score, averaged score, confidence, averaged confidence
#
# Confidence is read from the B-factor column of each residue's CA atom: AlphaFold pLDDT (0-100) or the SWISS-MODEL
# QMEANDisCo local score (0-1, multiplied by 100 here so both sources use one scale).

suppressPackageStartupMessages(library(parallel))

## arguments: where the results are, how many CPUs to use, the scoring settings, and where the structures are
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: Rscript pocketminer_scores.R <results_dir> <threads> <min_residues> <window> <hit_threshold> ",
       "<alphafold|swissmodel> <structures_dir> <min_conf>")
}
results_dir <- args[1]
threads <- as.integer(args[2])
min_residues <- as.integer(args[3])     # structures shorter than this are left out of the main table
window_size <- as.integer(args[4])      # residues averaged on either side of each residue
hit_threshold <- as.numeric(args[5])    # an averaged score at or above this is a hit
source_type <- args[6]                  # alphafold or swissmodel
structures_dir <- args[7]               # the structures PocketMiner was run on
min_conf <- as.numeric(args[8])         # a window is confident if its mean confidence (0-100 scale) is at least this
if (!source_type %in% c("alphafold", "swissmodel")) stop("source must be alphafold or swissmodel")
prefix <- if (source_type == "alphafold") "AF_" else "SM_"
conf_scale <- if (source_type == "alphafold") 1 else 100   # QMEANDisCo 0-1 -> 0-100

## the per-residue prediction files to score, and the folder for the averaged scores
txt_dir <- file.path(results_dir, "txt_results")
pocket_dir <- file.path(results_dir, "pocket_results")
dir.create(pocket_dir, showWarnings = FALSE)
files <- list.files(txt_dir, pattern = "-predictions\\.txt$")
if (length(files) == 0) stop("no prediction files in ", txt_dir)

## which structure file each prediction came from (same naming rule as pocketminer_predict.py)
pdb_files <- list.files(structures_dir, pattern = "\\.pdb$")
pdb_ids <- if (source_type == "alphafold") sub("^AF-([^-]+)-([^-]+)-.*$", "\\1-\\2", pdb_files) else sub("\\.pdb$", "", pdb_files)
ids <- sub("-predictions\\.txt$", "", files)
pdb_of <- pdb_files[match(ids, pdb_ids)]                    # structure file for each prediction file, in order
if (anyNA(pdb_of)) stop(sum(is.na(pdb_of)), " predictions have no structure in ", structures_dir,
                        ", starting with ", ids[is.na(pdb_of)][1])

## average each residue with the residues within the window on either side (fewer at the two ends of the chain)
window_mean <- function(x) {
  n <- length(x)
  vapply(seq_len(n), function(i) mean(x[max(1, i - window_size):min(n, i + window_size)]), numeric(1))
}

## score one structure
score_one <- function(i) {
  ## read the scores, one per residue. The protein ID is the part of the name before the first "-".
  file <- files[i]
  id <- ids[i]
  scores <- scan(file.path(txt_dir, file), quiet = TRUE)
  n_residues <- length(scores)
  uniprot_id <- strsplit(id, "-")[[1]][1]

  ## per-residue confidence: B-factor of each residue's CA atom, which must line up one to one with the scores
  lines <- readLines(file.path(structures_dir, pdb_of[i]))
  ca <- lines[startsWith(lines, "ATOM") & substr(lines, 13, 16) == " CA "]
  confidence <- round(as.numeric(substr(ca, 61, 66)) * conf_scale, 6)   # rounded so 0.29 x 100 is exactly 29
  if (length(confidence) != n_residues) {
    stop(id, ": ", n_residues, " predictions but ", length(confidence), " CA atoms in ", pdb_of[i])
  }

  ## smooth the scores: each residue becomes the average of itself and the residues within the window
  ## on either side. A cryptic pocket is a patch of neighbouring residues, so a single high residue on
  ## its own counts for less than a run of them. The confidence is averaged over the same window.
  averages <- window_mean(scores)
  confidence_avg <- window_mean(confidence)

  ## save the raw and averaged values per residue, then read them back so the summary below is
  ## worked out from exactly the numbers that are on disk
  out_file <- file.path(pocket_dir, file)
  write.table(data.frame(scores = scores, averages = averages, confidence = confidence, confidence_avg = confidence_avg),
              out_file, quote = FALSE, row.names = FALSE)
  saved <- read.table(out_file, header = TRUE)
  averages <- saved$averages
  confidence_avg <- saved$confidence_avg

  ## all windows:
  ## max_hit: the highest averaged score in the structure
  ## num_hits: how many separate stretches of the chain sit at or above the threshold
  stretches <- rle(averages >= hit_threshold)
  max_hit <- max(averages)
  num_hits <- sum(stretches$values)

  ## confident regions only (a window is confident if its mean confidence is at least min_conf):
  ## max_hit_conf: the highest averaged score among confident windows (NA if the structure has none)
  ## num_hits_conf: hit stretches that reach the threshold inside a confident window
  ## conf_at_max_hit: the mean confidence of the max_hit window
  ## frac_low_conf: share of residues whose own confidence is below min_conf
  confident <- confidence_avg >= min_conf
  max_hit_conf <- if (any(confident)) max(averages[confident]) else NA_real_
  ends <- cumsum(stretches$lengths)
  starts <- ends - stretches$lengths + 1
  num_hits_conf <- 0L
  for (k in which(stretches$values)) {
    if (any(confident[starts[k]:ends[k]])) num_hits_conf <- num_hits_conf + 1L
  }

  list(uniprot_id = uniprot_id, structure_id = id, n_residues = n_residues,
       frac_low_conf = mean(confidence < min_conf), max_hit = max_hit, num_hits = num_hits,
       conf_at_max_hit = confidence_avg[which.max(averages)], max_hit_conf = max_hit_conf, num_hits_conf = num_hits_conf)
}

## score all structures in parallel, a block at a time so progress can be shown
cat("  scoring", length(files), "structures\n")
blocks <- split(seq_along(files), ceiling(seq_along(files) / 500))
bar <- txtProgressBar(max = length(files), style = if (isatty(stdout())) 3 else 1)
rows <- list()
for (block in blocks) {
  rows <- c(rows, mclapply(block, score_one, mc.cores = threads))
  setTxtProgressBar(bar, length(rows))
}
close(bar)

## stop if any structure could not be scored, naming the first few and why
not_scored <- !vapply(rows, function(row) is.list(row) && !is.null(row$structure_id), logical(1))
if (any(not_scored)) {
  why <- vapply(rows[not_scored], function(row) as.character(row)[1], character(1))
  stop("scoring failed for ", sum(not_scored), " files, starting with: ", paste(head(why, 3), collapse = " | "))
}

## combine into one table: one row per structure
pocketminer_master <- data.frame(
  uniprot_id = vapply(rows, function(row) row$uniprot_id, character(1)),
  structure_id = vapply(rows, function(row) row$structure_id, character(1)),
  n_residues = vapply(rows, function(row) row$n_residues, integer(1)),
  frac_low_conf = vapply(rows, function(row) row$frac_low_conf, numeric(1)),
  max_hit = vapply(rows, function(row) row$max_hit, numeric(1)),
  num_hits = vapply(rows, function(row) row$num_hits, integer(1)),
  conf_at_max_hit = vapply(rows, function(row) row$conf_at_max_hit, numeric(1)),
  max_hit_conf = vapply(rows, function(row) row$max_hit_conf, numeric(1)),
  num_hits_conf = vapply(rows, function(row) row$num_hits_conf, integer(1)))
names(pocketminer_master)[-1] <- paste0(prefix, names(pocketminer_master)[-1])

## save the results: structures of at least min_residues residues in the main table, the shorter ones in
## small_structures.csv with the same columns. Nothing is moved, the prediction files all stay in txt_results/.
too_short <- pocketminer_master[[paste0(prefix, "n_residues")]] < min_residues
write.csv(pocketminer_master[!too_short, ], file.path(results_dir, "pocketminer_results.csv"), row.names = FALSE)
write.csv(pocketminer_master[too_short, ], file.path(results_dir, "small_structures.csv"), row.names = FALSE)
cat("  ", sum(!too_short), " structures scored, ", sum(too_short), " under ", min_residues,
    " residues written to small_structures.csv\n", sep = "")
