### Combine fpocket and PocketMiner into one druggability score per protein ###
#
# usage: Rscript combine_scores.R --af-fpocket <csv[,csv]> --sm-fpocket <csv[,csv]>
#                                 --af-pocketminer <csv[,csv]> --sm-pocketminer <csv[,csv]> --out <combined.csv>
#                                 [--fpocket-score druggability_conf|druggability] [--pocketminer-score max_hit_conf|max_hit]
#                                 [--min-residues <n>]
#   Inputs are the tables written by Fpocket_master.sh and pocketminer_master.sh (columns starting AF_ or SM_).
#   Several tables per input are given comma-separated (e.g. a main table and its small_structures.csv).
#   --fpocket-score      fpocket column to use: druggability_conf (default, confident pockets only) or druggability
#                        (all pockets)
#   --pocketminer-score  PocketMiner column to use: max_hit_conf (default, confident windows only) or max_hit (all
#                        windows)
#   --min-residues       structures shorter than this are not used for scores (default 80; 0 uses every structure)
#   The defaults are the settings of the bundled assets/druggability_scores_annot.csv.
#
# The rule:
#   - fpocket score per protein (druggability): the highest score over its AlphaFold fragments and its SWISS-MODEL
#     model; Fpocket_method says which source gave it, struct_score is that structure's whole-model confidence
#     (0-100)
#   - PocketMiner score per protein (CP_score): the highest score over the same; PocketMiner_method says which source
#   - a structure with no score (too short, or no confident pocket or window) is not used; a protein with no usable
#     score from a tool gets 0 for that tool, so every protein in the inputs is kept (NetIntegRank leaves a gene with
#     no drug_score out of its ranking altogether)
#   - highest_score = the higher of the two; tool and source say which tool and which structure gave it (fpocket wins
#     a tie between tools, AlphaFold a tie between sources)

## arguments: --name value pairs
args <- commandArgs(trailingOnly = TRUE)
opt <- list(fpocket_score = "druggability_conf", pocketminer_score = "max_hit_conf", min_residues = "80")
for (i in seq(1, length(args), by = 2)) {
  opt[[gsub("-", "_", sub("^--", "", args[i]))]] <- args[i + 1]
}
min_residues <- as.integer(opt$min_residues)
if (is.na(min_residues)) stop("--min-residues must be a whole number")
for (k in c("af_fpocket", "sm_fpocket", "af_pocketminer", "sm_pocketminer", "out")) {
  if (is.null(opt[[k]])) stop("missing --", gsub("_", "-", k), "; see the usage at the top of combine_scores.R")
}

## read one input: one or more tables, stacked
read_tables <- function(files) {
  do.call(rbind, lapply(strsplit(files, ",")[[1]], read.csv))
}

## check a table has the columns this script reads (tables from Fpocket_master.sh / pocketminer_master.sh do)
require_columns <- function(x, columns, files) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0) stop("no ", paste(missing, collapse = ", "), " column in ", files,
                                "; give the tables written by Fpocket_master.sh and pocketminer_master.sh")
}

## fpocket rows: uniprot_id, score, struct_score (0-100; SWISS-MODEL already x 100), method
fpocket_rows <- function(files, source) {
  x <- read_tables(files)
  prefix <- paste0(source, "_")
  columns <- paste0(prefix, c(opt$fpocket_score, "struct_score", "n_residues"))
  require_columns(x, c("uniprot_id", columns), files)
  score <- x[[columns[1]]]
  score[x[[columns[3]]] < min_residues] <- NA                  # too short: no score
  data.frame(uniprot_id = x$uniprot_id, score = score, struct_score = x[[columns[2]]], method = source)
}

## PocketMiner rows: uniprot_id, score, method
pocketminer_rows <- function(files, source) {
  x <- read_tables(files)
  prefix <- paste0(source, "_")
  columns <- paste0(prefix, c(opt$pocketminer_score, "n_residues"))
  require_columns(x, c("uniprot_id", columns), files)
  score <- x[[columns[1]]]
  score[x[[columns[2]]] < min_residues] <- NA                  # too short: no score
  data.frame(uniprot_id = x$uniprot_id, score = score, method = source)
}

fp <- rbind(fpocket_rows(opt$af_fpocket, "AF"), fpocket_rows(opt$sm_fpocket, "SM"))
pm <- rbind(pocketminer_rows(opt$af_pocketminer, "AF"), pocketminer_rows(opt$sm_pocketminer, "SM"))


## best score per protein for each tool. Structures with no score are skipped. The sort keeps the input order for
## ties, and AlphaFold rows come first, so AlphaFold wins a tie between sources.
fp_best <- fp[!is.na(fp$score), ]
fp_best <- fp_best[order(-fp_best$score), ]
fp_best <- fp_best[!duplicated(fp_best$uniprot_id), ]
colnames(fp_best) <- c("uniprot_id", "druggability", "struct_score", "Fpocket_method")

pm_best <- pm[!is.na(pm$score), ]
pm_best <- pm_best[order(-pm_best$score), ]
pm_best <- pm_best[!duplicated(pm_best$uniprot_id), ]
colnames(pm_best) <- c("uniprot_id", "CP_score", "PocketMiner_method")


## one row per protein: every protein in any input table
combined <- data.frame(uniprot_id = unique(c(fp$uniprot_id, pm$uniprot_id)))
combined <- merge(combined, fp_best, by = "uniprot_id", all.x = TRUE)
combined <- merge(combined, pm_best, by = "uniprot_id", all.x = TRUE)
combined$druggability[is.na(combined$druggability)] <- 0
combined$CP_score[is.na(combined$CP_score)] <- 0

# the higher of the two, and where it came from (fpocket wins a tie)
combined$highest_score <- pmax(combined$druggability, combined$CP_score)
fpocket_wins <- combined$highest_score == combined$druggability
combined$source <- ifelse(fpocket_wins, combined$Fpocket_method, combined$PocketMiner_method)
combined$tool <- ifelse(fpocket_wins, "fpocket", "PocketMiner")
combined$tool[is.na(combined$source)] <- NA          # no score from either tool

combined <- combined[order(-combined$highest_score), c("uniprot_id", "druggability", "struct_score", "Fpocket_method",
                                                        "CP_score", "PocketMiner_method", "highest_score", "source", "tool")]
write.csv(combined, opt$out, row.names = FALSE)

## summary
cat(sprintf("%d proteins (fpocket score for %d, PocketMiner score for %d, neither %d); PocketMiner sets highest_score for %d (%.1f%%)\n",
            nrow(combined), sum(!is.na(combined$Fpocket_method)), sum(!is.na(combined$PocketMiner_method)),
            sum(is.na(combined$tool)), sum(combined$tool == "PocketMiner", na.rm = TRUE),
            100 * mean(combined$tool == "PocketMiner", na.rm = TRUE)))
cat("written to", opt$out, "\n")
