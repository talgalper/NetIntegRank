### Combine fpocket and PocketMiner into one druggability score per protein ###
#
# usage: Rscript combine_scores.R --af-fpocket <csv[,csv]> --sm-fpocket <csv[,csv]>
#                                 --af-pocketminer <csv[,csv]> --sm-pocketminer <csv[,csv]> --out <combined.csv>
#                                 [--fpocket-score druggability|druggability_conf] [--pocketminer-score max_hit|max_hit_conf]
#                                 [--min-residues <n>]
#   Several tables per input are given comma-separated (e.g. a main table and its small_structures.csv).
#   Reads the 2024 tables (plain column names) and the 2026 tables (AF_ / SM_ column names) alike.
#   --min-residues: structures shorter than this are not used for scores (2026 tables only; they have residue counts).
#   Every protein in the input tables is kept: a protein with no usable score from either tool gets highest_score 0
#   (Tal, 2026-11-18), because NetIntegRank leaves genes without a drug_score out of its ranking altogether.
#
# The rule, from PhD-MOC Druggability_analysis/R/plot_drug_freq.R (2024), lines 1-42:
#   - fpocket score per protein: the highest druggability over its AlphaFold fragments and its SWISS-MODEL model
#   - PocketMiner score (CP_score): the highest max_hit over the same
#   - SWISS-MODEL struct_score on the 0-100 scale; missing scores are 0
#   - highest_score = max(fpocket, PocketMiner); source = AF or SM, the structure that gave highest_score
#     (fpocket wins a tie)
#
# Changes from 2024 (2026-11-18; highest_score is unchanged by them):
#   - Each tool's best score is found first, so druggability / Fpocket_method always hold the protein's best fpocket
#     score and CP_score / PocketMiner_method its best PocketMiner score.
#   - A tool column says which tool gave highest_score.
#   - A structure with no score (no confident pocket or window) is not used for the best score; the protein stays in
#     the table with 0 for that tool.

## arguments: --name value pairs
args <- commandArgs(trailingOnly = TRUE)
opt <- list(fpocket_score = "druggability", pocketminer_score = "max_hit")   # defaults: scores as run
for (i in seq(1, length(args), by = 2)) {
  opt[[gsub("-", "_", sub("^--", "", args[i]))]] <- args[i + 1]
}
min_residues <- if (is.null(opt$min_residues)) NA else as.integer(opt$min_residues)
for (k in c("af_fpocket", "sm_fpocket", "af_pocketminer", "sm_pocketminer", "out")) {
  if (is.null(opt[[k]])) stop("missing --", gsub("_", "-", k), "; see the usage at the top of combine_scores.R")
}

## read one input: one or more tables, stacked
read_tables <- function(files) {
  do.call(rbind, lapply(strsplit(files, ",")[[1]], read.csv))
}

## fpocket rows: uniprot_id, score, struct_score (0-100), method
fpocket_rows <- function(files, source) {
  x <- read_tables(files)
  prefix <- paste0(source, "_")
  if (paste0(prefix, "druggability") %in% names(x)) {          # 2026 tables: SWISS-MODEL struct_score already x 100
    score <- x[[paste0(prefix, opt$fpocket_score)]]
    struct_score <- x[[paste0(prefix, "struct_score")]]
    if (!is.na(min_residues)) score[x[[paste0(prefix, "n_residues")]] < min_residues] <- NA   # too short: no score
  } else {                                                     # 2024 tables
    if (opt$fpocket_score != "druggability") stop("the 2024 fpocket tables only have druggability")
    if (!is.na(min_residues)) stop("the 2024 fpocket tables have no residue counts, so --min-residues cannot be used")
    score <- x$druggability
    struct_score <- if (source == "SM") x$struct_score * 100 else x$struct_score
  }
  if (is.null(score)) stop("no ", opt$fpocket_score, " column in ", files)
  data.frame(uniprot_id = x$uniprot_id, score = score, struct_score = struct_score, method = source)
}

## PocketMiner rows: uniprot_id, score, method
pocketminer_rows <- function(files, source) {
  x <- read_tables(files)
  prefix <- paste0(source, "_")
  if (paste0(prefix, "max_hit") %in% names(x)) {               # 2026 tables
    score <- x[[paste0(prefix, opt$pocketminer_score)]]
    if (!is.na(min_residues)) score[x[[paste0(prefix, "n_residues")]] < min_residues] <- NA   # too short: no score
  } else {                                                     # 2024 tables
    if (opt$pocketminer_score != "max_hit") stop("the 2024 PocketMiner tables only have max_hit")
    if (!is.na(min_residues)) stop("the 2024 PocketMiner tables have no residue counts, so --min-residues cannot be used")
    score <- x$max_hit
  }
  if (is.null(score)) stop("no ", opt$pocketminer_score, " column in ", files)
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

# the higher of the two, and where it came from (fpocket wins a tie, as in 2024)
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
