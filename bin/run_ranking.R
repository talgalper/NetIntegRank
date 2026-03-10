#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--hhnet_metrics", type = "character",
              help = "Path to HHNet-derived network metrics table (TSV/CSV)."),
  make_option("--druggability", type = "character",
              help = "Path to druggability table (TSV/CSV). Required columns: uniprot_gn_id, highest_score."),
  make_option("--ml_scores", type = "character",
              help = "Path to ML scores table (TSV/CSV). Required columns: Protein, Prediction_Score_rf."),
  make_option("--citations", type = "character",
              help = "Path to citation counts table (TSV/CSV). Required columns: symbol, counts."),
  make_option("--gene_map", type = "character", default = NULL,
              help = "Optional mapping table (TSV/CSV). Recommended columns: ensembl_gene_id, uniprot_gn_id (external_gene_name optional)."),
  make_option("--ranking_features", type = "character",
              default = "degree,betweenness,closeness,eigen_centrality,page_rank,highest_score",
              help = "Comma-separated numeric feature columns to include in avg_rank. Default: degree,betweenness,closeness,eigen_centrality,page_rank,highest_score"),
  make_option("--out_tsv", type = "character",
              help = "Output TSV path for final ranking."),
  make_option("--out_rds", type = "character",
              help = "Output RDS path for final ranking."),
  make_option("--id_annot_cache", type = "character", default = NULL,
              help = "Optional cache path for biomaRt ID annotations. Reused across runs when provided."),
  make_option("--out_incomplete_tsv", type = "character", default = NULL,
              help = "Optional output TSV path for rows removed due to incomplete post-merge data."),
  make_option("--out_incomplete_rds", type = "character", default = NULL,
              help = "Optional output RDS path for rows removed due to incomplete post-merge data.")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

required_args <- c("hhnet_metrics", "druggability", "ml_scores", "citations", "out_tsv", "out_rds")
missing_args <- required_args[vapply(required_args, function(x) is.null(args[[x]]) || !nzchar(args[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(sprintf("--%s", missing_args), collapse = ", ")), call. = FALSE)
}

read_table_auto <- function(path) {
  if (!file.exists(path)) stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.table(path, header = TRUE, sep = "\t", quote = "", check.names = FALSE,
               stringsAsFactors = FALSE, comment.char = "")
  }
}

require_columns <- function(df, required_cols, table_label, table_path) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Input validation failed for %s (%s). Missing required columns: %s",
                 table_label, table_path, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
}

require_numeric <- function(df, cols, table_label, table_path) {
  bad_cols <- cols[!vapply(df[cols], is.numeric, logical(1))]
  if (length(bad_cols) > 0) {
    stop(sprintf("Input validation failed for %s (%s). Non-numeric columns where numeric expected: %s",
                 table_label, table_path, paste(bad_cols, collapse = ", ")), call. = FALSE)
  }
}

normalize_id_column <- function(df, colname) {
  if (!colname %in% colnames(df)) return(df)
  df[[colname]] <- trimws(as.character(df[[colname]]))
  df[[colname]][df[[colname]] == ""] <- NA_character_
  df
}

require_non_missing <- function(df, cols, table_label, table_path) {
  missing_counts <- vapply(cols, function(col) sum(is.na(df[[col]]) | trimws(as.character(df[[col]])) == ""), integer(1))
  bad_cols <- names(missing_counts)[missing_counts > 0]
  if (length(bad_cols) > 0) {
    bad_msg <- paste(sprintf("%s (%d missing)", bad_cols, missing_counts[bad_cols]), collapse = ", ")
    stop(sprintf("Input validation failed for %s (%s). Required identifier/value columns contain missing values: %s",
                 table_label, table_path, bad_msg), call. = FALSE)
  }
}

upsert_annotation_cache <- function(cache_path, input_type, convert_to, annotation_df) {
  cache_cols <- c(input_type, convert_to)
  cache_df <- annotation_df[, cache_cols, drop = FALSE]

  if (file.exists(cache_path)) {
    existing <- read_table_auto(cache_path)
    require_columns(existing, cache_cols, "id_annot_cache", cache_path)
    cache_df <- rbind(existing[, cache_cols, drop = FALSE], cache_df)
  }

  cache_df <- cache_df[!duplicated(cache_df[[input_type]]), , drop = FALSE]
  outdir <- dirname(cache_path)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.table(cache_df, file = cache_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cache_df
}

compute_avg_rank <- function(df, feature_cols) {
  ranked_features <- lapply(feature_cols, function(colname) {
    rank(-df[[colname]], ties.method = "average", na.last = "keep")
  })
  rank_matrix <- do.call(cbind, ranked_features)
  colnames(rank_matrix) <- paste0(feature_cols, "_rank")

  avg <- rowMeans(rank_matrix, na.rm = TRUE)
  avg[is.nan(avg)] <- Inf  # rows with all NA features go to bottom
  df$avg_rank <- avg
  cbind(df, as.data.frame(rank_matrix, check.names = FALSE))
}

# Original id_annot
id_annot <- function(data, input_type, convert_to, cache_path = NULL) {
  if (!input_type %in% colnames(data)) {
    stop(sprintf("id_annot input column '%s' not found in data.", input_type), call. = FALSE)
  }

  if (!all(convert_to %in% c("uniprot_gn_id", "external_gene_name"))) {
    stop("id_annot currently supports convert_to values: 'uniprot_gn_id', 'external_gene_name'.", call. = FALSE)
  }

  keys <- unique(as.character(data[[input_type]]))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  cache_cols <- c(input_type, convert_to)
  cached <- NULL

  if (!is.null(cache_path) && nzchar(cache_path) && file.exists(cache_path)) {
    cached <- read_table_auto(cache_path)
    require_columns(cached, cache_cols, "id_annot_cache", cache_path)
    cached <- cached[, cache_cols, drop = FALSE]
    cached <- cached[!duplicated(cached[[input_type]]), , drop = FALSE]
  }

  missing_keys <- keys
  if (!is.null(cached)) missing_keys <- setdiff(keys, as.character(cached[[input_type]]))

  fetched <- data.frame(setNames(vector("list", length(cache_cols)), cache_cols), stringsAsFactors = FALSE)

  if (length(missing_keys) > 0) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("Package 'biomaRt' is required when --gene_map is not supplied and cache is incomplete.", call. = FALSE)
    }

    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attrs <- unique(c(input_type, "uniprotswissprot", "uniprot_gn_id", convert_to))

    bm <- biomaRt::getBM(
      attributes = attrs,
      filters = input_type,
      values = missing_keys,
      mart = ensembl
    )

    by_key <- split(bm, bm[[input_type]])
    out_list <- lapply(by_key, function(df_sub) {
      swiss <- unique(stats::na.omit(ifelse(df_sub$uniprotswissprot == "", NA, df_sub$uniprotswissprot)))
      gids  <- unique(stats::na.omit(ifelse(df_sub$uniprot_gn_id == "", NA, df_sub$uniprot_gn_id)))
      selected_uniprot <- if (length(swiss) > 0) swiss else gids

      row <- setNames(as.list(rep(NA_character_, 1 + length(convert_to))), c(input_type, convert_to))
      row[[input_type]] <- df_sub[[input_type]][1]
      if ("uniprot_gn_id" %in% convert_to) {
        row[["uniprot_gn_id"]] <- if (length(selected_uniprot) > 0) paste(selected_uniprot, collapse = ";") else NA_character_
      }
      if ("external_gene_name" %in% convert_to) {
        names_vals <- unique(stats::na.omit(ifelse(df_sub$external_gene_name == "", NA, df_sub$external_gene_name)))
        row[["external_gene_name"]] <- if (length(names_vals) > 0) paste(names_vals, collapse = ";") else NA_character_
      }
      as.data.frame(row, stringsAsFactors = FALSE)
    })

    if (length(out_list) > 0) fetched <- do.call(rbind, out_list)

    if (!is.null(cache_path) && nzchar(cache_path) && nrow(fetched) > 0) {
      cached <- upsert_annotation_cache(cache_path, input_type, convert_to, fetched)
    }
  }

  annot <- if (!is.null(cached)) cached else fetched

  if (nrow(annot) == 0) {
    for (colname in convert_to) data[[colname]] <- NA_character_
    return(data)
  }

  annot <- annot[!duplicated(annot[[input_type]]), c(input_type, convert_to), drop = FALSE]
  merge(data, annot, by = input_type, all.x = TRUE, sort = FALSE)
}

primary_id <- function(x) {
  if (is.null(x)) return(character(0))

  x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- trimws(x)
  x[!nzchar(x)] <- NA_character_

  split_x <- strsplit(replace(x, is.na(x), ""), ";", fixed = TRUE)

  vapply(split_x, function(parts) {
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) NA_character_ else parts[1]
  }, character(1))
}

flag_incomplete_rows <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Cannot assess incomplete rows. Missing columns in merged table: %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  missing_matrix <- sapply(required_cols, function(col) {
    x <- df[[col]]
    if (is.character(x)) {
      is.na(x) | trimws(x) == ""
    } else {
      is.na(x)
    }
  }, simplify = "matrix")

  if (length(required_cols) == 1) {
    missing_matrix <- matrix(missing_matrix, ncol = 1)
    colnames(missing_matrix) <- required_cols
  }

  incomplete_idx <- apply(missing_matrix, 1, any)
  missing_fields <- apply(missing_matrix, 1, function(x) {
    paste(required_cols[which(x)], collapse = ";")
  })

  incomplete_rows <- df[incomplete_idx, , drop = FALSE]
  if (nrow(incomplete_rows) > 0) {
    incomplete_rows$missing_fields <- missing_fields[incomplete_idx]
  }

  complete_rows <- df[!incomplete_idx, , drop = FALSE]

  list(
    complete = complete_rows,
    incomplete = incomplete_rows
  )
}

# ------------------------
# Load inputs
# ------------------------

hhnet_metrics <- read_table_auto(args$hhnet_metrics)
druggability  <- read_table_auto(args$druggability)
ml_scores     <- read_table_auto(args$ml_scores)
citations     <- read_table_auto(args$citations)

# Expected HHNet-post schema (minimum required for downstream joins + defaults)
required_hh_cols <- c(
  "ensembl_gene_id","external_gene_name","description","gene_biotype",
  "degree","betweenness","closeness","eigen_centrality","page_rank",
  "cluster","source"
)
require_columns(hhnet_metrics, required_hh_cols, "hhnet_metrics", args$hhnet_metrics)
require_columns(druggability, c("uniprot_gn_id", "highest_score"), "druggability", args$druggability)
require_columns(ml_scores, c("Protein", "Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_columns(citations, c("symbol", "counts"), "citations", args$citations)

require_numeric(hhnet_metrics, c("degree","betweenness","closeness","eigen_centrality","page_rank"), "hhnet_metrics", args$hhnet_metrics)
require_numeric(druggability, c("highest_score"), "druggability", args$druggability)
require_numeric(ml_scores, c("Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_numeric(citations, c("counts"), "citations", args$citations)

# Normalise IDs
hhnet_metrics <- normalize_id_column(hhnet_metrics, "ensembl_gene_id")
hhnet_metrics <- normalize_id_column(hhnet_metrics, "external_gene_name")
druggability  <- normalize_id_column(druggability, "uniprot_gn_id")
ml_scores     <- normalize_id_column(ml_scores, "uniprot_gn_id")
citations     <- normalize_id_column(citations, "external_gene_name")

require_non_missing(hhnet_metrics, c("ensembl_gene_id"), "hhnet_metrics", args$hhnet_metrics)
require_non_missing(druggability, c("uniprot_gn_id","highest_score"), "druggability", args$druggability)
require_non_missing(ml_scores, c("Protein","Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_non_missing(citations, c("symbol","counts"), "citations", args$citations)

# Derive default filenames if the user does not supply them for incomplete cases file
if (is.null(args$out_incomplete_tsv) || !nzchar(args$out_incomplete_tsv)) {
  args$out_incomplete_tsv <- sub("\\.tsv$", "_incomplete.tsv", args$out_tsv, ignore.case = TRUE)
  if (identical(args$out_incomplete_tsv, args$out_tsv)) {
    args$out_incomplete_tsv <- paste0(args$out_tsv, "_incomplete.tsv")
  }
}
if (is.null(args$out_incomplete_rds) || !nzchar(args$out_incomplete_rds)) {
  args$out_incomplete_rds <- sub("\\.rds$", "_incomplete.rds", args$out_rds, ignore.case = TRUE)
  if (identical(args$out_incomplete_rds, args$out_rds)) {
    args$out_incomplete_rds <- paste0(args$out_rds, "_incomplete.rds")
  }
}

id_annot_cache_path <- args$id_annot_cache
if (is.null(id_annot_cache_path) || !nzchar(id_annot_cache_path)) {
  id_annot_cache_path <- file.path(dirname(args$out_tsv), "id_annot_cache.tsv")
}

rank_data <- hhnet_metrics

# ------------------------
# If gene_map provided, fill/override uniprot ids
# ------------------------
if (!is.null(args$gene_map) && nzchar(args$gene_map)) {
  gene_map <- read_table_auto(args$gene_map)

  if (all(c("ensembl_gene_id","uniprot_gn_id") %in% colnames(gene_map))) {
    gene_map <- normalize_id_column(gene_map, "ensembl_gene_id")
    gene_map <- normalize_id_column(gene_map, "uniprot_gn_id")
    require_non_missing(gene_map, c("ensembl_gene_id","uniprot_gn_id"), "gene_map", args$gene_map)

    rank_data <- merge(rank_data, gene_map[, c("ensembl_gene_id","uniprot_gn_id")],
                       by = "ensembl_gene_id", all.x = TRUE, sort = FALSE, suffixes = c("", ".map"))
    if ("uniprot_gn_id.map" %in% colnames(rank_data)) {
      rank_data$uniprot_gn_id <- ifelse(!is.na(rank_data$uniprot_gn_id.map) & nzchar(rank_data$uniprot_gn_id.map),
                                        rank_data$uniprot_gn_id.map, rank_data$uniprot_gn_id)
      rank_data$uniprot_gn_id.map <- NULL
    }
    rank_data$uniprot_gn_id_primary <- primary_id(rank_data$uniprot_gn_id)

  } else {
    stop("gene_map must contain columns: ensembl_gene_id and uniprot_gn_id.", call. = FALSE)
  }

} else {
  # Map ensembl_gene_id -> uniprot + external_gene_name if not provided
  rank_data <- id_annot(
    data = rank_data,
    input_type = "ensembl_gene_id",
    convert_to = c("uniprot_gn_id"),
    cache_path = id_annot_cache_path
  )
}

# Ensure expected ID columns exist after annotation/mapping
if (!"uniprot_gn_id" %in% colnames(rank_data)) {
  rank_data$uniprot_gn_id <- NA_character_
}
if (!"external_gene_name" %in% colnames(rank_data)) {
  rank_data$external_gene_name <- NA_character_
}

rank_data$uniprot_gn_id_primary <- primary_id(rank_data$uniprot_gn_id)
rank_data$external_gene_name_primary <- primary_id(rank_data$external_gene_name)

# ------------------------
# Merge external feature tables (always appended, not ranked unless user includes them)
# ------------------------

rank_data <- merge(rank_data, druggability[, c("uniprot_gn_id","highest_score")],
                   by.x = "uniprot_gn_id_primary", by.y = "uniprot_gn_id", all.x = TRUE, sort = FALSE)

rank_data <- merge(rank_data, ml_scores[, c("Protein","Prediction_Score_rf")],
                   by.x = "uniprot_gn_id_primary", by.y = "Protein", all.x = TRUE, sort = FALSE)

rank_data <- merge(rank_data, citations[, c("symbol","counts")],
                   by.x = "external_gene_name_primary", by.y = "symbol", all.x = TRUE, sort = FALSE)

# ------------------------
# Ranking (user-specified features; default excludes ML + citations)
# ------------------------

ranking_features <- trimws(strsplit(args$ranking_features, ",", fixed = TRUE)[[1]])
ranking_features <- ranking_features[nzchar(ranking_features)]
if (length(ranking_features) == 0) stop("--ranking_features resolved to an empty list.", call. = FALSE)

# Define which columns must be present for a gene to be retained.
# This is stricter than ranking-only completeness because it requires
# successful processing/merging across the appended evidence tables.
post_merge_required_cols <- unique(c(
  "ensembl_gene_id",
  "external_gene_name_primary",
  "uniprot_gn_id_primary",
  "highest_score",
  "Prediction_Score_rf",
  "counts",
  ranking_features
))

split_rank_data <- flag_incomplete_rows(rank_data, post_merge_required_cols)
rank_data_incomplete <- split_rank_data$incomplete
rank_data <- split_rank_data$complete

message(sprintf(
  "Post-merge completeness filter removed %d row(s); %d row(s) retained.",
  nrow(rank_data_incomplete), nrow(rank_data)
))

missing_feat <- setdiff(ranking_features, colnames(rank_data))
if (length(missing_feat) > 0) {
  stop(sprintf("Ranking feature(s) not found in merged table: %s", paste(missing_feat, collapse = ", ")), call. = FALSE)
}
non_numeric <- ranking_features[!vapply(rank_data[ranking_features], is.numeric, logical(1))]
if (length(non_numeric) > 0) {
  stop(sprintf("Ranking feature(s) must be numeric: %s", paste(non_numeric, collapse = ", ")), call. = FALSE)
}

rank_data <- compute_avg_rank(rank_data, ranking_features)
rank_data <- rank_data[order(rank_data$avg_rank, na.last = TRUE), ]
rownames(rank_data) <- NULL

# ------------------------
# Write outputs
# ------------------------
out_dirs <- unique(dirname(c(
  args$out_tsv,
  args$out_rds,
  args$out_incomplete_tsv,
  args$out_incomplete_rds
)))

for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write.table(rank_data, file = args$out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(rank_data, file = args$out_rds)

write.table(rank_data_incomplete, file = args$out_incomplete_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(rank_data_incomplete, file = args$out_incomplete_rds)