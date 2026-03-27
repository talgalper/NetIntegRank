#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--hhnet_metrics", type = "character",
              help = "Path to HHNet-derived network metrics table (TSV/CSV)."),
  make_option("--druggability", type = "character",
              help = "Path to druggability table (TSV/CSV). Required columns: external_gene_name, drug_score."),
  make_option("--ml_scores", type = "character",
              help = "Path to ML scores table (TSV/CSV). Required columns: Protein, Prediction_Score_rf."),
  make_option("--citations", type = "character",
              help = "Path to citation counts table (TSV/CSV). Required columns: symbol, counts."),
  make_option("--gene_map", type = "character", default = NULL,
              help = "Optional mapping table (TSV/CSV). Required columns: ensembl_gene_id, uniprot_gn_id."),
  make_option("--ranking_features", type = "character",
              default = "degree,betweenness,closeness,eigen_centrality,page_rank,drug_score",
              help = "Comma-separated feature columns to use in the sensitivity ranking."),
  make_option("--negative_features", type = "character", default = "",
              help = "Optional comma-separated ranking features to subtract rather than add."),
  make_option("--step", type = "double", default = 0.1,
              help = "Weight grid step for sensitivity ranking. Default: 0.1"),
  make_option("--node_id_type", type = "character", default = "ensembl_gene_id",
              help = "Primary node identifier in hhnet_metrics. Must be one of: ensembl_gene_id, external_gene_name. Default: ensembl_gene_id"),
  make_option("--out_tsv", type = "character",
              help = "Output TSV path for final ranking."),
  make_option("--out_rds", type = "character",
              help = "Output RDS path for final ranking."),
  make_option("--id_annot_cache", type = "character", default = NULL,
              help = "Optional cache path for biomaRt ID annotations."),
  make_option("--out_incomplete_tsv", type = "character", default = NULL,
              help = "Optional output TSV path for rows removed before ranking."),
  make_option("--out_incomplete_rds", type = "character", default = NULL,
              help = "Optional output RDS path for rows removed before ranking."),
  make_option("--out_missing_gene_name_csv", type = "character", default = NULL,
              help = "Optional output CSV path for hhnet_metrics rows missing external_gene_name."),
  make_option("--out_missing_annotations_tsv", type = "character", default = NULL,
              help = "Optional output TSV path for ranked genes that could not be matched to ML scores and/or citation counts."),
  make_option("--write_rds", type = "character", default = "TRUE",
              help = "Whether to write RDS outputs for ranked and incomplete tables. TRUE or FALSE. Default: TRUE")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

required_args <- c("hhnet_metrics", "druggability", "ml_scores", "citations", "out_tsv", "out_rds")
missing_args <- required_args[vapply(required_args, function(x) is.null(args[[x]]) || !nzchar(args[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop(sprintf("Missing required arguments: %s",
               paste(sprintf("--%s", missing_args), collapse = ", ")),
       call. = FALSE)
}

parse_bool <- function(x, arg_name) {
  if (is.logical(x) && length(x) == 1) return(x)
  x <- toupper(trimws(as.character(x)))
  if (x %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (x %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop(sprintf("Invalid value for %s: %s. Use TRUE or FALSE.", arg_name, as.character(x)),
       call. = FALSE)
}

write_rds <- parse_bool(args$write_rds, "--write_rds")

valid_node_id_types <- c("ensembl_gene_id", "external_gene_name")
node_id_type <- trimws(as.character(args$node_id_type))
if (!node_id_type %in% valid_node_id_types) {
  stop(sprintf(
    "Invalid --node_id_type: %s. Valid options are: %s",
    node_id_type,
    paste(valid_node_id_types, collapse = ", ")
  ), call. = FALSE)
}
primary_id_col <- node_id_type

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE)
  }
}

require_columns <- function(df, required_cols, table_label, table_path) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Input validation failed for %s (%s). Missing required columns: %s",
                 table_label, table_path, paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }
}

require_numeric <- function(df, cols, table_label, table_path) {
  cols <- intersect(cols, colnames(df))
  bad_cols <- cols[!vapply(df[cols], is.numeric, logical(1))]
  if (length(bad_cols) > 0) {
    stop(sprintf("Input validation failed for %s (%s). Non-numeric columns where numeric expected: %s",
                 table_label, table_path, paste(bad_cols, collapse = ", ")),
         call. = FALSE)
  }
}

normalise_id_column <- function(df, colname) {
  if (!colname %in% colnames(df)) return(df)
  df[[colname]] <- trimws(as.character(df[[colname]]))
  df[[colname]][df[[colname]] == ""] <- NA_character_
  df
}

require_non_missing <- function(df, cols, table_label, table_path) {
  missing_counts <- vapply(cols, function(col) {
    sum(is.na(df[[col]]) | trimws(as.character(df[[col]])) == "")
  }, integer(1))

  bad_cols <- names(missing_counts)[missing_counts > 0]
  if (length(bad_cols) > 0) {
    bad_msg <- paste(sprintf("%s (%d missing)", bad_cols, missing_counts[bad_cols]), collapse = ", ")
    stop(sprintf("Input validation failed for %s (%s). Required columns contain missing values: %s",
                 table_label, table_path, bad_msg),
         call. = FALSE)
  }
}

collapse_unique <- function(x) {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

first_id <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_

  parts <- strsplit(replace(x, is.na(x), ""), ";", fixed = TRUE)
  vapply(parts, function(vals) {
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0) NA_character_ else vals[1]
  }, character(1))
}

left_join_first <- function(x, y, by_x, by_y = by_x, cols_y = NULL) {
  if (is.null(cols_y)) cols_y <- setdiff(colnames(y), by_y)
  idx <- match(x[[by_x]], y[[by_y]])
  for (col in cols_y) {
    x[[col]] <- y[[col]][idx]
  }
  x
}

keep_max_by_key <- function(df, key_col, value_col) {
  df <- df[!is.na(df[[key_col]]) & nzchar(df[[key_col]]), , drop = FALSE]
  df <- df[order(df[[key_col]], -df[[value_col]]), , drop = FALSE]
  df[!duplicated(df[[key_col]]), , drop = FALSE]
}

minmax01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

flag_incomplete_rows <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Cannot assess incomplete rows. Missing columns in merged table: %s",
                 paste(missing_cols, collapse = ", ")),
         call. = FALSE)
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
  missing_fields <- apply(missing_matrix, 1, function(x) paste(required_cols[which(x)], collapse = ";"))

  incomplete_rows <- df[incomplete_idx, , drop = FALSE]
  if (nrow(incomplete_rows) > 0) {
    incomplete_rows$missing_fields <- missing_fields[incomplete_idx]
  }

  list(
    complete = df[!incomplete_idx, , drop = FALSE],
    incomplete = incomplete_rows
  )
}

upsert_annotation_cache <- function(cache_path, input_type, annotation_df) {
  cache_cols <- c(input_type, "uniprot_gn_id")
  cache_df <- annotation_df[, cache_cols, drop = FALSE]

  if (file.exists(cache_path)) {
    existing <- read_table_auto(cache_path)
    if (all(cache_cols %in% colnames(existing))) {
      cache_df <- rbind(existing[, cache_cols, drop = FALSE], cache_df)
    }
  }

  cache_df <- cache_df[!duplicated(cache_df[[input_type]]), , drop = FALSE]

  outdir <- dirname(cache_path)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  write.table(cache_df, file = cache_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cache_df
}

id_annot <- function(data, input_type = "external_gene_name", cache_path = NULL) {
  if (!input_type %in% colnames(data)) {
    stop(sprintf("id_annot input column '%s' not found in data.", input_type), call. = FALSE)
  }

  keys <- unique(as.character(data[[input_type]]))
  keys <- keys[!is.na(keys) & nzchar(keys)]

  cached <- NULL
  if (!is.null(cache_path) && nzchar(cache_path) && file.exists(cache_path)) {
    existing <- read_table_auto(cache_path)
    cache_cols <- c(input_type, "uniprot_gn_id")
    if (all(cache_cols %in% colnames(existing))) {
      cached <- existing[, cache_cols, drop = FALSE]
      cached <- cached[!duplicated(cached[[input_type]]), , drop = FALSE]
    }
  }

  missing_keys <- keys
  if (!is.null(cached)) {
    missing_keys <- setdiff(keys, cached[[input_type]])
  }

  fetched <- data.frame()
  if (length(missing_keys) > 0) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("Package 'biomaRt' is required when --gene_map is not supplied and cache is incomplete.",
           call. = FALSE)
    }

    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    bm <- biomaRt::getBM(
      attributes = c(input_type, "uniprotswissprot", "uniprot_gn_id"),
      filters = input_type,
      values = missing_keys,
      mart = ensembl
    )

    fetched <- data.frame(
      setNames(list(missing_keys, rep(NA_character_, length(missing_keys))),
               c(input_type, "uniprot_gn_id")),
      stringsAsFactors = FALSE
    )

    if (nrow(bm) > 0) {
      bm[[input_type]] <- trimws(as.character(bm[[input_type]]))
      bm$uniprotswissprot <- trimws(as.character(bm$uniprotswissprot))
      bm$uniprot_gn_id <- trimws(as.character(bm$uniprot_gn_id))

      by_key <- split(bm, bm[[input_type]])
      mapped <- lapply(by_key, function(df_sub) {
        swiss <- unique(df_sub$uniprotswissprot[!is.na(df_sub$uniprotswissprot) & nzchar(df_sub$uniprotswissprot)])
        gids  <- unique(df_sub$uniprot_gn_id[!is.na(df_sub$uniprot_gn_id) & nzchar(df_sub$uniprot_gn_id)])
        picked <- if (length(swiss) > 0) swiss else gids

        data.frame(
          setNames(list(df_sub[[input_type]][1], if (length(picked) > 0) paste(picked, collapse = ";") else NA_character_),
                   c(input_type, "uniprot_gn_id")),
          stringsAsFactors = FALSE
        )
      })

      mapped <- do.call(rbind, mapped)
      idx <- match(mapped[[input_type]], fetched[[input_type]])
      fetched$uniprot_gn_id[idx] <- mapped$uniprot_gn_id
    }

    if (!is.null(cache_path) && nzchar(cache_path)) {
      cached <- upsert_annotation_cache(cache_path, input_type, fetched)
    }
  }

  annot <- if (!is.null(cached)) cached else fetched
  if (nrow(annot) == 0) {
    data$uniprot_gn_id <- NA_character_
    return(data)
  }

  left_join_first(data, annot, by_x = input_type, by_y = input_type, cols_y = "uniprot_gn_id")
}

avg_rank_sensitivity <- function(input_data, features, negative_features = character(0), step = 0.1) {
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")

  grid_list <- lapply(seq_along(feature_names), function(i) seq(0, 1, by = step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)

  if (nrow(valid_combos) == 0) {
    stop("No valid weight combinations were generated. Check --step and the number of ranking features.",
         call. = FALSE)
  }

  num_combos <- nrow(valid_combos)
  num_genes <- nrow(input_data)

  rank_matrix <- matrix(NA_real_, nrow = num_genes, ncol = num_combos)
  colnames(rank_matrix) <- paste0("Combo_", seq_len(num_combos))

  for (i in seq_len(num_combos)) {
    current_weights <- valid_combos[i, ]
    combined_score <- rep(0, num_genes)

    for (feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- if (feature %in% negative_features) -1 else 1
      combined_score <- combined_score +
        sign_factor * current_weights[[weight_col]] * input_data[[features[[feature]]]]
    }

    order_idx <- order(combined_score, decreasing = TRUE)
    ranks <- integer(num_genes)
    ranks[order_idx] <- seq_along(combined_score)
    rank_matrix[, i] <- ranks
  }

  out <- data.frame(
    ensembl_gene_id = input_data$ensembl_gene_id,
    external_gene_name = input_data$external_gene_name,
    avg_rank = rowMeans(rank_matrix),
    rank_variance = apply(rank_matrix, 1, var),
    stringsAsFactors = FALSE
  )

  out[order(out$avg_rank), , drop = FALSE]
}

ranking_features <- trimws(strsplit(args$ranking_features, ",", fixed = TRUE)[[1]])
ranking_features <- ranking_features[nzchar(ranking_features)]
if (length(ranking_features) == 0) {
  stop("--ranking_features resolved to an empty list.", call. = FALSE)
}

negative_features <- trimws(strsplit(args$negative_features, ",", fixed = TRUE)[[1]])
negative_features <- negative_features[nzchar(negative_features)]
bad_negative <- setdiff(negative_features, ranking_features)
if (length(bad_negative) > 0) {
  stop(sprintf("--negative_features contains values not present in --ranking_features: %s",
               paste(bad_negative, collapse = ", ")),
       call. = FALSE)
}

if (args$step <= 0 || args$step > 1) {
  stop("--step must be > 0 and <= 1.", call. = FALSE)
}

hhnet_metrics <- read_table_auto(args$hhnet_metrics)
druggability  <- read_table_auto(args$druggability)
ml_scores     <- read_table_auto(args$ml_scores)
citations     <- read_table_auto(args$citations)

# Keep both columns required to preserve the existing output structure and downstream behaviour.
require_columns(hhnet_metrics, c("ensembl_gene_id", "external_gene_name"), "hhnet_metrics", args$hhnet_metrics)
require_columns(druggability, c("external_gene_name", "drug_score"), "druggability", args$druggability)
require_columns(ml_scores, c("Protein", "Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_columns(citations, c("symbol", "counts"), "citations", args$citations)

required_feature_sources <- setdiff(ranking_features, c("drug_score", "counts_norm", "counts", "Prediction_Score_rf"))
missing_hh_features <- setdiff(required_feature_sources, colnames(hhnet_metrics))
if (length(missing_hh_features) > 0) {
  stop(sprintf("Ranking feature(s) not found in hhnet_metrics: %s",
               paste(missing_hh_features, collapse = ", ")),
       call. = FALSE)
}

require_numeric(hhnet_metrics, intersect(required_feature_sources, colnames(hhnet_metrics)), "hhnet_metrics", args$hhnet_metrics)
require_numeric(druggability, c("drug_score"), "druggability", args$druggability)
require_numeric(ml_scores, c("Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_numeric(citations, c("counts"), "citations", args$citations)

hhnet_metrics <- normalise_id_column(hhnet_metrics, "ensembl_gene_id")
hhnet_metrics <- normalise_id_column(hhnet_metrics, "external_gene_name")
druggability  <- normalise_id_column(druggability, "external_gene_name")
ml_scores     <- normalise_id_column(ml_scores, "Protein")
citations     <- normalise_id_column(citations, "symbol")

message(sprintf(
  "Using primary node ID type: %s",
  primary_id_col
))

# Only the selected primary node ID is hard-required to be fully populated here.
# This is the key behavioural change.
require_non_missing(hhnet_metrics, c(primary_id_col), "hhnet_metrics", args$hhnet_metrics)
require_non_missing(druggability, c("drug_score"), "druggability", args$druggability)
require_non_missing(ml_scores, c("Protein", "Prediction_Score_rf"), "ml_scores", args$ml_scores)
require_non_missing(citations, c("symbol", "counts"), "citations", args$citations)

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
if (is.null(args$out_missing_gene_name_csv) || !nzchar(args$out_missing_gene_name_csv)) {
  args$out_missing_gene_name_csv <- sub("\\.tsv$", "_missing_external_gene_name.csv", args$out_tsv, ignore.case = TRUE)
  if (identical(args$out_missing_gene_name_csv, args$out_tsv)) {
    args$out_missing_gene_name_csv <- paste0(args$out_tsv, "_missing_external_gene_name.csv")
  }
}
if (is.null(args$out_missing_annotations_tsv) || !nzchar(args$out_missing_annotations_tsv)) {
  args$out_missing_annotations_tsv <- sub("\\.tsv$", "_missing_annotations.tsv", args$out_tsv, ignore.case = TRUE)
  if (identical(args$out_missing_annotations_tsv, args$out_tsv)) {
    args$out_missing_annotations_tsv <- paste0(args$out_tsv, "_missing_annotations.tsv")
  }
}

id_annot_cache_path <- args$id_annot_cache
if (is.null(id_annot_cache_path) || !nzchar(id_annot_cache_path)) {
  id_annot_cache_path <- file.path(dirname(args$out_tsv), "id_annot_cache.tsv")
}

missing_gene_name_idx <- is.na(hhnet_metrics$external_gene_name) |
  trimws(as.character(hhnet_metrics$external_gene_name)) == ""

missing_gene_name_rows <- hhnet_metrics[missing_gene_name_idx, , drop = FALSE]
if (nrow(missing_gene_name_rows) > 0) {
  missing_gene_name_rows$removal_reason <- "missing_external_gene_name"
}

hhnet_metrics <- hhnet_metrics[!missing_gene_name_idx, , drop = FALSE]

message(sprintf(
  "Removed %d row(s) from hhnet_metrics due to missing external_gene_name; %d row(s) remain for downstream ranking.",
  nrow(missing_gene_name_rows), nrow(hhnet_metrics)
))

rank_data <- unique(hhnet_metrics)

# merge citations by external_gene_name
citations <- stats::aggregate(citations$counts,
                              by = list(citations$symbol),
                              FUN = function(x) sum(as.numeric(x), na.rm = TRUE))
colnames(citations) <- c("external_gene_name", "counts")

rank_data <- left_join_first(rank_data, citations,
                             by_x = "external_gene_name", by_y = "external_gene_name",
                             cols_y = "counts")

# merge druggability by external_gene_name
druggability <- keep_max_by_key(
  druggability[, c("external_gene_name", "drug_score"), drop = FALSE],
  key_col = "external_gene_name",
  value_col = "drug_score"
)

rank_data <- left_join_first(rank_data, druggability,
                             by_x = "external_gene_name", by_y = "external_gene_name",
                             cols_y = "drug_score")

rank_data$counts_norm <- log10(pmax(rank_data$counts, 1))

missing_rank_features <- setdiff(ranking_features, colnames(rank_data))
if (length(missing_rank_features) > 0) {
  stop(sprintf("Ranking feature(s) not found after pre-ranking merges: %s",
               paste(missing_rank_features, collapse = ", ")),
       call. = FALSE)
}

non_numeric <- ranking_features[!vapply(rank_data[ranking_features], is.numeric, logical(1))]
if (length(non_numeric) > 0) {
  stop(sprintf("Ranking feature(s) must be numeric: %s", paste(non_numeric, collapse = ", ")),
       call. = FALSE)
}

# Only enforce completeness on the selected primary ID, not the alternate ID.
required_for_ranking <- unique(c(
  primary_id_col,
  "drug_score",
  ranking_features
))

split_rank_data <- flag_incomplete_rows(rank_data, required_for_ranking)
rank_data_incomplete <- split_rank_data$incomplete
rank_data_complete <- split_rank_data$complete

message(sprintf("Pre-ranking completeness filter removed %d row(s); %d row(s) retained.",
                nrow(rank_data_incomplete), nrow(rank_data_complete)))

if (nrow(rank_data_complete) == 0) {
  stop("No rows remained after the pre-ranking completeness filter.", call. = FALSE)
}

rank_input <- rank_data_complete
normalise_cols <- setdiff(ranking_features, c("counts_norm", "Prediction_Score_rf"))
for (col in intersect(normalise_cols, colnames(rank_input))) {
  rank_input[[col]] <- minmax01(rank_input[[col]])
}

feature_map <- stats::setNames(ranking_features, ranking_features)
RS_data <- avg_rank_sensitivity(rank_input,
                                features = feature_map,
                                negative_features = negative_features,
                                step = args$step)

# Reattach counts using the chosen primary identifier.
append_cols <- rank_data_complete[, c(primary_id_col, "counts"), drop = FALSE]
append_cols <- append_cols[!duplicated(append_cols[[primary_id_col]]), , drop = FALSE]

RS_data <- left_join_first(RS_data, append_cols,
                           by_x = primary_id_col, by_y = primary_id_col,
                           cols_y = "counts")

# add uniprot only after ranking
if (!is.null(args$gene_map) && nzchar(args$gene_map)) {
  gene_map <- read_table_auto(args$gene_map)
  require_columns(gene_map, c("ensembl_gene_id", "uniprot_gn_id"), "gene_map", args$gene_map)

  gene_map <- normalise_id_column(gene_map, "ensembl_gene_id")
  gene_map <- normalise_id_column(gene_map, "uniprot_gn_id")
  require_non_missing(gene_map, c("ensembl_gene_id", "uniprot_gn_id"), "gene_map", args$gene_map)

  gene_map <- stats::aggregate(gene_map$uniprot_gn_id,
                               by = list(gene_map$ensembl_gene_id),
                               FUN = collapse_unique)
  colnames(gene_map) <- c("ensembl_gene_id", "uniprot_gn_id")

  RS_data <- left_join_first(RS_data, gene_map,
                             by_x = "ensembl_gene_id", by_y = "ensembl_gene_id",
                             cols_y = "uniprot_gn_id")
} else {
  RS_data <- id_annot(RS_data, input_type = "external_gene_name", cache_path = id_annot_cache_path)
}

RS_data$join_uniprot_gn_id <- first_id(RS_data$uniprot_gn_id)

ml_scores <- keep_max_by_key(
  ml_scores[, c("Protein", "Prediction_Score_rf"), drop = FALSE],
  key_col = "Protein",
  value_col = "Prediction_Score_rf"
)

RS_data <- left_join_first(RS_data, ml_scores,
                           by_x = "join_uniprot_gn_id", by_y = "Protein",
                           cols_y = "Prediction_Score_rf")

RS_data <- RS_data[, c("external_gene_name", "ensembl_gene_id", "uniprot_gn_id",
                       "avg_rank", "rank_variance", "counts", "Prediction_Score_rf"),
                   drop = FALSE]
RS_data <- RS_data[order(RS_data$avg_rank), , drop = FALSE]
rownames(RS_data) <- NULL

out_dirs <- unique(dirname(c(
  args$out_tsv,
  args$out_rds,
  args$out_incomplete_tsv,
  args$out_incomplete_rds,
  args$out_missing_gene_name_csv,
  args$out_missing_annotations_tsv
)))
for (d in out_dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write.table(RS_data, file = args$out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

if (write_rds) {
  saveRDS(RS_data, file = args$out_rds)
}

write.table(rank_data_incomplete, file = args$out_incomplete_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

if (write_rds) {
  saveRDS(rank_data_incomplete, file = args$out_incomplete_rds)
}

write.csv(missing_gene_name_rows,
          file = args$out_missing_gene_name_csv,
          row.names = FALSE,
          quote = TRUE)

# -------------------------
# Post-ranking annotation audit
# Genes that were ranked but could not be matched to ML scores and/or citation
# counts.
# -------------------------

annotation_check_cols <- c("counts", "Prediction_Score_rf")
present_annot_cols    <- intersect(annotation_check_cols, colnames(RS_data))

if (length(present_annot_cols) > 0) {
  missing_annot_matrix <- sapply(present_annot_cols, function(col) {
    x <- RS_data[[col]]
    if (is.character(x)) is.na(x) | trimws(x) == "" else is.na(x)
  }, simplify = "matrix")

  if (length(present_annot_cols) == 1) {
    missing_annot_matrix <- matrix(missing_annot_matrix, ncol = 1)
    colnames(missing_annot_matrix) <- present_annot_cols
  }

  missing_annot_idx    <- apply(missing_annot_matrix, 1, any)
  missing_annot_fields <- apply(missing_annot_matrix, 1, function(x)
    paste(present_annot_cols[which(x)], collapse = ";"))

  missing_annot_rows <- RS_data[missing_annot_idx, , drop = FALSE]
  if (nrow(missing_annot_rows) > 0) {
    missing_annot_rows$missing_annotation_fields <- missing_annot_fields[missing_annot_idx]
  }
} else {
  missing_annot_rows <- RS_data[integer(0), , drop = FALSE]
  missing_annot_rows$missing_annotation_fields <- character(0)
}

message(sprintf(
  "%d ranked gene(s) could not be matched to one or more annotation sources (ML scores / citations).",
  nrow(missing_annot_rows)
))

write.table(missing_annot_rows,
            file = args$out_missing_annotations_tsv,
            sep  = "\t", row.names = FALSE, quote = FALSE)