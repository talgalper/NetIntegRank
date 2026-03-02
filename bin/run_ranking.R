### run ranking algorithm ###

# please ensure you have all the necessary packages before you start:
missing_pkgs <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) == 0) {
    message("All packages are installed.")
  } else {
    message("Missing packages: ", paste(missing, collapse = ", "))
  }
}

missing_pkgs(c("edgeR", "ggplot2", "ggrepel", "RColorBrewer", "biomaRt", "progress"))



source("R/functions.R") # load necessary functions
source("R/load_data.R") # load the dataset used in ranking (give it a min)


# add your extra data here


# run ranking function - indicate desired colnames to be included.
RS_data <- avg_rank_sensitivity(rank_data,
                                features = list(
                                  betweenness = "betweenness",
                                  centrality = "degree_centrality",
                                  druggability = "highest_score", # highest score between Fpocket and PocketMiner
                                  eigen_centrality = "eigen_centrality",
                                  closeness = "closeness",
                                  page_rank = "page_rank"), 
                                step = 0.1) # recommend keeping this default
rownames(RS_data) <- NULL # tidy rownames

# annotate ranked data with uniprot ids (requires BiomaRt package)
RS_data <- id_annot(data = RS_data, 
                    input_type = "external_gene_name", 
                    convert_to = "uniprot_gn_id")

# append ML scores and citation counts
pubtator_counts <- read.csv("data/PubTator3_counts.csv")
ML_data <- read.csv("data/ML_data.csv")
RS_data <- merge(RS_data, ML_data[, c(1,106)], by.x = "uniprot_gn_id", by.y = "Protein", all.x = T)
RS_data <- merge(RS_data, pubtator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
RS_data <- RS_data[order(RS_data$avg_rank), ]
rownames(RS_data) <- NULL

# subset low citation + high ML score genes (optional)
RS_data_subset <- RS_data[RS_data$counts < 5000, ]
RS_data_subset <- RS_data_subset[RS_data_subset$Prediction_Score_rf > 0.5, ]