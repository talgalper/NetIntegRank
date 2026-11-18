library(data.table)

#### STRING physical PPI network for NetIntegRank (assets/STRING_physical_ENSG.csv, assets/STRING_physical_geneSymbol.csv) ####

# where the downloads and the two output files go (relative to where the script is run)
data_dir <- "STRING_data"
dir.create(data_dir, showWarnings = FALSE)

# get STRING data: human (9606) physical links with subscores (v12.5), and the protein info and aliases files
# (STRING 12.5 reuses the 12.0 info and aliases files, so those are the 12.0 downloads)
options(timeout = 3600)
string_urls <- c(links   = "https://stringdb-downloads.org/download/protein.physical.links.full.v12.5/9606.protein.physical.links.full.v12.5.txt.gz",
                 info    = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
                 aliases = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz")
string_files <- file.path(data_dir, basename(string_urls))
names(string_files) <- names(string_urls)
for (i in seq_along(string_urls)) {
  if (!file.exists(string_files[i])) download.file(string_urls[i], string_files[i], mode = "wb")
}


STRING_edge <- fread(cmd = paste("gzip -dc", shQuote(string_files["links"])))
# STRING lists every interaction twice (A-B and B-A): keep one row per pair
STRING_edge <- STRING_edge[protein1 < protein2]
STRING_edge <- STRING_edge[, c(1,2,10)]
STRING_info <- fread(cmd = paste("gzip -dc", shQuote(string_files["info"])))
STRING_info <- STRING_info[, c(1,2)]
colnames(STRING_info) <- c("ENSP", "symbol")

STRING_edge_symbol <- merge.data.table(STRING_edge, STRING_info, by.x = "protein1", by.y = "ENSP", all.x = T)
colnames(STRING_edge_symbol)[4] <- "protein1_symbol"
STRING_edge_symbol <- merge.data.table(STRING_edge_symbol, STRING_info, by.x = "protein2", by.y = "ENSP", all.x = T, sort = F)
colnames(STRING_edge_symbol)[5] <- "protein2_symbol"
STRING_edge_symbol <- subset(STRING_edge_symbol, select = c("protein1_symbol", "protein2_symbol", "combined_score"))
STRING_edge_symbol <- STRING_edge_symbol[STRING_edge_symbol$combined_score >= 400, ]
STRING_edge_symbol$combined_score <- STRING_edge_symbol$combined_score / 1000

STRING_aliases <- fread(cmd = paste("gzip -dc", shQuote(string_files["aliases"])))
STRING_aliases <- STRING_aliases[STRING_aliases$source == "Ensembl_gene", ]
STRING_aliases <- STRING_aliases[, c(1,2)]

STRING_edge_ensemblGene <- merge.data.table(STRING_edge, STRING_aliases, by.x = "protein1", by.y = "#string_protein_id", all.x = T)
colnames(STRING_edge_ensemblGene)[4] <- "protein1_ENSG"
STRING_edge_ensemblGene <- merge.data.table(STRING_edge_ensemblGene, STRING_aliases, by.x = "protein2", by.y = "#string_protein_id", all.x = T, sort = F)
colnames(STRING_edge_ensemblGene)[5] <- "protein2_ENSG"
STRING_edge_ensemblGene <- subset(STRING_edge_ensemblGene, select = c("protein1_ENSG", "protein2_ENSG", "combined_score"))
STRING_edge_ensemblGene <- STRING_edge_ensemblGene[STRING_edge_ensemblGene$combined_score >= 400, ]
STRING_edge_ensemblGene$combined_score <- STRING_edge_ensemblGene$combined_score / 1000

# same layout as the files in assets/ (quoted IDs and header)
fwrite(STRING_edge_symbol, file.path(data_dir, "STRING_physical_geneSymbol.csv"), quote = TRUE)
fwrite(STRING_edge_ensemblGene, file.path(data_dir, "STRING_physical_ENSG.csv"), quote = TRUE)
