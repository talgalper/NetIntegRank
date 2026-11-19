library(data.table)

#### STRING physical PPI network for NetIntegRank (assets/STRING_physical_ENSG.csv, assets/STRING_physical_geneSymbol.csv) ####

# usage: Rscript update_STRING_data.R [--links-version 12.5] [--info-version 12.0] [--out-dir STRING_data]
#   --links-version  STRING release of the physical links file (default 12.5)
#   --info-version   STRING release of the protein info and aliases files (default 12.0: STRING 12.5 reuses the
#                    12.0 files; for a release that ships its own, give the same version as --links-version)
#   --out-dir        where the downloads, the two output files and the provenance record go (default STRING_data,
#                    relative to where the script is run). Downloads already there are not fetched again.
args <- commandArgs(trailingOnly = TRUE)
opt <- list(links_version = "12.5", info_version = "12.0", out_dir = "STRING_data")
if (length(args) %% 2 != 0) stop("arguments come in --name value pairs; see the usage at the top of the script")
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown option ", args[i], "; see the usage at the top of the script")
  opt[[key]] <- args[i + 1]
}
for (k in c("links_version", "info_version")) {
  if (!grepl("^[0-9]+(\\.[0-9]+)*$", opt[[k]])) stop("--", gsub("_", "-", k), " must be a version such as 12.5")
}

# where the downloads and the two output files go
data_dir <- opt$out_dir
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# get STRING data: human (9606) physical links with subscores, and the protein info and aliases files
options(timeout = 3600)
string_base <- "https://stringdb-downloads.org/download"
string_urls <- c(links   = sprintf("%s/protein.physical.links.full.v%s/9606.protein.physical.links.full.v%s.txt.gz",
                                   string_base, opt$links_version, opt$links_version),
                 info    = sprintf("%s/protein.info.v%s/9606.protein.info.v%s.txt.gz",
                                   string_base, opt$info_version, opt$info_version),
                 aliases = sprintf("%s/protein.aliases.v%s/9606.protein.aliases.v%s.txt.gz",
                                   string_base, opt$info_version, opt$info_version))
string_files <- file.path(data_dir, basename(string_urls))
names(string_files) <- names(string_urls)
for (i in seq_along(string_urls)) {
  if (!file.exists(string_files[i])) {
    # download beside the final name and rename at the end, so a broken download is never taken for a finished one
    download.file(string_urls[i], paste0(string_files[i], ".part"), mode = "wb")
    file.rename(paste0(string_files[i], ".part"), string_files[i])
  }
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


#### provenance record: sources, sizes, checksums and counts, next to the outputs ####
sha256 <- function(f) {
  tool <- if (nzchar(Sys.which("sha256sum"))) c("sha256sum") else c("shasum", "-a", "256")
  sub(" .*", "", system2(tool[1], c(tool[-1], shQuote(f)), stdout = TRUE))
}
file_line <- function(f) sprintf("  %-40s %12s bytes  SHA-256 %s", basename(f), format(file.size(f), big.mark = ","), sha256(f))
out_files <- file.path(data_dir, c("STRING_physical_ENSG.csv", "STRING_physical_geneSymbol.csv"))
genes <- unique(c(STRING_edge_ensemblGene$protein1_ENSG, STRING_edge_ensemblGene$protein2_ENSG))
writeLines(c(
  "STRING physical network: provenance record",
  sprintf("Written by update_STRING_data.R on %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z", tz = "UTC")),
  "",
  sprintf("data_version: STRING v%s", opt$links_version),
  sprintf("STRING release: links v%s; protein info and aliases v%s", opt$links_version, opt$info_version),
  "Cut-off: combined_score >= 400 (stored divided by 1000); each interaction listed once",
  "",
  "Downloads (organism 9606):",
  paste0("  ", string_urls),
  sapply(string_files, file_line),
  sprintf("  downloaded: %s", paste(format(file.mtime(string_files), "%Y-%m-%d %H:%M %Z", tz = "UTC"), collapse = "; ")),
  "",
  "Outputs:",
  sapply(out_files, file_line),
  sprintf("  STRING_physical_ENSG.csv: %s interactions, %s gene IDs (NA where STRING gives no Ensembl gene)",
          format(nrow(STRING_edge_ensemblGene), big.mark = ","), format(sum(!is.na(genes)), big.mark = ",")),
  sprintf("  STRING_physical_geneSymbol.csv: %s interactions", format(nrow(STRING_edge_symbol), big.mark = ",")),
  "",
  sprintf("Software: %s; data.table %s", R.version.string, as.character(packageVersion("data.table")))
), file.path(data_dir, "STRING.provenance.txt"))
cat("Wrote", out_files, "and", file.path(data_dir, "STRING.provenance.txt"), "\n")
