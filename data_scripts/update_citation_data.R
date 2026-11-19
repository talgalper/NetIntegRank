library(data.table)
library(httr)

#### raw (unfiltered) PubTator3 citation counts for NetIntegRank (assets/PubTator3_counts.csv) ####

# counts = the number of distinct PubMed articles in which PubTator3 found a mention of the gene. Every article in
# gene2pubtator3 counts, with no topic filter, and mentions linked to more than one gene (IDs with ";") are dropped.
# Gene names come from Ensembl 116, the release NetIntegRank names its own genes with (the ranking joins these
# counts on external_gene_name). Ensembl's human database only names human genes, so only human genes are kept.
#
# data
#   gene2pubtator3: https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/gene2pubtator3.gz (NCBI replaces this file with
#                   each PubTator3 release; the bundled counts use the file downloaded on 2026-11-13)
#   Ensembl 116 BioMart: https://jun2026.archive.ensembl.org/biomart/martservice (needs internet)
#
# usage: Rscript update_citation_data.R [--gene2pubtator3 <file>] [--out-dir citation_data] [--out <csv>]
#   --gene2pubtator3  use this copy of gene2pubtator3.gz instead of downloading NCBI's current file (for example a
#                     kept copy, to rebuild an earlier table exactly)
#   --out-dir         where the download, the checkpoint, the output and the provenance record go (default
#                     citation_data, relative to where the script is run). A download already there is used again.
#   --out             output file (default <out-dir>/PubTator3_counts_noFilt_<date>.csv)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(gene2pubtator3 = NULL, out_dir = "citation_data", out = NULL)
if (length(args) %% 2 != 0) stop("arguments come in --name value pairs; see the usage at the top of the script")
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown option ", args[i], "; see the usage at the top of the script")
  opt[[key]] <- args[i + 1]
}

# where the download, the checkpoint and the output go
data_dir <- opt$out_dir
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- if (is.null(opt$out)) file.path(data_dir, paste0("PubTator3_counts_noFilt_", format(Sys.Date(), "%Y.%m.%d"), ".csv")) else opt$out

# get pubtator3 gene data (~760 MB)
options(timeout = 3600)
gene2pubtator3_url <- "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/gene2pubtator3.gz"
server_date <- NA
if (!is.null(opt$gene2pubtator3)) {
  gene2pubtator3_file <- opt$gene2pubtator3
  if (!file.exists(gene2pubtator3_file)) stop("--gene2pubtator3 file not found: ", gene2pubtator3_file)
} else {
  gene2pubtator3_file <- file.path(data_dir, "gene2pubtator3.gz")
  server_date_file <- paste0(gene2pubtator3_file, ".server_date")   # kept beside the download for later runs
  if (!file.exists(gene2pubtator3_file)) {
    head_resp <- tryCatch(HEAD(gene2pubtator3_url, timeout(60)), error = function(e) NULL)
    if (!is.null(head_resp) && !is.null(headers(head_resp)[["last-modified"]])) {
      writeLines(headers(head_resp)[["last-modified"]], server_date_file)
    }
    # download beside the final name and rename at the end, so a broken download is never taken for a finished one
    download.file(gene2pubtator3_url, paste0(gene2pubtator3_file, ".part"), mode = "wb")
    file.rename(paste0(gene2pubtator3_file, ".part"), gene2pubtator3_file)
  }
  if (file.exists(server_date_file)) server_date <- readLines(server_date_file, n = 1)
}

# read in data (only PMID and Concept_ID, the NCBI gene ID, are needed: keeps memory down)
# (the file is plain tab-separated with no quoting; quote = "" stops fread reading quote marks inside mentions)
gene2pubtator3 <- fread(cmd = paste("gzip -dc", shQuote(gene2pubtator3_file)), sep = "\t", header = FALSE, quote = "",
                        select = c(1, 3), col.names = c("PMID", "Concept_ID"))
gene2pubtator3[, Concept_ID := as.character(Concept_ID)]

# Count number of unique PMIDs for each entrez_id (all species at this point)
citation_counts <- gene2pubtator3[, .(counts = uniqueN(PMID)), by = .(entrezgene_id = Concept_ID)]
rm(gene2pubtator3)
gc()

save(citation_counts, file = file.path(data_dir, "citation_counts_pubtator_noFilt.RData"))
load(file.path(data_dir, "citation_counts_pubtator_noFilt.RData"))


# separate multi mention-entries
multi_gene <- citation_counts[grepl(";", entrezgene_id)]
citation_counts <- citation_counts[!entrezgene_id %in% multi_gene$entrezgene_id]


# Ensembl 116 lookup, copied unchanged from NetIntegRank 2.0 (bin/run_ranking.R)
BIOMART_URL <- "https://jun2026.archive.ensembl.org/biomart/martservice"

query_biomart <- function(attributes, filter_name, values,
                          dataset = "hsapiens_gene_ensembl", url = BIOMART_URL, batch_size = 500) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for the BioMart lookup when --gene_map is not supplied and cache is incomplete.",
         call. = FALSE)
  }

  xml_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub("\"", "&quot;", x, fixed = TRUE)
  }

  empty <- as.data.frame(matrix(character(0), nrow = 0, ncol = length(attributes)),
                         stringsAsFactors = FALSE)
  colnames(empty) <- attributes

  values <- unique(as.character(values))
  if (length(values) == 0) return(empty)
  batches <- split(values, ceiling(seq_along(values) / batch_size))

  results <- lapply(batches, function(batch) {
    query <- paste0(
      "<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query>",
      "<Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"0\" uniqueRows=\"1\" count=\"\" ",
      "datasetConfigVersion=\"0.6\" completionStamp=\"1\">",
      "<Dataset name=\"", dataset, "\" interface=\"default\">",
      "<Filter name=\"", filter_name, "\" value=\"", paste(xml_escape(batch), collapse = ","), "\"/>",
      paste0("<Attribute name=\"", attributes, "\"/>", collapse = ""),
      "</Dataset></Query>"
    )

    resp <- tryCatch(
      httr::RETRY("POST", url, body = list(query = query), encode = "form",
                  httr::timeout(300), times = 3, pause_base = 5, quiet = TRUE),
      error = function(e) {
        stop(sprintf("BioMart request to %s failed: %s", url, conditionMessage(e)), call. = FALSE)
      }
    )
    body <- httr::content(resp, as = "text", encoding = "UTF-8")
    lines <- strsplit(body, "\r?\n")[[1]]
    lines <- lines[nzchar(lines)]

    # completionStamp="1" makes BioMart end a complete reply with "[success]". Anything else is an
    # HTTP error, a BioMart "Query ERROR" message or a reply that was cut off.
    complete <- length(lines) > 0 && identical(lines[length(lines)], "[success]")
    if (httr::http_error(resp) || !complete) {
      stop(sprintf("BioMart query to %s failed (HTTP %s). Reply began: %s",
                   url, httr::status_code(resp), substr(body, 1, 300)), call. = FALSE)
    }

    lines <- lines[-length(lines)]
    if (length(lines) == 0) return(empty)
    read.table(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, quote = "",
               comment.char = "", col.names = attributes, colClasses = "character",
               na.strings = character(0), check.names = FALSE, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

# every chromosome/scaffold name in Ensembl 116 (main chromosomes plus patches and haplotypes), so the lookup
# below returns every human gene with an NCBI gene ID in one request, instead of one request per 500 of the
# ~7 million IDs counted above
resp <- RETRY("GET", paste0(BIOMART_URL, "?type=filters&dataset=hsapiens_gene_ensembl"), timeout(120), times = 3)
stop_for_status(resp)
filter_list <- strsplit(content(resp, as = "text", encoding = "UTF-8"), "\r?\n")[[1]]
chromosome_line <- filter_list[startsWith(filter_list, "chromosome_name\t")]
chromosomes <- strsplit(gsub("^\\[|\\]$", "", strsplit(chromosome_line, "\t")[[1]][3]), ",")[[1]]
if (length(chromosomes) < 25) stop("Could not read the chromosome list from Ensembl 116")

entrezgene_ids <- query_biomart(attributes = c("entrezgene_id", "entrezgene_accession", "external_gene_name"),
                                filter_name = "chromosome_name",
                                values = chromosomes)
entrezgene_ids <- as.data.table(entrezgene_ids)
entrezgene_ids <- unique(entrezgene_ids[entrezgene_id != "" & external_gene_name != ""])

# where Ensembl links one NCBI gene to several genes (mostly near-identical paralogs, e.g. NCBI's DDTL to both DDT
# and DDTL), keep the gene whose name is NCBI's own symbol; if none matches, keep them all
entrezgene_ids[, ncbi_match := any(external_gene_name == entrezgene_accession), by = entrezgene_id]
entrezgene_ids <- entrezgene_ids[!ncbi_match | external_gene_name == entrezgene_accession]
entrezgene_ids <- unique(entrezgene_ids[, .(entrezgene_id, external_gene_name)])


# merge with count data (IDs Ensembl does not know, including other species, drop out here)
citation_counts <- merge.data.table(entrezgene_ids, citation_counts, by = "entrezgene_id", allow.cartesian = TRUE)

# sum counts per gene name
citation_counts <- citation_counts[, .(counts = sum(counts)), by = .(symbol = external_gene_name)]
citation_counts <- citation_counts[order(-counts, symbol)]

# same layout as assets/PubTator3_counts.csv ("symbol","counts"); copy it there under that name
fwrite(citation_counts, out_file, quote = TRUE)


# quick checks
print(nrow(citation_counts))
print(summary(citation_counts$counts))
print(citation_counts[symbol %in% c("TP53", "TNF", "EGFR", "BRCA1")])
if (interactive()) hist(log10(citation_counts$counts))


#### provenance record: source, size, checksum and counts, next to the output ####
sha256 <- function(f) {
  tool <- if (nzchar(Sys.which("sha256sum"))) c("sha256sum") else c("shasum", "-a", "256")
  sub(" .*", "", system2(tool[1], c(tool[-1], shQuote(f)), stdout = TRUE))
}
file_line <- function(f) sprintf("  %-40s %14s bytes  SHA-256 %s", basename(f), format(file.size(f), big.mark = ","), sha256(f))
download_date <- format(file.mtime(gene2pubtator3_file), "%Y-%m-%d", tz = "UTC")
writeLines(c(
  "PubTator3 citation counts: provenance record",
  sprintf("Written by update_citation_data.R on %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z", tz = "UTC")),
  "",
  sprintf("data_version: PubTator3 gene2pubtator3 of %s", download_date),
  sprintf("Source: %s", if (is.null(opt$gene2pubtator3)) gene2pubtator3_url else
            paste0("a local copy given with --gene2pubtator3 (", normalizePath(gene2pubtator3_file), ")")),
  "  NCBI replaces this file with each PubTator3 release; keep the download to rebuild this table exactly.",
  if (!is.na(server_date)) sprintf("  File date on server: %s", server_date),
  sprintf("  File date here (download date): %s", format(file.mtime(gene2pubtator3_file), "%Y-%m-%d %H:%M %Z", tz = "UTC")),
  file_line(gene2pubtator3_file),
  sprintf("Gene names: Ensembl 116 BioMart, %s", BIOMART_URL),
  "Method: distinct PubMed articles per NCBI gene ID, every article, mentions of more than one gene dropped,",
  "  human genes only (named by Ensembl 116), summed per gene name",
  "",
  "Output:",
  file_line(out_file),
  sprintf("  %s gene names; counts median %s, max %s", format(nrow(citation_counts), big.mark = ","),
          format(median(citation_counts$counts)), format(max(citation_counts$counts), big.mark = ",")),
  "",
  sprintf("Software: %s; data.table %s; httr %s", R.version.string, as.character(packageVersion("data.table")),
          as.character(packageVersion("httr")))
), file.path(data_dir, "citations.provenance.txt"))
cat("Wrote", out_file, "and", file.path(data_dir, "citations.provenance.txt"), "\n")
