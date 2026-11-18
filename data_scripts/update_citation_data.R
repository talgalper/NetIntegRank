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


# where the download, the checkpoint and the output go (relative to where the script is run)
data_dir <- "citation_data"
dir.create(data_dir, showWarnings = FALSE)

# get pubtator3 gene data (~760 MB)
options(timeout = 3600)
gene2pubtator3_file <- file.path(data_dir, "gene2pubtator3.gz")
if (!file.exists(gene2pubtator3_file)) {
  download.file("https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/gene2pubtator3.gz", gene2pubtator3_file, mode = "wb")
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
fwrite(citation_counts, file.path(data_dir, paste0("PubTator3_counts_noFilt_", format(Sys.Date(), "%Y.%m.%d"), ".csv")), quote = TRUE)


# quick checks
print(nrow(citation_counts))
print(summary(citation_counts$counts))
print(citation_counts[symbol %in% c("TP53", "TNF", "EGFR", "BRCA1")])
hist(log10(citation_counts$counts))
