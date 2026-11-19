### Gene annotation of a combined druggability table: UniProt ID -> gene name ###
#
# usage: Rscript annotate_scores.R --in <combined.csv> --out <annotated.csv> [--url <martservice>]
#   --in     a table from combine_scores.R: first column = UniProt ID
#   --url    BioMart service; default the Ensembl 116 archive, the last Ensembl release with BioMart and the release
#            NetIntegRank names its genes with
#
# The lookup (gene name, description and biotype for each UniProt ID):
#   1. look up every protein by Swiss-Prot accession (filter uniprotswissprot)
#   2. look up every protein by any UniProt accession (filter uniprot_gn_id, includes TrEMBL)
#   3. keep all Swiss-Prot matches; use lookup 2 only for proteins with no Swiss-Prot match
#   4. remove the "[Source:...]" part of each description
#   5. merge onto the table, keeping every protein (no match: empty gene name; several genes: several rows)
# query_biomart() is the same function as in NetIntegRank's bin/run_ranking.R.


# Ensembl 116 (June 2026) is the last Ensembl release that includes BioMart, and it is only served
# from the archive site. The biomaRt package versions available through conda cannot reach it, so
# the query is sent to the BioMart web service directly. Returns one character column per
# attribute, with "" where BioMart has no value (the table biomaRt::getBM() used to return).
BIOMART_URL <- "https://jun2026.archive.ensembl.org/biomart/martservice"

query_biomart <- function(attributes, filter_name, values,
                          dataset = "hsapiens_gene_ensembl", url = BIOMART_URL, batch_size = 500) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for the BioMart lookup when the annotation cache is incomplete.",
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


#### arguments ####
args <- commandArgs(trailingOnly = TRUE)
opt <- list(url = BIOMART_URL)
for (i in seq(1, length(args), by = 2)) {
  opt[[sub("^--", "", args[i])]] <- args[i + 1]
}
if (is.null(opt[["in"]]) || is.null(opt$out)) stop("usage: see the top of annotate_scores.R (--in and --out are needed)")
input_type <- "uniprot_gn_id"
convert_to <- c("external_gene_name", "description", "gene_biotype")

# which Ensembl release the service is (for the record)
registry <- httr::content(httr::GET(paste0(opt$url, "?type=registry"), httr::timeout(60)), as = "text", encoding = "UTF-8")
release <- regmatches(registry, regexpr("Ensembl Genes [0-9]+", registry))
cat("BioMart:", opt$url, "|", if (length(release)) release else "release not found", "\n")


#### the table and the proteins to look up ####
df <- read.csv(opt[["in"]])
lookup_keys <- df[[1]]


#### lookups ####
strip_desc <- function(d) {
  d$description <- gsub("\\[.*?\\]", "", d$description)
  d
}

# 1. Swiss-Prot accessions
bm1 <- query_biomart(attributes = c("uniprotswissprot", convert_to), filter_name = "uniprotswissprot",
                     values = lookup_keys, url = opt$url)
bm1 <- bm1[bm1$uniprotswissprot != "" & !is.na(bm1$uniprotswissprot), , drop = FALSE]
names(bm1)[1] <- input_type
bm1 <- strip_desc(bm1)

# 2. any UniProt accession
bm2 <- query_biomart(attributes = c("uniprot_gn_id", convert_to), filter_name = "uniprot_gn_id",
                     values = lookup_keys, url = opt$url)
bm2 <- bm2[bm2$uniprot_gn_id != "" & !is.na(bm2$uniprot_gn_id), , drop = FALSE]
names(bm2)[1] <- input_type
bm2 <- strip_desc(bm2)

# 3. all Swiss-Prot matches, then lookup 2 for proteins without one
annot <- rbind(bm1, bm2[!bm2[[input_type]] %in% bm1[[input_type]], , drop = FALSE])

# 5. merge back onto the table, keeping every protein
annotated <- merge(annot, df, by.x = input_type, by.y = colnames(df)[1], all.y = TRUE, sort = FALSE)
write.csv(annotated, opt$out, row.names = FALSE)


#### summary ####
named <- !is.na(annotated$external_gene_name) & annotated$external_gene_name != ""
genes_per_protein <- tapply(annotated$external_gene_name[named], annotated[[input_type]][named], function(x) length(unique(x)))
proteins_per_gene <- table(annotated$external_gene_name[named])
cat(sprintf("%d rows for %d proteins: %d proteins with a gene name, %d without; %d proteins with more than one gene name; %d gene names shared by more than one protein\n",
            nrow(annotated), length(unique(annotated[[input_type]])), length(unique(annotated[[input_type]][named])),
            length(setdiff(unique(annotated[[input_type]]), annotated[[input_type]][named])),
            sum(genes_per_protein > 1), sum(proteins_per_gene > 1)))
cat("written to", opt$out, "\n")
