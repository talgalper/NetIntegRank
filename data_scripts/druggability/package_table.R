### Package an annotated druggability table for NetIntegRank ###
# usage: Rscript package_table.R --in <annotated.csv from annotate_scores.R> --out <druggability_scores_annot.csv>
#
# Gives the columns their final names and order, trims the trailing space left on descriptions by removing
# "[Source:...]", sorts the rows (highest drug_score first, then UniProt ID and gene name, so the file is the same
# whatever order BioMart answered in) and checks the two columns NetIntegRank needs.
#
# Columns:
#   uniprot_id          UniProt accession
#   external_gene_name  Ensembl gene name (NetIntegRank joins drug_score to its genes on this)
#   description         Ensembl gene description
#   gene_biotype        Ensembl gene biotype
#   FP_score            best fpocket druggability over the protein's structures, as chosen in combine_scores.R
#                       (confident pockets only with its defaults; 0 if none)
#   FP_source           AF or SM: the structure that gave FP_score
#   FP_struct_score     whole-model confidence of that structure (pLDDT, or QMEANDisCo x 100)
#   CP_score            best PocketMiner max_hit over the protein's structures, as chosen in combine_scores.R
#                       (confident windows only with its defaults; 0 if none)
#   CP_source           AF or SM: the structure that gave CP_score
#   drug_score          the higher of FP_score and CP_score
#   drug_score_tool     fpocket or PocketMiner: which gave drug_score (empty if the protein has no score)
#   drug_score_source   AF or SM: the structure that gave drug_score
# A protein with several Ensembl genes has one row per gene; a protein with no Ensembl match has an empty gene name.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
if (is.null(opt[["in"]]) || is.null(opt$out)) stop("usage: Rscript package_table.R --in <annotated.csv> --out <table.csv>")

x <- read.csv(opt[["in"]])

## new names, in their final order (old name = new name)
renames <- c(uniprot_gn_id = "uniprot_id", external_gene_name = "external_gene_name", description = "description",
             gene_biotype = "gene_biotype", druggability = "FP_score", Fpocket_method = "FP_source",
             struct_score = "FP_struct_score", CP_score = "CP_score", PocketMiner_method = "CP_source",
             highest_score = "drug_score", tool = "drug_score_tool", source = "drug_score_source")
if (!all(names(renames) %in% names(x))) stop("columns missing from ", opt[["in"]], ": ",
                                              paste(setdiff(names(renames), names(x)), collapse = ", "))
x <- x[, names(renames)]
names(x) <- unname(renames)

## trailing space left on descriptions by removing "[Source:...]"
x$description <- trimws(x$description, which = "right")

## sort: highest drug_score first, then UniProt ID, then gene name
x <- x[order(-x$drug_score, x$uniprot_id, x$external_gene_name), ]

## checks: a numeric drug_score everywhere and no fully repeated row. (A protein can map to two Ensembl genes with the
## same name, or both unnamed, that differ in biotype or description, e.g. SNHG32 as lncRNA and as protein_coding;
## those rows are kept.)
if (!is.numeric(x$drug_score) || anyNA(x$drug_score)) stop("drug_score must be numeric and never missing")
if (any(duplicated(x))) stop("fully repeated rows in ", opt[["in"]])

write.csv(x, opt$out, row.names = FALSE)
named <- !is.na(x$external_gene_name) & x$external_gene_name != ""
cat(sprintf("%d rows, %d proteins, %d gene names; drug_score 0 for %d proteins; written to %s\n",
            nrow(x), length(unique(x$uniprot_id)), length(unique(x$external_gene_name[named])),
            length(unique(x$uniprot_id[x$drug_score == 0])), opt$out))
