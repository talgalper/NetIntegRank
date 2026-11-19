# usage: Rscript swissmodel_select.R <SWISS-MODEL_Repository/INDEX> <selected_models.tsv>
#
# Lists, for each protein, the one model that swissmodel_sort.R copies, so that download_swissmodel_structures.sh
# can unpack just those models instead of the whole repository archive. The rule and the reading of INDEX are
# swissmodel_sort.R's, line for line: per UniProt accession (isoform suffix removed), the model with sequence
# identity >= 60 and the highest QMEANDisCo global score (the first one on a tie), looked for in the folder
# SWISS-MODEL_Repository/<characters 1-2>/<3-4>/<the rest of the accession>/swissmodel/ (P04637 -> P0/46/37,
# A0A075B734 -> A0/A0/75B734). swissmodel_sort.R then runs on what was
# unpacked, so it still makes the final choice; this list only decides what to unpack.
#
# Output columns: uniprot_id, dir, coordinate_id

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript swissmodel_select.R <INDEX> <selected_models.tsv>")

meta <- read.table(args[1], sep = "\t", header = T)
meta$uniprot_id <- sub("-\\d+$", "", meta$UniProtKB_ac)

# rows swissmodel_sort.R can choose from: sequence identity of at least 60 (experimental PDB entries have none)
cand <- meta[!is.na(meta$seqid) & meta$seqid >= 60, ]
rows <- split(seq_len(nrow(cand)), factor(cand$uniprot_id, levels = unique(cand$uniprot_id)))
best <- vapply(rows, function(r) {
  i <- which.max(cand$qmeandisco_global[r])
  if (length(i) == 0) stop("no QMEANDisCo score for any model of ", cand$uniprot_id[r[1]],
                           " (swissmodel_sort.R stops on this too)")
  r[i]
}, integer(1))
sel <- cand[best, ]

# folder as swissmodel_sort.R builds it: characters 1-2, 3-4 and the rest of the accession
sel$dir <- file.path("SWISS-MODEL_Repository", substr(sel$uniprot_id, 1, 2), substr(sel$uniprot_id, 3, 4),
                     substr(sel$uniprot_id, 5, nchar(sel$uniprot_id)), "swissmodel")

write.table(sel[c("uniprot_id", "dir", "coordinate_id")], args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d proteins in INDEX, %d with a model of sequence identity >= 60\n",
            length(unique(meta$uniprot_id)), nrow(sel)))
