#!/bin/bash
# usage: build_provenance.sh <build_dir> <druggability_scores_annot.csv> <combine --min-residues> [record files...]
#
# Writes (to stdout) the provenance record of a druggability table built by update_data.nf: the table's size,
# checksum and counts, the settings, the structure downloads and checks (the record files), and the run_info.txt
# of every step, which hold the tool versions and the checksums of the scripts and inputs.
set -euo pipefail

build_dir="$1"; table="$2"; min_residues="$3"; shift 3
records=("$@")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

versions=""
for r in "${records[@]}"; do
  v="$(sed -n 's/^data_version: //p' "$r" | head -n 1)"
  if [ -z "$v" ]; then v="$(sed -n 's/^AlphaFold DB version *\(v[0-9]*\).*/AlphaFold DB \1/p' "$r" | head -n 1)"; fi
  [ -n "$v" ] && versions="${versions:+$versions, }$v"
done

stats="$(Rscript -e '
  t <- read.csv(commandArgs(TRUE)[1], stringsAsFactors = FALSE)
  named <- !is.na(t$external_gene_name) & t$external_gene_name != ""
  cat(sprintf("%s rows, %s proteins, %s gene names; %s proteins at drug_score 0; %s proteins without an Ensembl 116 gene name",
      format(nrow(t), big.mark = ","), format(length(unique(t$uniprot_id)), big.mark = ","),
      format(length(unique(t$external_gene_name[named])), big.mark = ","),
      format(length(unique(t$uniprot_id[t$drug_score == 0])), big.mark = ","),
      format(length(setdiff(unique(t$uniprot_id), unique(t$uniprot_id[named]))), big.mark = ",")))' "$table")"

cat <<EOF
druggability_scores_annot.csv: provenance record
Written by update_data.nf on $(date -u '+%Y-%m-%d %H:%M UTC')

data_version: druggability from ${versions:-unrecorded structures}, built $(date -u '+%Y-%m-%d')

TABLE
  $(wc -c < "$table" | tr -d ' ') bytes, SHA-256 $(sha256sum "$table" | cut -d' ' -f1)
  $stats
  drug_score = the higher of FP_score (fpocket druggability of the best confident pocket) and CP_score (PocketMiner
  max_hit over confident windows), each the best over the protein's AlphaFold fragments and SWISS-MODEL model.
  combine_scores.R --min-residues $min_residues --fpocket-score druggability_conf --pocketminer-score max_hit_conf
  Gene names: annotate_scores.R, Ensembl 116 BioMart (https://jun2026.archive.ensembl.org/biomart/martservice)

SCRIPTS (SHA-256)
$(cd "$script_dir" && sha256sum combine_scores.R annotate_scores.R package_table.R | sed 's/^/  /')

STRUCTURES
EOF
for r in "${records[@]}"; do
  echo "  ---- $(basename "$r")"
  sed 's/^/  /' "$r"
done
echo
echo "STEPS (run_info.txt of each step: settings, tool versions, script and input checksums)"
for source in alphafold swissmodel; do
  for step in structures_filtered fpocket_results pocketminer_results; do
    f="$build_dir/$source/$step/run_info.txt"
    [ -f "$f" ] || continue
    echo "  ---- $source/$step/run_info.txt"
    sed 's/^/  /' "$f"
  done
done
