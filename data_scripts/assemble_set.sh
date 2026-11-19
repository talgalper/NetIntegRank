#!/bin/bash
# usage: bash assemble_set.sh <set_dir> <reference_dir> <built_dir>
#
# Used by update_data.nf (ASSEMBLE_SET). Makes the new reference data set complete, so that it works as
# main.nf --reference_dir: for each file the pipeline reads, the copy built in this run (in <built_dir>) is put
# in the current folder for publishing into <set_dir>; a file not built in this run is kept if <set_dir> already
# holds it (from an earlier update_data.nf run into the same folder), and copied from <reference_dir> otherwise.
# Writes MANIFEST.tsv: each file's size, SHA-256 and origin.
set -euo pipefail

set_dir="$1"; ref_dir="$2"; built="$3"
files=(STRING_physical_ENSG.csv STRING_physical_geneSymbol.csv PubTator3_counts.csv druggability_scores_annot.csv ML_data.csv)
old_manifest="$set_dir/provenance/MANIFEST.tsv"
now="$(date -u '+%Y-%m-%d %H:%M UTC')"
sha() { sha256sum "$1" | cut -d' ' -f1; }

printf 'file\tbytes\tsha256\torigin\n' > MANIFEST.tsv.new
for f in "${files[@]}"; do
  if [ -f "$built/$f" ]; then
    cp "$built/$f" "$f"; src="$f"; origin="built by update_data.nf on $now"
  elif [ -f "$set_dir/$f" ]; then
    src="$set_dir/$f"; origin=""
    if [ -f "$old_manifest" ]; then
      origin="$(awk -F'\t' -v f="$f" -v s="$(sha "$src")" '$1 == f && $3 == s {print $4}' "$old_manifest")"
    fi
    [ -n "$origin" ] || origin="already in this folder (not in an earlier manifest, or changed since)"
  elif [ -f "$ref_dir/$f" ]; then
    cp "$ref_dir/$f" "$f"; src="$f"; origin="copied from $ref_dir on $now"
  else
    echo "No $f in $ref_dir to complete the set" >&2; exit 1
  fi
  printf '%s\t%s\t%s\t%s\n' "$f" "$(wc -c < "$src" | tr -d ' ')" "$(sha "$src")" "$origin" >> MANIFEST.tsv.new
done
mv MANIFEST.tsv.new MANIFEST.tsv
cat MANIFEST.tsv
