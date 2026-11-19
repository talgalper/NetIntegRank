#!/bin/bash


# example commands: ./promote_data.sh ../reference_data/2027-03-01 --datasets string
#                   ./promote_data.sh ../reference_data/2027-03-01 --datasets string,citations,druggability --yes
#
# For the maintainer, when releasing: makes datasets built by update_data.nf the new bundled files in assets/.
# Users never need this: they use a new set with main.nf --reference_dir <set folder>.
#
# options: --datasets  which datasets to promote: string, citations, druggability (comma-separated) or all
#                      (default: every dataset the set folder's MANIFEST.tsv says was built by update_data.nf)
#          --label <dataset>=<label>
#                      required for each dataset promoted: the version of the bundled file being replaced, used to
#                      name it in legacy_data/, e.g. --label string=v12.5 --label citations=2026-11-13
#                      (see "Bundled reference data" in the main README for the versions now bundled)
#          --yes       do not ask before changing assets/
#          --accept-changes
#                      promote even if a table lost more than 20% of its rows or 5 points of gene-name coverage
#                      (without it, such a table stops the promotion: it usually means a broken build)
#
# What it does, for each dataset:
#   1. Checks the new table against the bundled one: same columns, row counts, gene-name coverage. Prints a summary,
#      and stops if a table lost more than 20% of its rows or 5 points of gene-name coverage (--accept-changes).
#   2. After you confirm: moves the bundled file to assets/legacy_data/<name>_<label>.<ext> and adds a row to
#      legacy_data/README.md.
#   3. Copies the new file into assets/. The provenance records stay in the set folder, outside git: they hold
#      local paths of the machine that built the set. The README's "Bundled reference data" section is the public
#      description of the bundled files.
#   4. Prints what is left to do by hand: the README's "Bundled reference data" section, the tests and test
#      snapshot, the version history, the commit.
# It stops before changing anything if a check fails or a file it would write is already there.



### promote a reference data set into assets/ ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
assets="$(cd "$script_dir/../assets" && pwd)"
legacy="$assets/legacy_data"
die() { echo "ERROR: $*" >&2; exit 1; }

set_dir="${1:-}"
[ -n "$set_dir" ] && [ "${set_dir#--}" = "$set_dir" ] || die "usage: $0 <set folder> [--datasets ...] [--label dataset=label] [--yes]"
shift
[ -d "$set_dir" ] || die "$set_dir is not a folder"
set_dir="$(cd "$set_dir" && pwd)"
manifest="$set_dir/provenance/MANIFEST.tsv"
[ -f "$manifest" ] || die "$set_dir has no provenance/MANIFEST.tsv: is it a folder made by update_data.nf?"

datasets=""; yes=0; accept=0; declare -A label
while [ $# -gt 0 ]; do
  case "$1" in
    --datasets) datasets="$2"; shift 2 ;;
    --label) [[ "$2" == *=* ]] || die "--label takes dataset=label"; label["${2%%=*}"]="${2#*=}"; shift 2 ;;
    --yes) yes=1; shift ;;
    --accept-changes) accept=1; shift ;;
    *) die "unknown option $1" ;;
  esac
done

# the files of each dataset and the name of its provenance record
declare -A files=( [string]="STRING_physical_ENSG.csv STRING_physical_geneSymbol.csv"
                   [citations]="PubTator3_counts.csv" [druggability]="druggability_scores_annot.csv" )
declare -A record=( [string]="STRING.provenance.txt" [citations]="citations.provenance.txt"
                    [druggability]="druggability.provenance.txt" )
built_in_set() { awk -F'\t' -v f="$1" '$1 == f && $4 ~ /^built by update_data.nf/ {found = 1} END {exit !found}' "$manifest"; }
if [ -z "$datasets" ]; then
  for d in string citations druggability; do
    first="${files[$d]%% *}"
    if built_in_set "$first"; then datasets="${datasets:+$datasets,}$d"; fi
  done
  [ -n "$datasets" ] || die "$manifest lists no dataset built by update_data.nf"
fi
[ "$datasets" = all ] && datasets="string,citations,druggability"
IFS=',' read -r -a todo <<< "$datasets"


## 1. checks and summary (nothing is changed yet)
echo "Promoting from $set_dir into $assets"
plan=(); suspicious=0
for d in "${todo[@]}"; do
  [ -n "${files[$d]:-}" ] || die "unknown dataset '$d'"
  new_rec="$set_dir/provenance/${record[$d]}"
  [ -f "$new_rec" ] || die "no ${record[$d]} in $set_dir/provenance"
  old_label="${label[$d]:-}"
  [ -n "$old_label" ] || die "give the version of the bundled $d files, used to name them in legacy_data/: --label $d=<label>, e.g. --label string=v12.5"
  [[ "$old_label" =~ ^[A-Za-z0-9._-]+$ ]] || die "--label $d=$old_label: use letters, digits, '.', '_' or '-' only (it becomes part of a file name)"
  new_version="$(sed -n 's/^data_version: //p' "$new_rec" | head -n 1)"
  echo
  echo "== $d: bundled ($old_label) -> new ($new_version)"
  for f in ${files[$d]}; do
    built_in_set "$f" || die "$f in $set_dir was not built by update_data.nf (see MANIFEST.tsv)"
    [ -f "$assets/$f" ] || die "no bundled $f to replace"
    retired="$legacy/${f%.*}_${old_label}.${f##*.}"
    [ ! -e "$retired" ] || die "$retired already exists"
    [ "$(head -n 1 "$assets/$f")" = "$(head -n 1 "$set_dir/$f")" ] \
      || die "$f: the new file's columns differ from the bundled file's: $(head -n 1 "$set_dir/$f")"
    Rscript --vanilla -e '
      a <- commandArgs(TRUE); old <- read.csv(a[1], stringsAsFactors = FALSE); new <- read.csv(a[2], stringsAsFactors = FALSE)
      g <- intersect(c("external_gene_name", "symbol", "protein1_ENSG", "protein1_symbol"), names(new))[1]
      named <- function(t) if (is.na(g) || nrow(t) == 0) NA else round(100 * mean(!is.na(t[[g]]) & t[[g]] != ""), 1)
      change <- 100 * (nrow(new) / nrow(old) - 1)
      bad <- change < -20 || isTRUE(named(new) < named(old) - 5)
      cat(sprintf("   %-32s rows %s -> %s (%+.1f%%); rows with %s: %s%% -> %s%%%s\n", a[3],
          format(nrow(old), big.mark = ","), format(nrow(new), big.mark = ","), change,
          g, named(old), named(new), if (bad) "   <-- CHECK THIS" else ""))
      quit(status = if (bad) 3 else 0)' "$assets/$f" "$set_dir/$f" "$f" || suspicious=1
    plan+=("$d|$f|$retired")
  done
done

echo
if [ "$suspicious" -eq 1 ] && [ "$accept" -ne 1 ]; then
  die "a table above lost more than 20% of its rows or 5 points of gene-name coverage. Check the build (its
       provenance records and logs); if the change is real, run again with --accept-changes. Nothing changed."
fi
if [ "$yes" -ne 1 ]; then
  read -r -p "Move the bundled files to legacy_data/ and copy the new ones into assets/? [y/N] " answer
  [ "$answer" = y ] || [ "$answer" = Y ] || { echo "Nothing changed."; exit 0; }
fi

## 2-3. retire the bundled files, then copy the new ones in
for item in "${plan[@]}"; do
  IFS='|' read -r d f retired <<< "$item"
  mv "$assets/$f" "$retired"
  cp "$set_dir/$f" "$assets/$f"
  contents="$d, ${label[$d]}"
  # the new row goes after the last row of the README's table
  row="$(printf '| `%s` | `assets/%s` (until %s) | %s |' "$(basename "$retired")" "$f" "$(date -u +%Y-%m-%d)" "$contents")"
  awk -v row="$row" '{lines[NR] = $0; if ($0 ~ /^\|/) last = NR}
       END {for (i = 1; i <= NR; i++) {print lines[i]; if (i == last) print row}; if (!last) print row}' \
    "$legacy/README.md" > "$legacy/README.md.new" && mv "$legacy/README.md.new" "$legacy/README.md"
  echo "  $f: bundled file -> legacy_data/$(basename "$retired"); new file copied in"
done

## 4. what is left to do by hand
cat <<EOF

Done. Still to do by hand:
  1. README.md, "Bundled reference data": the new versions, dates, counts and sources (from the provenance
     records in $set_dir/provenance/).
  2. legacy_data/README.md: fill in the Contents of the rows added at the end of its table.
  3. Tests: druggability and citation files feed the end-to-end test, so its snapshot will differ. Compare the
     outputs, then re-record tests/main.nf.test.snap (nf-test test tests/main.nf.test --profile conda
     --update-snapshot). STRING does not affect the tests (they use tests/fixtures/PPI_test.csv).
  4. README "Version history": a line for the data update.
  5. Commit assets/ and assets/legacy_data/ together. Keep the set folder: its provenance records are the full
     record of how the new files were built.
EOF
