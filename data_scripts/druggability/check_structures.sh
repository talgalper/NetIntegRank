#!/bin/bash


# example commands: ./check_structures.sh alphafold alphafold/structures_raw
#                   ./check_structures.sh swissmodel my_models/ -t 16
#
# Checks a folder of structures before it is scored, so that a structure set from elsewhere (your own, or a new
# release) is in the form the scoring steps expect:
#   - file names: AlphaFold  AF-<UniProt ID>-F<fragment>-model_v<version>.pdb, one AlphaFold DB version only;
#                 SWISS-MODEL <UniProt ID>.pdb, one structure per protein
#   - every file is a PDB file with ATOM records
#   - per-residue confidence in the B-factor column of the first chain (the one scored), on the scale the source
#     uses: AlphaFold pLDDT 0-100, SWISS-MODEL QMEANDisCo local score 0-1 (scaled x 100 in the scoring steps)
# SWISS-MODEL models with several chains or insertion codes are fine: filter_structures.sh keeps the first chain
# and removes insertion codes. Prints a summary; exits with an error if any check fails. Writes nothing.



### check a structure set ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

usage() { echo "Usage: $0 <alphafold|swissmodel> <structures_dir> [-t threads]"; }
source_type="${1:-}"; structures_dir="${2:-}"
case "$source_type" in alphafold|swissmodel) ;; *) usage; exit 1 ;; esac
[ -n "$structures_dir" ] || { usage; exit 1; }
shift 2
threads="$(nproc)"
while getopts "t:" opt; do
  case "$opt" in t) threads="$OPTARG" ;; *) usage; exit 1 ;; esac
done
[ -d "$structures_dir" ] || { echo "$structures_dir does not exist"; exit 1; }
structures_dir="$(cd "$structures_dir" && pwd)"

failed=0
fail() { echo "FAIL: $*"; failed=1; }

## file names
if [ "$source_type" = alphafold ]; then
  name_re='^AF-[A-Z0-9]+(-[0-9]+)?-F[0-9]+-model_v[0-9]+\.pdb$'
else
  name_re='^[A-Z0-9]{6}([A-Z0-9]{4})?\.pdb$'
fi
n_pdb=$(find "$structures_dir" -maxdepth 1 -name '*.pdb' | wc -l)
n_other=$(find "$structures_dir" -maxdepth 1 -type f ! -name '*.pdb' | wc -l)
echo "structures_dir: $structures_dir"
echo "source: $source_type"
echo "pdb_files: $n_pdb (other files, ignored: $n_other)"
[ "$n_pdb" -gt 0 ] || { fail "no .pdb files"; exit 1; }
bad_names=$(find "$structures_dir" -maxdepth 1 -name '*.pdb' -printf '%f\n' | grep -E -v "$name_re" || true)
if [ -n "$bad_names" ]; then
  fail "$(echo "$bad_names" | wc -l) file names do not follow the $source_type naming, e.g. $(echo "$bad_names" | head -n 3 | tr '\n' ' ')"
fi
if [ "$source_type" = alphafold ]; then
  versions=$(find "$structures_dir" -maxdepth 1 -name 'AF-*.pdb' -printf '%f\n' | sed -nE 's/.*-model_(v[0-9]+)\.pdb$/\1/p' | sort | uniq -c)
  echo "alphafold_versions: $(echo "$versions" | awk '{printf "%s%s (%s)", (NR > 1 ? ", " : ""), $2, $1}')"
  [ "$(echo "$versions" | grep -c .)" -le 1 ] || fail "more than one AlphaFold DB version in one folder"
fi

## contents: ATOM records and the confidence scale (B-factor column, characters 61-66 of ATOM lines), in the
## first chain only, the one that is scored (SWISS-MODEL files can add chains copied from the template, such as
## bound peptides, which keep the template's crystallographic B-factors)
summary=$(find "$structures_dir" -maxdepth 1 -name '*.pdb' -print0 \
  | xargs -0 -r -n 500 -P "$threads" awk '
      # one line per file, flushed at once: parallel workers share the output, and whole short lines do not mix
      FNR == 1 {if (f != "") {print f, n, lo, hi; fflush()}; f = FILENAME; n = 0; lo = ""; hi = ""; chain = ""}
      /^ATOM  / {if (chain == "") chain = substr($0, 22, 1); if (substr($0, 22, 1) != chain) next
                 b = substr($0, 61, 6) + 0; n++; if (lo == "" || b < lo) lo = b; if (hi == "" || b > hi) hi = b}
      END {if (f != "") {print f, n, lo, hi; fflush()}}')
no_atoms=$(echo "$summary" | awk '$2 == 0' | wc -l)
[ "$no_atoms" -eq 0 ] || fail "$no_atoms files have no ATOM records, e.g. $(echo "$summary" | awk '$2 == 0 {print $1; exit}')"
if [ "$source_type" = alphafold ]; then
  out_of_scale=$(echo "$summary" | awk '$2 > 0 && ($3 < 0 || $4 > 100 || $4 <= 1)' | wc -l)
  scale="pLDDT, 0-100"
else
  out_of_scale=$(echo "$summary" | awk '$2 > 0 && ($3 < 0 || $4 > 1 || $4 == 0)' | wc -l)
  scale="QMEANDisCo local score, 0-1"
fi
echo "confidence_scale_expected: $scale"
echo "$summary" | awk '$2 > 0 {if (lo == "" || $3 < lo) lo = $3; if (hi == "" || $4 > hi) hi = $4}
                       END {printf "confidence_range_found: %s to %s\n", lo, hi}'
[ "$out_of_scale" -eq 0 ] || fail "$out_of_scale files have B-factors outside the $scale scale (or none above 0), e.g. $(echo "$summary" | awk -v s="$source_type" '$2 > 0 && ((s == "alphafold" && ($3 < 0 || $4 > 100 || $4 <= 1)) || (s == "swissmodel" && ($3 < 0 || $4 > 1 || $4 == 0))) {print $1; exit}')"

if [ "$failed" -eq 0 ]; then echo "result: all checks passed"; else echo "result: checks failed"; exit 1; fi
