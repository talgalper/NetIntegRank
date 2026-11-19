#!/bin/bash


# example commands: ./filter_structures.sh alphafold
#                   ./filter_structures.sh swissmodel -i path/to/structures_raw/ -o path/to/structures_filtered/ -t 16
#
# options: -i  raw structures dir (default: <this dir>/<source>/structures_raw)
#          -o  filtered dir       (default: <this dir>/<source>/structures_filtered)
#          -c  minimum whole-model confidence (default: MIN_STRUCT_CONF_AF or MIN_STRUCT_CONF_SM below)
#          -t  threads            (default: all CPUs)
#
# Step 1 of the druggability scoring. Scores each structure's overall confidence and writes the structures
# to score to <filtered dir>/pdb/, which both Fpocket_master.sh and pocketminer_master.sh read.
# Running it again rebuilds the filtered set from scratch. Raw structure files are never moved or changed.



### structure confidence filter ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

## whole-model confidence filter: a structure is scored only if the mean confidence of all its atoms is at
## least this. The default 0 keeps every readable structure, because the scoring steps already leave out
## pockets and windows in low-confidence regions (MIN_CONF in Fpocket_master.sh and pocketminer_master.sh).
MIN_STRUCT_CONF_AF=0     # mean pLDDT, 0-100 scale
MIN_STRUCT_CONF_SM=0     # mean QMEANDisCo local score, 0-1 scale

usage() {
  echo "Usage: $0 <alphafold|swissmodel> [-i structures_dir] [-o filtered_dir] [-c min_conf] [-t threads]"
}

## arguments: source sets the confidence scale, paths, threshold and threads are optional
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_name="$(basename "${BASH_SOURCE[0]}")"   # this script, whatever it is called
source_type="${1:-}"
case "$source_type" in
  alphafold)  min_conf="$MIN_STRUCT_CONF_AF" ;;
  swissmodel) min_conf="$MIN_STRUCT_CONF_SM" ;;
  *) usage; exit 1 ;;
esac
shift

structures_dir="$script_dir/$source_type/structures_raw"
filtered_dir="$script_dir/$source_type/structures_filtered"
threads="$(nproc)"
while getopts "i:o:c:t:" opt; do
  case "$opt" in
    i) structures_dir="$OPTARG" ;;
    o) filtered_dir="$OPTARG" ;;
    c) min_conf="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done
if ! [[ "$min_conf" =~ ^[0-9]+([.][0-9]+)?$ ]]; then echo "-c must be a number, e.g. 50 (AlphaFold) or 0.5 (SWISS-MODEL)"; exit 1; fi

## check there are structures to filter, then set up the filtered dir
if [ ! -d "$structures_dir" ]; then
  echo "$structures_dir does not exist, give the raw structures dir with -i"
  exit 1
fi
structures_dir="$(cd "$structures_dir" && pwd)"
n_structures=$(find "$structures_dir" -maxdepth 1 -name '*.pdb' | wc -l)
if [ "$n_structures" -eq 0 ]; then
  echo "$structures_dir does not contain any .pdb files, please add at least 1 to begin"
  exit 1
fi
mkdir -p "$filtered_dir/scripts_as_run"
filtered_dir="$(cd "$filtered_dir" && pwd)"

## refuse to mix sources in one filtered dir
run_info="$filtered_dir/run_info.txt"
if [ -f "$run_info" ] && ! grep -qxF "source: $source_type" "$run_info"; then
  echo "$filtered_dir holds structures from a different source, use a new filtered dir"
  exit 1
fi

## record provenance: settings, tool versions, script checksums and checksums of the raw structures
{
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "source: $source_type"
  echo "structures_dir: $structures_dir"
  echo "n_structures: $n_structures"
  echo "min_conf: $min_conf"
  echo "threads: $threads"
  echo "R: $(Rscript -e 'cat(R.version.string)')"
  echo "bio3d: $(Rscript -e 'cat(as.character(packageVersion("bio3d")))')"
  sha256sum "$script_dir/$script_name" "$script_dir/src/struct_conf.R" | sed 's/^/sha256: /'
} > "$run_info"
find "$structures_dir" -maxdepth 1 -name '*.pdb' -printf '%f\0' \
  | (cd "$structures_dir" && xargs -0 -n 500 -P "$threads" sha256sum) \
  | sort -k 2 > "$filtered_dir/input_sha256.txt"

## keep a copy of the scripts exactly as they were for this run, next to the filtered set
cp "$script_dir/$script_name" "$script_dir/src/struct_conf.R" "$filtered_dir/scripts_as_run/"

## score and filter: each structure's confidence is the mean of the
## per-atom confidence values stored in the file (pLDDT for AlphaFold, QMEAN local for SWISS-MODEL).
## Structures at or above min_conf go into pdb/: AlphaFold as links to the raw files, SWISS-MODEL
## as trimmed copies (protein atoms only, first chain, no insertion codes, no header). The scores of
## every structure go to struct_score.csv (passed) and low_conf_struct.csv (below the threshold).
## Takes about a minute for the full AlphaFold set.
echo "Scoring structure confidence for $n_structures structures ($threads threads)"
Rscript "$script_dir/src/struct_conf.R" "$source_type" "$structures_dir" "$filtered_dir" "$min_conf" "$threads"

## summary
n_pass=$(find "$filtered_dir/pdb" -maxdepth 1 -name '*.pdb' | wc -l)
if [ "$n_pass" -eq 0 ]; then
  echo "No structures passed the confidence filter (min_conf $min_conf), nothing for fpocket or PocketMiner to do"
  exit 1
fi
echo "Finished in $(( SECONDS / 60 )) min $(( SECONDS % 60 )) s. $n_pass structures ready for fpocket and PocketMiner in $filtered_dir/pdb"
