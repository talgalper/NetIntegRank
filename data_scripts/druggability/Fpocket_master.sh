#!/bin/bash


# example commands: ./Fpocket_master.sh alphafold
#                   ./Fpocket_master.sh swissmodel -i path/to/structures_filtered/pdb/ -o path/to/results/ -t 16
#
# options: -i  structures dir (default: <this dir>/<source>/structures_filtered/pdb, made by filter_structures.sh)
#          -o  results dir    (default: <this dir>/<source>/fpocket_results)
#          -t  threads        (default: all CPUs)
#          -m  minimum chain length for the main table (default MIN_RESIDUES below)
#          -c  confidence cut-off for a confident pocket (default MIN_CONF below)
#
# Step 2 of the druggability scoring (fpocket); run filter_structures.sh first. The confidence scores are
# read from struct_score.csv in the folder above the structures dir.
# fpocket: set FPOCKET_BIN to an fpocket installed by conda (it must be release FPOCKET_VERSION below);
# otherwise that release is built from source next to this script on the first run.
# Re-runnable: structures that already have fpocket output in the results dir are skipped,
# so a stopped run carries on where it left off. A results dir is only reused for exactly the same structure
# files (compared by checksum). Outputs packed by pack_results.sh are unpacked first. Structure files are never
# moved or changed.



### master script for fpocket druggability scoring ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

## fpocket release to use
FPOCKET_VERSION="4.2.3"

## table settings (defaults; -m and -c set them for one run)
MIN_RESIDUES=80      # structures shorter than this are scored but left out of the main table (small_structures.csv)
MIN_CONF=70          # a pocket is confident if its lining residues' mean pLDDT (AlphaFold) or QMEANDisCo x 100 (SWISS-MODEL) is at least this

usage() {
  echo "Usage: $0 <alphafold|swissmodel> [-i structures_dir] [-o results_dir] [-t threads] [-m min_residues] [-c min_conf]"
}

## arguments: source sets the file naming and default folders, paths and threads are optional
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_name="$(basename "${BASH_SOURCE[0]}")"   # this script, whatever it is called
source_type="${1:-}"
case "$source_type" in
  alphafold)  ;;   # files named AF-<protein>-<fragment>-model_<version>.pdb
  swissmodel) ;;   # files named <protein>.pdb
  *) usage; exit 1 ;;
esac
shift

structures_dir="$script_dir/$source_type/structures_filtered/pdb"
results_dir="$script_dir/$source_type/fpocket_results"
threads="$(nproc)"
while getopts "i:o:t:m:c:" opt; do
  case "$opt" in
    i) structures_dir="$OPTARG" ;;
    o) results_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    m) MIN_RESIDUES="$OPTARG" ;;
    c) MIN_CONF="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done
if ! [[ "$MIN_RESIDUES" =~ ^[0-9]+$ ]]; then echo "-m must be a whole number, e.g. 80"; exit 1; fi
if ! [[ "$MIN_CONF" =~ ^[0-9]+([.][0-9]+)?$ ]]; then echo "-c must be a number on the 0-100 scale, e.g. 70"; exit 1; fi
source "$script_dir/src/packed_dirs.sh"

## check there are structures to score and their confidence scores, then set up the results dir
if [ ! -d "$structures_dir" ]; then
  echo "$structures_dir does not exist, run filter_structures.sh first or give the structures dir with -i"
  exit 1
fi
structures_dir="$(cd "$structures_dir" && pwd)"
n_structures=$(find "$structures_dir" -maxdepth 1 -name '*.pdb' | wc -l)
if [ "$n_structures" -eq 0 ]; then
  echo "$structures_dir does not contain any .pdb files, please add at least 1 to begin"
  exit 1
fi
struct_score_csv="$(dirname "$structures_dir")/struct_score.csv"   # written by filter_structures.sh next to pdb/
if [ ! -f "$struct_score_csv" ]; then
  echo "No struct_score.csv next to $structures_dir, run filter_structures.sh first"
  exit 1
fi
mkdir -p "$results_dir"
results_dir="$(cd "$results_dir" && pwd)"
unpack_dir "$results_dir" structural   # outputs packed by pack_results.sh after an earlier run
mkdir -p "$results_dir"/{scores,structural,logs,tmp,scripts_as_run}

## fpocket setup: an fpocket from conda (FPOCKET_BIN), checked to be the pinned release, or else the pinned
## release built locally from source (no sudo, leaves any system fpocket alone)
fpocket_dir="$script_dir/fpocket-$FPOCKET_VERSION"
if [ -n "${FPOCKET_BIN:-}" ]; then
  fpocket_bin="$FPOCKET_BIN"
  [ -x "$fpocket_bin" ] || { echo "FPOCKET_BIN is not an executable file: $fpocket_bin"; exit 1; }
  conda_meta="$(ls "$(dirname "$(dirname "$fpocket_bin")")"/conda-meta/fpocket-*.json 2> /dev/null | head -n 1 || true)"
  fpocket_build="$(basename "${conda_meta:-unknown}" .json)"   # fpocket-<version>-<build>
  if [ "$(echo "$fpocket_build" | cut -d- -f2)" != "$FPOCKET_VERSION" ]; then
    echo "FPOCKET_BIN must be fpocket $FPOCKET_VERSION installed by conda; found '$fpocket_build' for $fpocket_bin"
    exit 1
  fi
  fpocket_build="conda $fpocket_build"
else
  fpocket_bin="$fpocket_dir/bin/fpocket"
  if [ ! -x "$fpocket_bin" ]; then
    echo "Building fpocket $FPOCKET_VERSION in $fpocket_dir"
    if ! dpkg -s libnetcdf-dev > /dev/null 2>&1; then
      echo "fpocket needs libnetcdf-dev: sudo apt-get install -y libnetcdf-dev"
      exit 1
    fi
    [ -d "$fpocket_dir" ] || git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$FPOCKET_VERSION" \
      https://github.com/Discngine/fpocket.git "$fpocket_dir"
    if ! make -C "$fpocket_dir" > "$fpocket_dir/build.log" 2>&1; then
      echo "fpocket build failed, see $fpocket_dir/build.log"
      exit 1
    fi
  fi
  fpocket_build="source, commit $(git -C "$fpocket_dir" rev-parse HEAD)"
fi

## refuse to mix sources, fpocket versions or structure sets in one results dir
## method: a results dir is only reused for the same source and exactly the same structure files, compared by
## checksum (not by file name: SWISS-MODEL files are named <protein>.pdb in every release, so a name check
## alone would keep old outputs for new models). Moving or renaming the structures folder is fine.
run_info="$results_dir/run_info.txt"
new_sums="$results_dir/logs/input_sha256.new"
find "$structures_dir" -maxdepth 1 -name '*.pdb' -printf '%f\0' \
  | (cd "$structures_dir" && xargs -0 -n 500 -P "$threads" sha256sum) \
  | sort -k 2 > "$new_sums"
if [ -f "$run_info" ]; then
  for line in "source: $source_type" "fpocket: $FPOCKET_VERSION"; do
    if ! grep -qxF "$line" "$run_info"; then
      echo "$results_dir holds results from a different run ('$line' differs), use a new results dir"
      rm -f "$new_sums"; exit 1
    fi
  done
  if [ -f "$results_dir/input_sha256.txt" ] && \
     ! cmp -s <(LC_ALL=C sort "$new_sums") <(LC_ALL=C sort "$results_dir/input_sha256.txt"); then
    echo "$results_dir holds results for a different set of structures, use a new results dir"
    rm -f "$new_sums"; exit 1
  fi
fi
mv "$new_sums" "$results_dir/input_sha256.txt"

## record provenance: settings, tool versions, script checksums (input checksums are in input_sha256.txt)
{
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "source: $source_type"
  echo "structures_dir: $structures_dir"
  echo "n_structures: $n_structures"
  echo "threads: $threads"
  echo "fpocket: $FPOCKET_VERSION"
  echo "fpocket_build: $fpocket_build"
  echo "R: $(Rscript -e 'cat(R.version.string)')"
  echo "readr: $(Rscript -e 'cat(as.character(packageVersion("readr")))')"
  echo "min_residues: $MIN_RESIDUES"
  echo "min_conf: $MIN_CONF"
  sha256sum "$script_dir/$script_name" "$script_dir/src/druggability_scores.R" "$struct_score_csv" | sed 's/^/sha256: /'
} > "$run_info"

## keep a copy of the scripts exactly as they were for this run, next to the results
cp "$script_dir/$script_name" "$script_dir/src/druggability_scores.R" "$results_dir/scripts_as_run/"

## 1. fpocket on one structure: outputs move to structural/, _info.txt to scores/
## method: fpocket finds pockets as clusters of empty spheres between the atoms and scores each one,
## including a druggability score (default settings, fpocket -f). fpocket writes its output folder next
## to the file it is given, so it is given a temporary link in tmp/ inside the results dir. That keeps
## its output out of the shared structures_filtered folder.
run_fpocket() {
  local pdb="$1" name link out
  name="$(basename "$pdb" .pdb)"
  [ -d "$results_dir/structural/${name}_out" ] && return 0   # done in an earlier run
  link="$results_dir/tmp/$name.pdb"
  out="$results_dir/tmp/${name}_out"
  rm -rf "$out"                                                # partial output from an interrupted run
  ln -sfn "$pdb" "$link"
  if "$fpocket_bin" -f "$link" > "$results_dir/logs/$name.log" 2>&1 && [ -f "$out/${name}_info.txt" ]; then
    mv "$out/${name}_info.txt" "$results_dir/scores/"         # info first, so a half-moved structure is redone
    mv "$out" "$results_dir/structural/"
    rm -f "$results_dir/logs/$name.log"
  else
    echo "$name" >> "$results_dir/fpocket_failed.txt"
  fi
  rm -f "$link"
}
export -f run_fpocket
export fpocket_bin results_dir

## 1. progress: structures finished (output folders + failures) against the total, with rate and ETA
fmt_secs() {   # 40 -> "40 s", 75 -> "1 min", 4000 -> "1 h 06 min"
  if [ "$1" -ge 3600 ]; then printf '%d h %02d min' $(( $1 / 3600 )) $(( $1 % 3600 / 60 ))
  elif [ "$1" -ge 60 ]; then printf '%d min' $(( $1 / 60 ))
  else printf '%d s' "$1"; fi
}
show_progress() {
  local n n_failed elapsed rate eta fill
  n_failed=$(wc -l < "$results_dir/fpocket_failed.txt")
  n=$(( $(find "$results_dir/structural" -mindepth 1 -maxdepth 1 -type d | wc -l) + n_failed ))
  elapsed=$(( $(date +%s) - start_time ))
  if [ "$n" -gt "$n_start" ] && [ "$elapsed" -gt 0 ]; then
    rate="$(( (n - n_start) * 60 / elapsed ))/min"
    eta="$(fmt_secs $(( (n_input - n) * elapsed / (n - n_start) ))) left"
  else
    rate="-/min"; eta="ETA pending"
  fi
  fill=$(( 30 * n / n_input ))
  printf '[%s%s] %d/%d (%d%%), %s, %d failed, %s' \
    "$(printf '%*s' "$fill" '' | tr ' ' '#')" "$(printf '%*s' $(( 30 - fill )) '' | tr ' ' '-')" \
    "$n" "$n_input" $(( 100 * n / n_input )) "$rate" "$n_failed" "$eta"
}

## 1. run fpocket over all structures in parallel, with a progress readout
n_input=$n_structures
: > "$results_dir/fpocket_failed.txt"
n_start=$(find "$results_dir/structural" -mindepth 1 -maxdepth 1 -type d | wc -l)   # done in an earlier run
echo "fpocket_outputs_from_earlier_runs: $n_start of $n_input" >> "$run_info"           # reused, not recomputed
start_time=$(date +%s)
echo "Starting fpocket $FPOCKET_VERSION on $n_input structures ($threads threads, $n_start already done)"
if [ -t 1 ]; then live=1; else live=0; fi   # terminal: one line redrawn every 5 s; log file: a line every 5 min
( tick=0
  while sleep 5; do   # short sleeps in both cases, so nothing is left waiting once the run ends
    tick=$(( tick + 1 ))
    if [ "$live" -eq 1 ]; then printf '\r\033[K  %s' "$(show_progress)"
    elif [ $(( tick % 60 )) -eq 0 ]; then echo "  $(date -u +%T) UTC: $(show_progress)"; fi
  done ) &
ticker=$!
trap 'kill $ticker 2> /dev/null || true' EXIT
find "$structures_dir" -maxdepth 1 -name '*.pdb' -print0 \
  | xargs -0 -r -n 1 -P "$threads" bash -c 'run_fpocket "$1"' _
kill $ticker 2> /dev/null || true
wait $ticker 2> /dev/null || true
if [ "$live" -eq 1 ]; then printf '\r\033[K'; fi
rmdir "$results_dir/tmp" 2> /dev/null || true   # left in place only if a failed structure left output there
echo "  fpocket finished in $(fmt_secs $(( $(date +%s) - start_time ))): $(show_progress | sed 's/, [^,]*$//')"

## 2. format scores into one table per source, adding each structure's confidence score
echo "Summarising fpocket outputs"
## method: the best pocket's scores and the number of druggable pockets per structure, then the same
## counting only confident pockets (pocket confidence = mean confidence of the residues lining it). Structures under
## MIN_RESIDUES residues go to small_structures.csv; every pocket is listed in fpocket_pockets.csv.
Rscript "$script_dir/src/druggability_scores.R" "$source_type" "$results_dir" "$threads" "$struct_score_csv" \
  "$structures_dir" "$MIN_RESIDUES" "$MIN_CONF"

## summary
n_failed=$(wc -l < "$results_dir/fpocket_failed.txt")
if [ "$n_failed" -gt 0 ]; then
  echo "fpocket failed on $n_failed structures: see fpocket_failed.txt and logs/ in $results_dir"
fi
echo "Finished in $(fmt_secs "$SECONDS"). Results saved as fpocket_druggability_full.csv, small_structures.csv and fpocket_pockets.csv in $results_dir"
