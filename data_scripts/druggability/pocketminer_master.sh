#!/bin/bash


# example commands: ./pocketminer_master.sh alphafold
#                   ./pocketminer_master.sh swissmodel -i path/to/structures/ -o path/to/results/ -t 16
#
# options: -i  structures dir (default: <this dir>/<source>/structures_filtered/pdb, made by filter_structures.sh)
#          -o  results dir    (default: <this dir>/<source>/pocketminer_results)
#          -t  threads        (default: all CPUs)
#          -m  minimum chain length for the main table   (default MIN_RESIDUES below)
#          -w  residues averaged on either side          (default WINDOW below)
#          -s  hit threshold for an averaged score       (default HIT_THRESHOLD below)
#          -c  confidence cut-off for a confident window (default MIN_CONF below)
#
# Re-runnable: structures that already have a PocketMiner prediction in the results dir are skipped,
# so a stopped run carries on where it left off. A results dir is only reused for exactly the same structure
# files (compared by checksum). Outputs packed by pack_results.sh are unpacked first. Structure files are never
# moved or changed.



### master script for PocketMiner cryptic pocket scoring ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

## scoring settings (defaults; -m, -w, -s and -c set them for one run)
MIN_RESIDUES=80      # structures shorter than this are scored but left out of the main table (small_structures.csv)
WINDOW=10            # each residue is averaged with this many residues on either side
HIT_THRESHOLD=0.7    # an averaged score at or above this counts as a cryptic pocket hit
MIN_CONF=70          # a window is confident if its mean pLDDT (AlphaFold) or QMEANDisCo x 100 (SWISS-MODEL) is at least this

## structures handed to a worker at a time (the model is loaded once per batch, not once per structure)
BATCH=100

usage() {
  echo "Usage: $0 <alphafold|swissmodel> [-i structures_dir] [-o results_dir] [-t threads] [-m min_residues]"
  echo "       [-w window] [-s hit_threshold] [-c min_conf]"
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
results_dir="$script_dir/$source_type/pocketminer_results"
threads="$(nproc)"
while getopts "i:o:t:m:w:s:c:" opt; do
  case "$opt" in
    i) structures_dir="$OPTARG" ;;
    o) results_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    m) MIN_RESIDUES="$OPTARG" ;;
    w) WINDOW="$OPTARG" ;;
    s) HIT_THRESHOLD="$OPTARG" ;;
    c) MIN_CONF="$OPTARG" ;;
    *) usage; exit 1 ;;
  esac
done
if ! [[ "$MIN_RESIDUES" =~ ^[0-9]+$ ]]; then echo "-m must be a whole number, e.g. 80"; exit 1; fi
if ! [[ "$WINDOW" =~ ^[0-9]+$ ]]; then echo "-w must be a whole number, e.g. 10"; exit 1; fi
if ! [[ "$HIT_THRESHOLD" =~ ^(0([.][0-9]+)?|1([.]0+)?)$ ]]; then echo "-s must be a number from 0 to 1, e.g. 0.7"; exit 1; fi
if ! [[ "$MIN_CONF" =~ ^[0-9]+([.][0-9]+)?$ ]]; then echo "-c must be a number on the 0-100 scale, e.g. 70"; exit 1; fi
source "$script_dir/src/packed_dirs.sh"

## check there are structures to score, then set up the results dir
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
mkdir -p "$results_dir"
results_dir="$(cd "$results_dir" && pwd)"
unpack_dir "$results_dir" txt_results      # outputs packed by pack_results.sh after an earlier run
unpack_dir "$results_dir" pocket_results
mkdir -p "$results_dir"/{txt_results,pocket_results,logs}

## PocketMiner setup: the authors' code and trained model sit in pocketminer_src/upstream/, and
## Python comes from the "pocketminer" conda env (Python 3.7, TensorFlow 2.1.0, mdtraj).
## The model is given by its prefix "pocketminer", with no .index on the end. With .index TensorFlow
## loads nothing and says nothing, so pocketminer_predict.py checks every weight after loading.
pm_dir="$script_dir/pocketminer_src"
pm_python="${POCKETMINER_PYTHON:-$HOME/miniconda3/envs/pocketminer/bin/python}"
model_prefix="$pm_dir/upstream/models/pocketminer"
if [ ! -x "$pm_python" ]; then
  echo "Python for PocketMiner not found at $pm_python"
  echo "It needs the pocketminer conda env, or set POCKETMINER_PYTHON to the right python"
  exit 1
fi
for file in "$model_prefix.index" "$model_prefix.data-00000-of-00001" "$pm_dir/upstream/src/models.py" \
            "$pm_dir/pocketminer_predict.py" "$pm_dir/pocketminer_scores.R"; do
  if [ ! -f "$file" ]; then echo "PocketMiner file missing: $file"; exit 1; fi
done
export PYTHONDONTWRITEBYTECODE=1   # stops Python leaving cache files in the code folder
if ! pm_versions="$("$pm_python" "$pm_dir/pocketminer_predict.py" versions 2> "$results_dir/logs/env_check.log")"; then
  echo "The PocketMiner Python environment did not start, see $results_dir/logs/env_check.log"
  exit 1
fi

## refuse to mix sources or structure sets in one results dir
## method: a results dir is only reused for the same source and exactly the same structure files.
## The files are compared by checksum, not by folder path, so moving or renaming the structures folder is fine.
new_sums="$results_dir/logs/input_sha256.new"
find "$structures_dir" -maxdepth 1 -name '*.pdb' -printf '%f\0' \
  | (cd "$structures_dir" && xargs -0 -n 500 -P "$threads" sha256sum) \
  | sort -k 2 > "$new_sums"
run_info="$results_dir/run_info.txt"
if [ -f "$run_info" ]; then
  if ! grep -qxF "source: $source_type" "$run_info"; then
    echo "$results_dir holds results for a different source, use a new results dir"
    rm -f "$new_sums"; exit 1
  fi
  if [ -f "$results_dir/input_sha256.txt" ] && \
     ! cmp -s <(LC_ALL=C sort "$new_sums") <(LC_ALL=C sort "$results_dir/input_sha256.txt"); then
    echo "$results_dir holds results for a different set of structures, use a new results dir"
    rm -f "$new_sums"; exit 1
  fi
fi
mv "$new_sums" "$results_dir/input_sha256.txt"

## record provenance: settings, tool versions, script and model checksums (input checksums are in input_sha256.txt)
{
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "source: $source_type"
  echo "structures_dir: $structures_dir"
  echo "n_structures: $n_structures"
  echo "threads: $threads"
  echo "min_residues: $MIN_RESIDUES"
  echo "window_each_side: $WINDOW"
  echo "hit_threshold: $HIT_THRESHOLD"
  echo "min_conf: $MIN_CONF"
  echo "model: $model_prefix"
  echo "$pm_versions"
  echo "R: $(Rscript -e 'cat(R.version.string)')"
  sha256sum "$script_dir/$script_name" "$pm_dir/pocketminer_predict.py" "$pm_dir/pocketminer_scores.R" \
    "$pm_dir"/upstream/src/*.py "$model_prefix".* | sed 's/^/sha256: /'
} > "$run_info"

## keep a copy of the three scripts exactly as they were for this run, next to the results
mkdir -p "$results_dir/scripts_as_run"
cp "$script_dir/$script_name" "$pm_dir/pocketminer_predict.py" "$pm_dir/pocketminer_scores.R" \
  "$results_dir/scripts_as_run/"

## 1. list the structures that still need a prediction (those done in an earlier run are skipped)
todo="$results_dir/logs/todo.txt"
"$pm_python" "$pm_dir/pocketminer_predict.py" todo "$source_type" "$structures_dir" "$results_dir" > "$todo"
n_todo=$(wc -l < "$todo")
n_start=$(( n_structures - n_todo ))
echo "predictions_from_earlier_runs: $n_start of $n_structures" >> "$run_info"   # reused, not recomputed
: > "$results_dir/pocketminer_failed.txt"
: > "$results_dir/logs/predict.log"

## 1. progress: structures finished (prediction files + failures) against the total, with rate and ETA
fmt_secs() {   # 40 -> "40 s", 75 -> "1 min", 4000 -> "1 h 06 min"
  if [ "$1" -ge 3600 ]; then printf '%d h %02d min' $(( $1 / 3600 )) $(( $1 % 3600 / 60 ))
  elif [ "$1" -ge 60 ]; then printf '%d min' $(( $1 / 60 ))
  else printf '%d s' "$1"; fi
}
count_done() {
  find "$results_dir/txt_results" -maxdepth 1 -name '*-predictions.txt' | wc -l
}
show_progress() {
  local n n_failed elapsed rate eta fill
  n_failed=$(wc -l < "$results_dir/pocketminer_failed.txt")
  n=$(( $(count_done) + n_failed ))
  elapsed=$(( $(date +%s) - start_time ))
  if [ "$n" -gt "$n_start" ] && [ "$elapsed" -gt 0 ]; then
    rate="$(( (n - n_start) * 60 / elapsed ))/min"
    eta="$(fmt_secs $(( (n_structures - n) * elapsed / (n - n_start) ))) left"
  else
    rate="-/min"; eta="ETA pending"
  fi
  fill=$(( 30 * n / n_structures ))
  printf '[%s%s] %d/%d (%d%%), %s, %d failed, %s' \
    "$(printf '%*s' "$fill" '' | tr ' ' '#')" "$(printf '%*s' $(( 30 - fill )) '' | tr ' ' '-')" \
    "$n" "$n_structures" $(( 100 * n / n_structures )) "$rate" "$n_failed" "$eta"
}

## 1. run PocketMiner over the structures in parallel, with a progress readout
## method: the list is cut into batches and each batch goes to a worker (one CPU each). A worker
## builds the network, loads the trained weights, checks they loaded, then gives every residue of
## each structure a score from 0 to 1. One file per structure is saved in txt_results/.
start_time=$(date +%s)
echo "Starting PocketMiner on $n_structures structures ($threads threads, $n_start already done)"
if [ -t 1 ]; then live=1; else live=0; fi   # terminal: one line redrawn every 5 s; log file: a line every 5 min
( tick=0
  while sleep 5; do   # short sleeps in both cases, so nothing is left waiting once the run ends
    tick=$(( tick + 1 ))
    if [ "$live" -eq 1 ]; then printf '\r\033[K  %s' "$(show_progress)"
    elif [ $(( tick % 60 )) -eq 0 ]; then echo "  $(date -u +%T) UTC: $(show_progress)"; fi
  done ) &
ticker=$!
trap 'kill $ticker 2> /dev/null || true' EXIT
xargs -r -a "$todo" -d '\n' -n "$BATCH" -P "$threads" \
  "$pm_python" "$pm_dir/pocketminer_predict.py" run "$source_type" "$results_dir" "$model_prefix" \
  2>> "$results_dir/logs/predict.log" || true
kill $ticker 2> /dev/null || true
wait $ticker 2> /dev/null || true
if [ "$live" -eq 1 ]; then printf '\r\033[K'; fi
echo "  PocketMiner finished in $(fmt_secs $(( $(date +%s) - start_time ))): $(show_progress | sed 's/, [^,]*$//')"

## 1. stop here if any structure was left without a prediction (for example the weights did not load,
## or a worker was killed), so the final table is never made from a partial run
n_failed=$(wc -l < "$results_dir/pocketminer_failed.txt")
n_missing=$(( n_structures - $(count_done) - n_failed ))
if [ "$n_missing" -gt 0 ]; then
  echo "PocketMiner stopped early: $n_missing structures have no prediction yet"
  grep -h "FATAL" "$results_dir/logs/predict.log" | sort -u | head -n 3 || true
  echo "See $results_dir/logs/predict.log, then run the same command again to carry on"
  exit 1
fi

## 2. format the predictions into one table per source
## method: each residue's score is averaged with its neighbours along the chain, then each structure
## gets its highest averaged score (max_hit) and its number of stretches at or above the threshold (num_hits).
## The same is done counting only confident windows (max_hit_conf, num_hits_conf), with the confidence read
## from the structures. Structures under MIN_RESIDUES residues go to small_structures.csv.
echo "Scoring PocketMiner predictions"
Rscript "$pm_dir/pocketminer_scores.R" "$results_dir" "$threads" "$MIN_RESIDUES" "$WINDOW" "$HIT_THRESHOLD" \
  "$source_type" "$structures_dir" "$MIN_CONF"

## summary
if [ "$n_failed" -gt 0 ]; then
  echo "PocketMiner failed on $n_failed structures: see pocketminer_failed.txt in $results_dir"
fi
echo "Finished in $(fmt_secs "$SECONDS"). Results saved as pocketminer_results.csv and small_structures.csv in $results_dir"
