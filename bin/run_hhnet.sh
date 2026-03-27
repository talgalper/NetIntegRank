#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_hhnet.sh --scores <scores.tsv> --ppi <ppi_edge_list.tsv> [options]

Required:
  --scores            Two columns: gene_id, score.
  --ppi               Edge list with gene names in columns 1-2.
                      Optional 3rd+ columns allowed (e.g. interaction score) and will be removed.

Optional:
  --run_id            If provided, write outputs under <hhnet_dir>/runs/<run_id>/{data,intermediate,results}.
                      If omitted, outputs go directly under <hhnet_dir>/{data,intermediate,results} and may overwrite prior runs.
  --outdir            Directory to stage clusters output for Nextflow (default: current directory).
  --network_name      Network label (default: STRING)
  --score_name        Score label (default: DE)
  --num_permutations  Number of score permutations (default: 100)
  --num_cores         Number of cores for GNU parallel and HHNet post-processing (default: 1)
  --compile_fortran        Fortran compilation mode: auto|always|never (default: auto)
  --hhnet_dir              Path to hierarchical_hotnet directory (default: <script_dir>/hierarchical_hotnet
                           fallback: <script_dir>/hierarchical-hotnet)
  --intermediate_cache_dir Directory to cache intermediate HHNet files keyed by (PPI, score) hash.
                           If omitted, no intermediate caching is performed.
  -h, --help               Show this help message

Outputs:
  Canonical HHNet outputs are written under either:
    (a) <hhnet_dir>/runs/<run_id>/{data,intermediate,results}     (if --run_id provided)
    (b) <hhnet_dir>/{data,intermediate,results}                  (if --run_id omitted)

  Downstream file:
    .../results/clusters_<NETWORK>_<SCORE>.tsv

  Staged copy for Nextflow:
    <outdir>/clusters_<NETWORK>_<SCORE>.tsv

Caching:
  Similarity matrix:
    Cached globally under <hhnet_dir>/cache/similarity_matrices/ using a hash of the
    indexed PPI edge list. Reusing the same PPI network reuses cached similarity_matrix.h5
    and beta.txt in later runs.

  Intermediate files (score_bins, permuted scores, hierarchies):
    If --intermediate_cache_dir is provided, these files are cached under:
      <intermediate_cache_dir>/<ppi_hash>_<score_hash>/
    The cache key is a SHA-256 hash of (indexed_edge_list, scores_0.tsv).
    A cache hit skips find_permutation_bins, permute_scores, and construct_hierarchy entirely,
    and links the cached files into the current intermediate directory before process_hierarchies.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SCORES=""
PPI=""
OUTDIR="."
RUN_ID=""

NETWORK="STRING"
SCORE="DE"
NUM_PERMUTATIONS=100
NUM_CORES=1
COMPILE_FORTRAN="auto"
INTERMEDIATE_CACHE_DIR=""

# Default HHNET_DIR resolution (prefer underscore, fallback to hyphen)
HHNET_DIR_DEFAULT_UNDERSCORE="${SCRIPT_DIR}/hierarchical_hotnet"
HHNET_DIR_DEFAULT_HYPHEN="${SCRIPT_DIR}/hierarchical-hotnet"
HHNET_DIR="$HHNET_DIR_DEFAULT_UNDERSCORE"

while [[ $# -gt 0 ]]; do
  case "$1" in
      --scores) SCORES="$2"; shift 2 ;;
    --ppi) PPI="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --run_id) RUN_ID="$2"; shift 2 ;;
    --network_name) NETWORK="$2"; shift 2 ;;
    --score_name) SCORE="$2"; shift 2 ;;
    --num_permutations) NUM_PERMUTATIONS="$2"; shift 2 ;;
    --num_cores) NUM_CORES="$2"; shift 2 ;;
    --compile_fortran) COMPILE_FORTRAN="$2"; shift 2 ;;
    --hhnet_dir) HHNET_DIR="$2"; shift 2 ;;
    --intermediate_cache_dir) INTERMEDIATE_CACHE_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$SCORES" ]] || { echo "ERROR: --scores is required" >&2; exit 2; }
[[ -n "$PPI" ]] || { echo "ERROR: --ppi is required" >&2; exit 2; }
[[ "$COMPILE_FORTRAN" =~ ^(auto|always|never)$ ]] || { echo "ERROR: --compile_fortran must be one of auto|always|never" >&2; exit 2; }

command -v python >/dev/null 2>&1 || { echo "ERROR: python not found in PATH" >&2; exit 2; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel is required but not found" >&2; exit 2; }

# Resolve HHNET_DIR robustly + support underscore/hyphen defaults.
if [[ "$HHNET_DIR" != /* ]]; then
  if [[ -d "${HHNET_DIR}/src" ]]; then
    HHNET_DIR="$(cd "$HHNET_DIR" && pwd)"
  elif [[ -d "${SCRIPT_DIR}/${HHNET_DIR}/src" ]]; then
    HHNET_DIR="$(cd "${SCRIPT_DIR}/${HHNET_DIR}" && pwd)"
  else
    echo "ERROR: --hhnet_dir does not contain src/: ${HHNET_DIR}" >&2
    exit 2
  fi
else
  if [[ ! -d "${HHNET_DIR}/src" ]]; then
    if [[ "$HHNET_DIR" == "$HHNET_DIR_DEFAULT_UNDERSCORE" ]] && [[ -d "${HHNET_DIR_DEFAULT_HYPHEN}/src" ]]; then
      HHNET_DIR="$HHNET_DIR_DEFAULT_HYPHEN"
    else
      echo "ERROR: HHNet src not found at: ${HHNET_DIR}/src" >&2
      exit 2
    fi
  fi
fi

# Sanitise RUN_ID (avoid traversal / weird chars)
SAFE_RUN_ID=""
if [[ -n "$RUN_ID" ]]; then
  SAFE_RUN_ID="${RUN_ID//[^A-Za-z0-9._-]/_}"
  if [[ "$SAFE_RUN_ID" != "$RUN_ID" ]]; then
    echo "NOTE: --run_id contained unsupported characters; using sanitised run_id: '$SAFE_RUN_ID' (original: '$RUN_ID')"
  fi
  if [[ -z "$SAFE_RUN_ID" ]]; then
    echo "ERROR: --run_id sanitised to empty; choose a different run_id" >&2
    exit 2
  fi
fi

# Decide output root
OUTPUT_ROOT="$HHNET_DIR"
if [[ -n "$SAFE_RUN_ID" ]]; then
  OUTPUT_ROOT="${HHNET_DIR}/runs/${SAFE_RUN_ID}"
  if [[ -d "$OUTPUT_ROOT" ]]; then
    echo "WARNING: run_id directory already exists:"
    echo "  $OUTPUT_ROOT"
    echo "Files in this directory may be overwritten. If you did not intend this, press Ctrl+C now to cancel."
    sleep 5
  fi
else
  echo "WARNING: No --run_id provided."
  echo "Outputs will be written into:"
  echo "  ${HHNET_DIR}/data"
  echo "  ${HHNET_DIR}/intermediate"
  echo "  ${HHNET_DIR}/results"
  echo "These may overwrite outputs from a previous run. If you did not intend this, press Ctrl+C now to cancel."
  sleep 5
fi

DATA_DIR="${OUTPUT_ROOT}/data"
INTERMEDIATE_DIR="${OUTPUT_ROOT}/intermediate"
RESULTS_DIR="${OUTPUT_ROOT}/results"

mkdir -p "$DATA_DIR" "$INTERMEDIATE_DIR" "$RESULTS_DIR"
mkdir -p "${INTERMEDIATE_DIR}/${NETWORK}" "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}"
mkdir -p "$OUTDIR"

EDGE_LIST_FILE="${DATA_DIR}/edge_list.tsv"      # indexed
INDEX_GENE_FILE="${DATA_DIR}/index_gene.tsv"    # index -> gene_id  (IMPORTANT)
GENES_FILE="${DATA_DIR}/genes_in_ppi.txt"       # one gene_id per line
SCORES0_FILE="${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv"  # gene_id, score
SIMILARITY_MATRIX_FILE="${INTERMEDIATE_DIR}/${NETWORK}/similarity_matrix.h5"
BETA_FILE="${INTERMEDIATE_DIR}/${NETWORK}/beta.txt"
SIMILARITY_CACHE_ROOT="${HHNET_DIR}/cache/similarity_matrices"

prepare_ppi_indexed() {
  local in_path="$1"
  local edge_out="$2"
  local index_out="$3"
  local genes_out="$4"

  python - "$in_path" "$edge_out" "$index_out" "$genes_out" <<'PY'
import csv
import pathlib
import sys

in_path = pathlib.Path(sys.argv[1])
edge_out = pathlib.Path(sys.argv[2])
index_out = pathlib.Path(sys.argv[3])
genes_out = pathlib.Path(sys.argv[4])

if not in_path.exists():
    print("ERROR: PPI file not found: {}".format(in_path), file=sys.stderr)
    sys.exit(2)

text = in_path.read_text(errors="replace")
if not text.strip():
    print("ERROR: PPI file is empty: {}".format(in_path), file=sys.stderr)
    sys.exit(2)

sample = text[:4096]
try:
    dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
except csv.Error:
    dialect = csv.excel_tab

rows = csv.reader(text.splitlines(), dialect)

def looks_like_header(row):
    if not row or len(row) < 2:
        return True
    a = row[0].strip().lower()
    b = row[1].strip().lower()
    header_tokens = {
        "node1","node2","node","gene","gene1","gene2",
        "protein1","protein2","source","target","from","to","a","b"
    }
    return (a in header_tokens) or (b in header_tokens)

directed_seen = set()
undirected_seen = set()
unique_edges = []

total_rows = 0
kept_rows = 0
extra_col_rows = 0
self_loop_rows = 0
directed_dup_rows = 0
reverse_dup_rows = 0
undirected_dup_rows = 0
header_skipped = False

for row in rows:
    if not row:
        continue
    if row[0].lstrip().startswith("#"):
        continue

    if total_rows == 0 and looks_like_header(row):
        header_skipped = True
        total_rows += 1
        continue

    total_rows += 1

    if len(row) >= 3 and any(x.strip() for x in row[2:]):
        extra_col_rows += 1

    a = row[0].strip()
    b = row[1].strip() if len(row) > 1 else ""
    if not a or not b:
        continue

    if a == b:
        self_loop_rows += 1
        continue

    if (a, b) in directed_seen:
        directed_dup_rows += 1
    if (b, a) in directed_seen:
        reverse_dup_rows += 1

    directed_seen.add((a, b))

    x, y = (a, b) if a <= b else (b, a)
    key = (x, y)
    if key in undirected_seen:
        undirected_dup_rows += 1
        continue

    undirected_seen.add(key)
    unique_edges.append(key)
    kept_rows += 1

if kept_rows == 0:
    print("ERROR: No valid edges found after processing PPI.", file=sys.stderr)
    sys.exit(2)

genes = sorted({g for e in unique_edges for g in e})

edge_out.parent.mkdir(parents=True, exist_ok=True)
index_out.parent.mkdir(parents=True, exist_ok=True)
genes_out.parent.mkdir(parents=True, exist_ok=True)

gene_to_idx = {g: i+1 for i, g in enumerate(genes)}

with index_out.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    for g in genes:
        w.writerow((gene_to_idx[g], g))

with genes_out.open("w") as f:
    for g in genes:
        f.write(g + "\n")

with edge_out.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    for a, b in unique_edges:
        ia = gene_to_idx[a]
        ib = gene_to_idx[b]
        if ia <= ib:
            w.writerow((ia, ib))
        else:
            w.writerow((ib, ia))

if header_skipped:
    print("PPI: detected a header row; it was skipped.")
if extra_col_rows > 0:
    print("PPI: detected {} row(s) with extra columns (e.g. interaction scores). Extra columns will be removed (keeping columns 1-2).".format(extra_col_rows))
if reverse_dup_rows > 0:
    print("PPI: detected {} reverse-duplicate edge row(s) (A B and B A). Reverse duplicates will be removed (keeping one per undirected pair).".format(reverse_dup_rows))
if directed_dup_rows > 0:
    print("PPI: detected {} exact duplicate directed edge row(s); duplicates will be removed.".format(directed_dup_rows))
if undirected_dup_rows > 0:
    print("PPI: removed {} duplicate undirected edge row(s) in total (includes reverse and exact duplicates).".format(undirected_dup_rows))
if self_loop_rows > 0:
    print("PPI: detected {} self-loop row(s) (A A); these were removed.".format(self_loop_rows))

print("PPI: kept {} unique undirected edges across {} genes.".format(kept_rows, len(genes)))
PY
}

normalize_scores_filter_to_ppi() {
  local scores_path="$1"
  local genes_file="$2"
  local out_path="$3"

  python - "$scores_path" "$genes_file" "$out_path" <<'PY'
import csv
import pathlib
import sys

scores_path = pathlib.Path(sys.argv[1])
genes_file = pathlib.Path(sys.argv[2])
out_path = pathlib.Path(sys.argv[3])

if not scores_path.exists():
    print("ERROR: score file not found: {}".format(scores_path), file=sys.stderr)
    sys.exit(2)
if not genes_file.exists():
    print("ERROR: genes_in_ppi file not found: {}".format(genes_file), file=sys.stderr)
    sys.exit(2)

genes_in_ppi = set(g.strip() for g in genes_file.read_text().splitlines() if g.strip())
if not genes_in_ppi:
    print("ERROR: genes_in_ppi is empty; cannot filter score file.", file=sys.stderr)
    sys.exit(2)

text = scores_path.read_text(errors="replace")
if not text.strip():
    print("ERROR: score file is empty: {}".format(scores_path), file=sys.stderr)
    sys.exit(2)

sample = text[:4096]
try:
    dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
except csv.Error:
    dialect = csv.excel_tab

rows = csv.reader(text.splitlines(), dialect)

def is_number(x):
    try:
        float(x)
        return True
    except Exception:
        return False

def looks_like_header(row):
    if not row or len(row) < 2:
        return True
    a = row[0].strip().lower()
    b = row[1].strip().lower()
    return (a in {"gene","gene_id","id","node","symbol"}) or (b in {"score","logfc","stat","value","hhnet_score"})

total = 0
kept = 0
dropped_not_in_ppi = 0
invalid = 0
header_skipped = False

out_path.parent.mkdir(parents=True, exist_ok=True)

with out_path.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")

    for row in rows:
        if not row:
            continue
        if row[0].lstrip().startswith("#"):
            continue

        if total == 0 and looks_like_header(row):
            header_skipped = True
            total += 1
            continue

        total += 1
        if len(row) < 2:
            invalid += 1
            continue

        gid = row[0].strip()
        val = row[1].strip()
        if not gid or not val or (not is_number(val)):
            invalid += 1
            continue

        if gid not in genes_in_ppi:
            dropped_not_in_ppi += 1
            continue

        w.writerow((gid, "{0}".format(float(val))))
        kept += 1

if header_skipped:
    print("Scores: detected a header row; it was skipped.")
if dropped_not_in_ppi > 0:
    print("Scores: dropped {} row(s) because gene_id was not present in the PPI-derived index.".format(dropped_not_in_ppi))
if invalid > 0:
    print("Scores: skipped {} invalid row(s) (missing columns or non-numeric score).".format(invalid))

if kept == 0:
    print("ERROR: After filtering to PPI genes, 0 score rows remained. Check gene_id naming consistency between the score file and PPI.", file=sys.stderr)
    sys.exit(2)

print("Scores: kept {} row(s) after filtering to PPI genes.".format(kept))
PY
}

compute_file_sha256() {
  local in_path="$1"
  python - "$in_path" <<'PY'
import hashlib
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
h = hashlib.sha256()
with path.open('rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        h.update(chunk)
print(h.hexdigest())
PY
}

stage_cached_similarity_matrix() {
  local edge_list_file="$1"
  local sim_out="$2"
  local beta_out="$3"
  local cache_root="$4"

  local ppi_hash
  ppi_hash="$(compute_file_sha256 "$edge_list_file")"

  local cache_dir="${cache_root}/${ppi_hash}"
  local cache_sim="${cache_dir}/similarity_matrix.h5"
  local cache_beta="${cache_dir}/beta.txt"
  local lock_dir="${cache_dir}.lock"

  mkdir -p "$cache_root"

  if [[ -s "$cache_sim" && -s "$cache_beta" ]]; then
    echo "Using cached similarity matrix for PPI hash: ${ppi_hash}"
    ln -sfn "$cache_sim" "$sim_out"
    ln -sfn "$cache_beta" "$beta_out"
    return 0
  fi

  local have_lock=0
  if mkdir "$lock_dir" 2>/dev/null; then
    have_lock=1
    trap '[[ "$have_lock" -eq 1 ]] && rm -rf "$lock_dir"' RETURN
  else
    echo "Waiting for existing similarity-matrix cache build to finish for PPI hash: ${ppi_hash}"
    while [[ -d "$lock_dir" ]]; do
      sleep 5
    done

    if [[ -s "$cache_sim" && -s "$cache_beta" ]]; then
      echo "Using cached similarity matrix for PPI hash: ${ppi_hash}"
      ln -sfn "$cache_sim" "$sim_out"
      ln -sfn "$cache_beta" "$beta_out"
      return 0
    fi

    if mkdir "$lock_dir" 2>/dev/null; then
      have_lock=1
      trap '[[ "$have_lock" -eq 1 ]] && rm -rf "$lock_dir"' RETURN
    else
      echo "ERROR: Failed to acquire similarity-matrix cache lock for PPI hash: ${ppi_hash}" >&2
      exit 2
    fi
  fi

  mkdir -p "$cache_dir"

  echo "Constructing similarity matrix (cache miss) for PPI hash: ${ppi_hash}"
  python src/construct_similarity_matrix.py \
    -i   "$edge_list_file" \
    -o   "$cache_sim" \
    -bof "$cache_beta"

  [[ -s "$cache_sim" ]]  || { echo "ERROR: Cached similarity matrix was not created: $cache_sim" >&2; exit 2; }
  [[ -s "$cache_beta" ]] || { echo "ERROR: Cached beta file was not created: $cache_beta" >&2; exit 2; }

  ln -sfn "$cache_sim" "$sim_out"
  ln -sfn "$cache_beta" "$beta_out"
}

# ---------------------------------------------------------------------------
# Intermediate-file cache: keyed on SHA-256(indexed_edge_list + scores_0.tsv)
# ---------------------------------------------------------------------------
# Returns 0 (cache hit) or 1 (cache miss).
# On a hit: symlinks every file from the cache dir into INTERMEDIATE_DIR/<NETWORK>_<SCORE>/
# On a miss: after the caller computes the intermediates it should call
#            store_cached_intermediates to populate the cache.
check_cached_intermediates() {
  local edge_list_file="$1"
  local scores0_file="$2"
  local intermediate_dir="$3"   # e.g. .../intermediate/<NETWORK>_<SCORE>
  local cache_root="$4"

  local ppi_hash score_hash cache_key cache_dir stamp
  ppi_hash="$(compute_file_sha256 "$edge_list_file")"
  score_hash="$(compute_file_sha256 "$scores0_file")"
  cache_key="${ppi_hash}_${score_hash}"
  cache_dir="${cache_root}/${cache_key}"
  stamp="${cache_dir}/.complete"

  # Export so store_cached_intermediates can use the same values
  INTERMEDIATE_CACHE_DIR_ACTIVE="$cache_root"
  INTERMEDIATE_CACHE_KEY="$cache_key"

  if [[ -f "$stamp" && -d "$cache_dir" ]]; then
    echo "[Intermediate cache] HIT for key: ${cache_key}"
    echo "[Intermediate cache] Linking cached files into: ${intermediate_dir}"
    for cached_file in "${cache_dir}"/*; do
      [[ -f "$cached_file" ]] || continue
      local fname
      fname="$(basename "$cached_file")"
      ln -sfn "$cached_file" "${intermediate_dir}/${fname}"
    done
    return 0
  fi

  echo "[Intermediate cache] MISS for key: ${cache_key} — will compute intermediates"
  return 1
}

store_cached_intermediates() {
  local intermediate_dir="$1"   # source: .../intermediate/<NETWORK>_<SCORE>
  local num_permutations="$2"
  local network="$3"
  local score="$4"

  [[ -n "${INTERMEDIATE_CACHE_DIR_ACTIVE:-}" ]] || return 0
  [[ -n "${INTERMEDIATE_CACHE_KEY:-}" ]] || return 0

  local cache_dir="${INTERMEDIATE_CACHE_DIR_ACTIVE}/${INTERMEDIATE_CACHE_KEY}"
  local stamp="${cache_dir}/.complete"

  mkdir -p "$cache_dir"

  echo "[Intermediate cache] Storing intermediates to: ${cache_dir}"

  # Copy score_bins
  [[ -f "${intermediate_dir}/score_bins.tsv" ]] && \
    cp "${intermediate_dir}/score_bins.tsv" "${cache_dir}/score_bins.tsv"

  # Copy permuted scores (1..N)
  for i in $(seq "$num_permutations"); do
    local f="${intermediate_dir}/scores_${i}.tsv"
    [[ -f "$f" ]] && cp "$f" "${cache_dir}/scores_${i}.tsv"
  done

  # Copy hierarchy files (0..N)
  for i in $(seq 0 "$num_permutations"); do
    local elf="${intermediate_dir}/hierarchy_edge_list_${i}.tsv"
    local igf="${intermediate_dir}/hierarchy_index_gene_${i}.tsv"
    [[ -f "$elf" ]] && cp "$elf" "${cache_dir}/hierarchy_edge_list_${i}.tsv"
    [[ -f "$igf" ]] && cp "$igf" "${cache_dir}/hierarchy_index_gene_${i}.tsv"
  done

  touch "$stamp"
  echo "[Intermediate cache] Store complete — stamp written"
}

echo "Preparing PPI: generating index_gene.tsv + indexed edge_list.tsv (and removing duplicates / extra columns as needed)..."
prepare_ppi_indexed "$PPI" "$EDGE_LIST_FILE" "$INDEX_GENE_FILE" "$GENES_FILE"

echo "Preparing scores: converting score to abs(score) and filtering to PPI genes..."
normalize_scores_filter_to_ppi "$SCORES" "$GENES_FILE" "$SCORES0_FILE"

pushd "$HHNET_DIR" >/dev/null

if [[ "$COMPILE_FORTRAN" != "never" ]]; then
  if [[ "$COMPILE_FORTRAN" == "always" ]] || ! compgen -G "src/fortran_module*.so" > /dev/null; then
    if command -v f2py >/dev/null 2>&1; then
      echo "Compiling Fortran module..."
      if ! (cd src && f2py -c fortran_module.f95 -m fortran_module > /dev/null); then
        if [[ "$COMPILE_FORTRAN" == "always" ]]; then
          echo "ERROR: Fortran compilation failed and --compile_fortran=always was requested." >&2
          exit 2
        else
          echo "WARN: Fortran compilation failed; continuing with Python fallback." >&2
        fi
      fi
    elif [[ "$COMPILE_FORTRAN" == "always" ]]; then
      echo "ERROR: --compile_fortran always requested, but f2py is not available" >&2
      exit 2
    else
      echo "WARN: f2py not available; continuing with Python fallback" >&2
    fi
  fi
fi

echo "Preparing similarity matrix cache..."
stage_cached_similarity_matrix \
  "$EDGE_LIST_FILE" \
  "$SIMILARITY_MATRIX_FILE" \
  "$BETA_FILE" \
  "$SIMILARITY_CACHE_ROOT"

_INTERMEDIATE_CACHE_HIT=0
if [[ -n "$INTERMEDIATE_CACHE_DIR" ]]; then
  if check_cached_intermediates \
      "$EDGE_LIST_FILE" \
      "$SCORES0_FILE" \
      "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}" \
      "$INTERMEDIATE_CACHE_DIR"; then
    _INTERMEDIATE_CACHE_HIT=1
  fi
fi

if [[ "$_INTERMEDIATE_CACHE_HIT" -eq 0 ]]; then
  echo "Finding permutation bins..."
  python src/find_permutation_bins.py \
    -gsf "$SCORES0_FILE" \
    -igf "$INDEX_GENE_FILE" \
    -elf "$EDGE_LIST_FILE" \
    -ms  1000 \
    -o   "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/score_bins.tsv"

  echo "Permuting scores (${NUM_PERMUTATIONS})..."
  parallel -u -j "$NUM_CORES" --bar \
    python src/permute_scores.py \
      -i  "$SCORES0_FILE" \
      -bf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/score_bins.tsv" \
      -s  {} \
      -o  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_{}.tsv" \
    ::: $(seq "$NUM_PERMUTATIONS")

  echo "Constructing hierarchies..."
  parallel -u -j "$NUM_CORES" --bar \
    python src/construct_hierarchy.py \
      -smf  "$SIMILARITY_MATRIX_FILE" \
      -igf  "$INDEX_GENE_FILE" \
      -gsf  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_{}.tsv" \
      -helf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_edge_list_{}.tsv" \
      -higf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_index_gene_{}.tsv" \
    ::: $(seq 0 "$NUM_PERMUTATIONS")

  if [[ -n "$INTERMEDIATE_CACHE_DIR" ]]; then
    store_cached_intermediates \
      "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}" \
      "$NUM_PERMUTATIONS" \
      "$NETWORK" \
      "$SCORE"
  fi
fi

echo "Processing hierarchies..."
python src/process_hierarchies.py \
  -oelf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_edge_list_0.tsv" \
  -oigf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_index_gene_0.tsv" \
  -pelf $(for i in $(seq "$NUM_PERMUTATIONS"); do echo "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_edge_list_${i}.tsv"; done) \
  -pigf $(for i in $(seq "$NUM_PERMUTATIONS"); do echo "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_index_gene_${i}.tsv"; done) \
  -lsb  10 \
  -cf   "${RESULTS_DIR}/clusters_${NETWORK}_${SCORE}.tsv" \
  -pl   "$NETWORK" "$SCORE" \
  -pf   "${RESULTS_DIR}/sizes_${NETWORK}_${SCORE}.pdf" \
  -nc   "$NUM_CORES"

popd >/dev/null

CLUSTERS_CANON="${RESULTS_DIR}/clusters_${NETWORK}_${SCORE}.tsv"
[[ -s "$CLUSTERS_CANON" ]] || { echo "ERROR: clusters file not created: $CLUSTERS_CANON" >&2; exit 2; }

cp "$CLUSTERS_CANON" "${OUTDIR}/clusters_${NETWORK}_${SCORE}.tsv"

echo "HHNet complete"
echo "  Canonical: $CLUSTERS_CANON"
echo "  Staged:    ${OUTDIR}/clusters_${NETWORK}_${SCORE}.tsv"