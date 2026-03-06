#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_hhnet.sh --de <de_scores.tsv> --ppi <ppi_edge_list.tsv> [options]

Required:
  --de                Two columns: gene_id, score (score converted to abs()).
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
  --compile_fortran   Fortran compilation mode: auto|always|never (default: auto)
  --hhnet_dir         Path to hierarchical_hotnet directory (default: <script_dir>/hierarchical_hotnet
                      fallback: <script_dir>/hierarchical-hotnet)
  -h, --help          Show this help message

Outputs:
  Canonical HHNet outputs are written under either:
    (a) <hhnet_dir>/runs/<run_id>/{data,intermediate,results}     (if --run_id provided)
    (b) <hhnet_dir>/{data,intermediate,results}                  (if --run_id omitted)

  Downstream file:
    .../results/clusters_<NETWORK>_<SCORE>.tsv

  Staged copy for Nextflow:
    <outdir>/clusters_<NETWORK>_<SCORE>.tsv
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DE=""
PPI=""
OUTDIR="."
RUN_ID=""

NETWORK="STRING"
SCORE="DE"
NUM_PERMUTATIONS=100
NUM_CORES=1
COMPILE_FORTRAN="auto"

# Default HHNET_DIR resolution (prefer underscore, fallback to hyphen)
HHNET_DIR_DEFAULT_UNDERSCORE="${SCRIPT_DIR}/hierarchical_hotnet"
HHNET_DIR_DEFAULT_HYPHEN="${SCRIPT_DIR}/hierarchical-hotnet"
HHNET_DIR="$HHNET_DIR_DEFAULT_UNDERSCORE"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --de) DE="$2"; shift 2 ;;
    --ppi) PPI="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --run_id) RUN_ID="$2"; shift 2 ;;
    --network_name) NETWORK="$2"; shift 2 ;;
    --score_name) SCORE="$2"; shift 2 ;;
    --num_permutations) NUM_PERMUTATIONS="$2"; shift 2 ;;
    --num_cores) NUM_CORES="$2"; shift 2 ;;
    --compile_fortran) COMPILE_FORTRAN="$2"; shift 2 ;;
    --hhnet_dir) HHNET_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$DE" ]] || { echo "ERROR: --de is required" >&2; exit 2; }
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
SCORES0_FILE="${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv"  # gene_id, abs(score)

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
    print(f"ERROR: PPI file not found: {in_path}", file=sys.stderr)
    sys.exit(2)

text = in_path.read_text(errors="replace")
if not text.strip():
    print(f"ERROR: PPI file is empty: {in_path}", file=sys.stderr)
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

    # header skip: only for the first real line
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

    # undirected canonicalisation (AB == BA)
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

# deterministic index: alphabetical order, 1..N
gene_to_idx = {g: i+1 for i, g in enumerate(genes)}

# IMPORTANT FORMAT: index -> gene (col1=index, col2=gene)
with index_out.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    for g in genes:
        w.writerow((gene_to_idx[g], g))

with genes_out.open("w") as f:
    for g in genes:
        f.write(g + "\n")

# Indexed edge list, undirected, written as i<=j
with edge_out.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    for a, b in unique_edges:
        ia = gene_to_idx[a]
        ib = gene_to_idx[b]
        if ia <= ib:
            w.writerow((ia, ib))
        else:
            w.writerow((ib, ia))

# User-facing messages (requested)
if header_skipped:
    print("PPI: detected a header row; it was skipped.")
if extra_col_rows > 0:
    print(f"PPI: detected {extra_col_rows} row(s) with extra columns (e.g. interaction scores). Extra columns will be removed (keeping columns 1-2).")
if reverse_dup_rows > 0:
    print(f"PPI: detected {reverse_dup_rows} reverse-duplicate edge row(s) (A B and B A). Reverse duplicates will be removed (keeping one per undirected pair).")
if directed_dup_rows > 0:
    print(f"PPI: detected {directed_dup_rows} exact duplicate directed edge row(s); duplicates will be removed.")
if undirected_dup_rows > 0:
    print(f"PPI: removed {undirected_dup_rows} duplicate undirected edge row(s) in total (includes reverse and exact duplicates).")
if self_loop_rows > 0:
    print(f"PPI: detected {self_loop_rows} self-loop row(s) (A A); these were removed.")

print(f"PPI: kept {kept_rows} unique undirected edges across {len(genes)} genes.")
PY
}

normalize_de_abs_filter_to_ppi() {
  local de_path="$1"
  local genes_file="$2"
  local out_path="$3"

  python - "$de_path" "$genes_file" "$out_path" <<'PY'
import csv
import pathlib
import sys

de_path = pathlib.Path(sys.argv[1])
genes_file = pathlib.Path(sys.argv[2])
out_path = pathlib.Path(sys.argv[3])

if not de_path.exists():
    print(f"ERROR: DE file not found: {de_path}", file=sys.stderr)
    sys.exit(2)
if not genes_file.exists():
    print(f"ERROR: genes_in_ppi file not found: {genes_file}", file=sys.stderr)
    sys.exit(2)

genes_in_ppi = set(g.strip() for g in genes_file.read_text().splitlines() if g.strip())
if not genes_in_ppi:
    print("ERROR: genes_in_ppi is empty; cannot filter DE.", file=sys.stderr)
    sys.exit(2)

text = de_path.read_text(errors="replace")
if not text.strip():
    print(f"ERROR: DE file is empty: {de_path}", file=sys.stderr)
    sys.exit(2)

sample = text[:4096]
try:
    dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
except csv.Error:
    dialect = csv.excel_tab

rows = csv.reader(text.splitlines(), dialect)

def is_number(x: str) -> bool:
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

        w.writerow((gid, f"{abs(float(val))}"))
        kept += 1

if header_skipped:
    print("DE: detected a header row; it was skipped.")
if dropped_not_in_ppi > 0:
    print(f"DE: dropped {dropped_not_in_ppi} row(s) because gene_id was not present in the PPI-derived index.")
if invalid > 0:
    print(f"DE: skipped {invalid} invalid row(s) (missing columns or non-numeric score).")

if kept == 0:
    print("ERROR: After filtering to PPI genes, 0 DE rows remained. Check gene_id naming consistency between DE and PPI.", file=sys.stderr)
    sys.exit(2)

print(f"DE: kept {kept} row(s) after abs(score) + filtering to PPI genes.")
PY
}

echo "Preparing PPI: generating index_gene.tsv + indexed edge_list.tsv (and removing duplicates / extra columns as needed)..."
prepare_ppi_indexed "$PPI" "$EDGE_LIST_FILE" "$INDEX_GENE_FILE" "$GENES_FILE"

echo "Preparing DE: converting score to abs(score) and filtering to PPI genes..."
normalize_de_abs_filter_to_ppi "$DE" "$GENES_FILE" "$SCORES0_FILE"

pushd "$HHNET_DIR" >/dev/null

if [[ "$COMPILE_FORTRAN" != "never" ]]; then
  if [[ "$COMPILE_FORTRAN" == "always" ]] || [[ ! -f "src/fortran_module"*.so ]]; then
    if command -v f2py >/dev/null 2>&1; then
      echo "Compiling Fortran module..."
      (cd src && f2py -c fortran_module.f95 -m fortran_module > /dev/null)
    elif [[ "$COMPILE_FORTRAN" == "always" ]]; then
      echo "ERROR: --compile_fortran always requested, but f2py is not available" >&2
      exit 2
    else
      echo "WARN: f2py not available; continuing with Python fallback" >&2
    fi
  fi
fi

echo "Construct similarity matrix..."
python src/construct_similarity_matrix.py \
  -i   "$EDGE_LIST_FILE" \
  -o   "${INTERMEDIATE_DIR}/${NETWORK}/similarity_matrix.h5" \
  -bof "${INTERMEDIATE_DIR}/${NETWORK}/beta.txt"

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
    -smf  "${INTERMEDIATE_DIR}/${NETWORK}/similarity_matrix.h5" \
    -igf  "$INDEX_GENE_FILE" \
    -gsf  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_{}.tsv" \
    -helf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_edge_list_{}.tsv" \
    -higf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/hierarchy_index_gene_{}.tsv" \
  ::: $(seq 0 "$NUM_PERMUTATIONS")

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