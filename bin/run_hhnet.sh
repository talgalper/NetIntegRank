#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_hhnet.sh --de <de_scores.tsv> --ppi <ppi_edge_list.tsv> --outdir <output_dir> [options]

Required:
  --de                Gene score file for HHNet (typically differential expression scores)
  --ppi               PPI edge list file (HHNet-compatible edge list)
  --outdir            Output directory for Nextflow process outputs

Optional:
  --index_gene        Precomputed index-gene file; if omitted, an identity mapping is generated
  --network_name      Network label (default: STRING)
  --score_name        Score label (default: DE)
  --num_permutations  Number of score permutations (default: 100)
  --num_cores         Number of cores for GNU parallel and HHNet post-processing (default: 1)
  --compile_fortran   Fortran compilation mode: auto|always|never (default: auto)
  --hhnet_dir         Path to hierarchical-hotnet directory (default: script_dir/hierarchical-hotnet)
  -h, --help          Show this help message
USAGE
}

DE=""
PPI=""
OUTDIR=""
INDEX_GENE=""
NETWORK="STRING"
SCORE="DE"
NUM_PERMUTATIONS=100
NUM_CORES=1
COMPILE_FORTRAN="auto"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HHNET_DIR="${SCRIPT_DIR}/hierarchical-hotnet"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --de) DE="$2"; shift 2 ;;
    --ppi) PPI="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --index_gene) INDEX_GENE="$2"; shift 2 ;;
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
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }
[[ -d "$HHNET_DIR/src" ]] || { echo "ERROR: HHNet src not found at: $HHNET_DIR/src" >&2; exit 2; }
command -v python >/dev/null 2>&1 || { echo "ERROR: python not found in PATH" >&2; exit 2; }
command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel is required but not found" >&2; exit 2; }
[[ "$COMPILE_FORTRAN" =~ ^(auto|always|never)$ ]] || { echo "ERROR: --compile_fortran must be one of auto|always|never" >&2; exit 2; }

DATA_DIR="${OUTDIR}/data"
INTERMEDIATE_DIR="${OUTDIR}/intermediate"
RESULTS_DIR="${OUTDIR}/results"

mkdir -p "$DATA_DIR" "$INTERMEDIATE_DIR" "$RESULTS_DIR"
mkdir -p "${INTERMEDIATE_DIR}/${NETWORK}" "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}"

EDGE_LIST_FILE="edge_list.tsv"
INDEX_GENE_FILE="index_gene.tsv"

cp "$PPI" "${DATA_DIR}/${EDGE_LIST_FILE}"
cp "$DE" "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv"

if [[ -n "$INDEX_GENE" ]]; then
  cp "$INDEX_GENE" "${DATA_DIR}/${INDEX_GENE_FILE}"
else
  awk 'NR>1{print $1; print $2}' "${DATA_DIR}/${EDGE_LIST_FILE}" | awk 'NF>0' | sort -u | awk 'BEGIN{OFS="\t"} {print $1,$1}' > "${DATA_DIR}/${INDEX_GENE_FILE}"
fi

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
  -i   "${DATA_DIR}/${EDGE_LIST_FILE}" \
  -o   "${INTERMEDIATE_DIR}/${NETWORK}/similarity_matrix.h5" \
  -bof "${INTERMEDIATE_DIR}/${NETWORK}/beta.txt"

echo "Finding permutation bins..."
python src/find_permutation_bins.py \
  -gsf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv" \
  -igf "${DATA_DIR}/${INDEX_GENE_FILE}" \
  -elf "${DATA_DIR}/${EDGE_LIST_FILE}" \
  -ms  1000 \
  -o   "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/score_bins.tsv"

echo "Permuting scores (${NUM_PERMUTATIONS})..."
parallel -u -j "$NUM_CORES" --bar \
  python src/permute_scores.py \
    -i  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv" \
    -bf "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/score_bins.tsv" \
    -s  {} \
    -o  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_{}.tsv" \
  ::: $(seq "$NUM_PERMUTATIONS")

echo "Constructing hierarchies..."
parallel -u -j "$NUM_CORES" --bar \
  python src/construct_hierarchy.py \
    -smf  "${INTERMEDIATE_DIR}/${NETWORK}/similarity_matrix.h5" \
    -igf  "${DATA_DIR}/${INDEX_GENE_FILE}" \
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

# Nextflow-facing outputs
cp "${RESULTS_DIR}/clusters_${NETWORK}_${SCORE}.tsv" "${OUTDIR}/subnetwork.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1{print "node","hhnet_score"; next} NF>0{print $1,$2}' \
  "${INTERMEDIATE_DIR}/${NETWORK}_${SCORE}/scores_0.tsv" > "${OUTDIR}/metrics.tsv"

test -s "${OUTDIR}/metrics.tsv"
echo "HHNet run complete: ${OUTDIR}"