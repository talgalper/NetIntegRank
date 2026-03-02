#!/usr/bin/env bash
set -euo pipefail

while [[ $# -gt 0 ]]; do
  case "$1" in
    --de) DE="$2"; shift 2 ;;
    --ppi) PPI="$2"; shift 2 ;;
    --num_cores) NUM_CORES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${DE:-}" ]]; then
  echo "Missing required argument: --de" >&2
  exit 2
fi

if [[ -z "${PPI:-}" ]]; then
  echo "Missing required argument: --ppi" >&2
  exit 2
fi

if [[ -z "${OUTDIR:-}" ]]; then
  echo "Missing required argument: --outdir" >&2
  exit 2
fi

if [[ -z "${NUM_CORES:-}" ]]; then
  echo "Missing required argument: --num_cores" >&2
  exit 2
fi

mkdir -p "${OUTDIR}"


# Placeholder: replace with your real HHnet invocation.
# Example outputs expected by Nextflow:
#   ${OUTDIR}/metrics.tsv
#   ${OUTDIR}/subnetwork.tsv (optional)

echo -e "node\tdegree\tbetweenness\tcloseness\teigenvector" > "${OUTDIR}/metrics.tsv"
