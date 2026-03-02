# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

## Overview
NetIntegRank runs two stages:
1) HHnet subnetwork inference from user-provided differential expression results and a PPI network
2) Integrated ranking using network metrics + druggability + ML scores + citation evidence

## Inputs
- DE results (user-supplied)
- PPI network (default provided; optional user override)
- Druggability table
- ML score table
- Citation score table
- Gene ID mapping table (precomputed; no biomaRt dependency)

## HHNet Python environment (pinned)
HHNet depends on older Python-era packages. This repo now includes a pinned Conda environment at:

- `envs/hhnet.yml`

The HHNET process is configured to use that environment via:

- `params.hhnet_conda_env` in `nextflow.config`
- `process.withName: HHNET { conda = params.hhnet_conda_env }`

### Run with Conda
Use the `conda` profile so Nextflow creates/uses the HHNet env automatically:

```bash
nextflow run main.nf -profile conda \
  --de_results path/to/de.tsv \
  --ppi_network path/to/ppi.tsv \
  --druggability path/to/druggability.tsv \
  --ml_scores path/to/ml.tsv \
  --citations path/to/citations.tsv \
  --gene_map path/to/gene_map.tsv
```

### Fortran compile behavior
`bin/run_hhnet.sh` supports:
- `--compile_fortran auto` (default): compile only when module is missing and `f2py` is available
- `--compile_fortran always`: force compile each task invocation
- `--compile_fortran never`: skip compile and use Python fallback

You can control this via `params.hhnet_compile_fortran` in `nextflow.config`.