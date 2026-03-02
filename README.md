# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

## Overview
NetIntegRank runs two stages:
1) HHnet subnetwork inference from user-provided differential expression results and a PPI network
2) Integrated ranking using network metrics + druggability + ML scores + citation evidence

## Inputs
- DE results (user-supplied)
- PPI network (`assets/STRING_physical_ENSG.tsv` is used by default; optional user override with `--ppi_network`)
- Druggability table
- ML score table
- Citation score table
- Gene ID mapping table (`--gene_map`, optional)
  - If not provided, ranking falls back to biomaRt annotation (`external_gene_name -> uniprot_gn_id`).
  - biomaRt lookup results are cached for reuse across reruns.

## PPI network format (exact)
The HHnet wrapper expects a **tab-delimited edge list** with:
- **No header row**
- **Column 1:** node ID (integer)
- **Column 2:** node ID (integer)

Example:

```tsv
1	2
1	3
2	3
```

See https://github.com/raphael-group/hierarchical-hotnet.git for more details

## Ranking script behavior
`bin/run_ranking.R` is now argument-driven and writes deterministic outputs:

Required arguments:
- `--hhnet_metrics`
- `--druggability`
- `--ml_scores`
- `--citations`
- `--out_tsv`
- `--out_rds`

Optional arguments:
- `--gene_map` (skip biomaRt when available)
- `--id_annot_cache` (cache file for biomaRt conversions; defaults to `<dirname(out_tsv)>/id_annot_cache.tsv`)

The Nextflow ranking module passes a stable cache path under `${params.outdir}/ranking/id_annot_cache.tsv` so cached conversions can be reused between runs.

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

If you omit `--gene_map`, biomaRt annotation + cache will be used:


### Fortran compile behavior
`bin/run_hhnet.sh` supports:
- `--compile_fortran auto` (default): compile only when module is missing and `f2py` is available
- `--compile_fortran always`: force compile each task invocation
- `--compile_fortran never`: skip compile and use Python fallback

You can control this via `params.hhnet_compile_fortran` in `nextflow.config`.