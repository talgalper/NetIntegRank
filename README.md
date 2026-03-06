# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

This pipeline is based on the work done in my PhD. I wouldnt recomend trawling through the repo but in case you really want to, all the data for the project can be found at: https://github.com/talgalper/PhD-MOC.git

## Overview

NetIntegRank runs three stages:

1) **HHNet**: infer dysregulated clusters from differential expression (DE) scores and a PPI network  
2) **HHNet post-processing**: build (a) subnet and (b) neighbour-enriched subnet, export GraphML + plots, compute network metrics  
3) **Ranking**: integrate network metrics + druggability (and optionally user-added features) into a final ranked table  
   - ML scores and citation counts are appended to the output but are **not used for ranking by default**

## Requirements

- Nextflow (and Java 17+)
- Conda (recommended) OR container runtime (Docker/Apptainer) if you set containers

This repo expects two conda environments (provided):
- `envs/hhnet.yaml` (Python stack for HHNet)
- `envs/r.yaml` (R stack for HHNet post-processing + ranking)

### Run with conda
```bash
nextflow run main.nf -profile conda \
  --de_results path/to/de.tsv \
  --ppi_network path/to/ppi.tsv \
  --druggability path/to/druggability.tsv \
  --ml_scores path/to/ml_scores.tsv \
  --citations path/to/citations.tsv
````

## Input formats

### 1) Differential expression scores (`--de_results`)

Differential expression scores (logFC) should be as absolute values in order to get the best result. However, You will lose the directionality of up and down regulation and should be aware of this when interpreting results. I would not recomend using signed logFC scores unless you absolutely know what you are doing.

Two columns:

* col1: gene_id (must match the ID namespace used in the PPI edge list)
* col2: score (converted to abs(score) inside HHNet wrapper)

Header row is allowed.

### 2) PPI edge list (`--ppi_network`)

The PPI network provided wis a human network extracted from STRINGdb (last updated: Nov 2024).
A custom network can be provided using the following format:

* col1: gene_id
* col2: gene_id

Optional:

* header row (allowed)
* third+ columns (e.g. interaction scores) are detected and dropped (col1/col2 are kept)

The HHNet wrapper will:

* remove self-loops (A A)
* collapse AB/BA duplicates to one undirected edge
* generate an index mapping and indexed edge list internally for HHNet

### 3) Druggability (`--druggability`)

Druggability scores were obtained using a combination of Fpocket and PocketMiner tools on both AlphaFold2 and SWISS-MODEL datasets. The score represents the likeliness that there is a small molecule drug binding pocket.

Required columns:

* `uniprot_gn_id`
* `highest_score`

### 4) ML scores (`--ml_scores`)

The Machine Learning (ML) scores were obtained by running a random forest model on a large set of features extarcted from a previous study as well as some custom features input by me. The paper the model was based off can be found here: https://pubmed.ncbi.nlm.nih.gov/32171238/ 

Required columns:

* `uniprot_gn_id`
* `Prediction_Score_rf`

> Note: ML scores are appended to the final table but are not used for ranking unless included in `--ranking_features`.

### 5) Citations (`--citations`)

Required columns:

* `external_gene_name`
* `counts`

> Note: citation counts are appended to the final table but are not used for ranking unless included in `--ranking_features`.

### 6) Optional gene mapping (`--gene_map`)

Optional TSV/CSV used to map Ensembl → UniProt and avoid biomaRt calls.

Recommended columns:

* `ensembl_gene_id`
* `uniprot_gn_id`

If not provided, `run_ranking.R` will use biomaRt and cache results.

## Key parameters

### HHNet

* `--hhnet_num_permutations` (default: 100)
* `--hhnet_network_name` (default: STRING)
* `--hhnet_score_name` (default: DE)
* `--hhnet_compile_fortran` (default: auto)
* `--hhnet_run_id` (default: null)

  * if set, HHNet writes canonical outputs under `bin/hierarchical_hotnet/runs/<run_id>/...`
  * if not set, HHNet overwrites outputs under `bin/hierarchical_hotnet/...` (wrapper warns and pauses 5 seconds)

### HHNet post-processing

* `--hhnet_min_cluster_size` (default: 2)
* `--hhnet_n_clusters` (default: 0, meaning "all clusters with size >= min_cluster_size")
* `--hhnet_node_id_type` (default: ensembl_gene_id)
* `--hhnet_plot_max_nodes` (default: 5000)
* `--hhnet_extra_colours` (default: "")

  * colours for clusters 11+; clusters 1–10 are fixed to:
    tomato, springgreen, royalblue, maroon1, gold, orchid, cyan, yellowgreen, mediumseagreen, saddlebrown

Outputs are written to:

* `${params.outdir}/hhnet_post/` (GraphML, plots, metrics TSVs)

### Ranking

* `--ranking_network` (default: neighbours)

  * `neighbours` uses the neighbour-enriched network metrics
  * `subnet` uses the cluster-only subnet metrics
* `--ranking_features` (default: `degree,betweenness,closeness,eigen_centrality,page_rank,highest_score`)

  * Comma-separated column names in the HHNet metrics table to include in the averaged rank
  * Users can append custom numeric columns to the metrics table and include them here

Final outputs:

* `${params.outdir}/ranking/final_ranked.tsv`
* `${params.outdir}/ranking/final_ranked.rds`

## Example “starter” run command

See `examples/run_example.sh` for a copy/paste template.

## Notes

* If biomaRt is used, it requires internet access. To avoid this, provide `--gene_map` (or a pre-populated cache).