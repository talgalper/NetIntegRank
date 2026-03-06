#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf -profile conda \
  --outdir results_example \
  --de_results example/logFC_scores_abs.tsv \
  --ppi_network assets/STRING_physical_ENSG.csv \
  --druggability assets/druggability_scores_annot2.0.csv \
  --ml_scores assets/ML_data.csv \
  --citations assets/PubTator3_counts.csv \
  --hhnet_run_id example_run \
  --hhnet_num_permutations 100 \
  --hhnet_min_cluster_size 2 \
  --hhnet_n_clusters 0 \
  --ranking_network neighbours \
  --ranking_features "degree,betweenness,closeness,eigen_centrality,page_rank,highest_score"