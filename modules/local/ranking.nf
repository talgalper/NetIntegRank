process RANKING {
  tag "ranking"
  publishDir "${params.outdir}/ranking", mode: 'copy'

  input:
    path hhnet_metrics
    path druggability
    path ml_scores
    path citations
    path gene_map

  output:
    path "ranking/final_ranked.tsv"
    path "ranking/final_ranked.rds", optional: true

  script:
    """
    mkdir -p ranking

    Rscript run_ranking.R \
      --hhnet_metrics ${hhnet_metrics} \
      --druggability ${druggability} \
      --ml_scores ${ml_scores} \
      --citations ${citations} \
      --gene_map ${gene_map} \
      --out_tsv ranking/final_ranked.tsv \
      --out_rds ranking/final_ranked.rds

    test -s ranking/final_ranked.tsv
    """
}
