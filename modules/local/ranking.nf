// modules/local/ranking.nf
process RANKING {
  tag "ranking"
  publishDir "${params.outdir}/ranking", mode: 'copy', overwrite: true

  input:
    path subnet_metrics
    path neighbour_metrics
    path druggability
    path ml_scores
    path citations
    path gene_map

  output:
    path "final_ranked.tsv"
    path "final_ranked.rds", optional: true

  script:
    def useNeighbour = (params.ranking_network ?: "neighbours") == "neighbours"
    def hhnet_metrics = useNeighbour ? neighbour_metrics : subnet_metrics

    def geneMapArg = gene_map ? "--gene_map ${gene_map}" : ""
    def featArg    = params.ranking_features ? "--ranking_features ${params.ranking_features}" : ""

    """
    Rscript run_ranking.R \
      --hhnet_metrics ${hhnet_metrics} \
      --druggability ${druggability} \
      --ml_scores ${ml_scores} \
      --citations ${citations} \
      ${geneMapArg} \
      ${featArg} \
      --id_annot_cache id_annot_cache.tsv \
      --out_tsv final_ranked.tsv \
      --out_rds final_ranked.rds

    test -s final_ranked.tsv
    """

  stub:
    """
    cat > final_ranked.tsv <<'EOF'
gene_symbol	final_score
GENE1	0.91
GENE2	0.73
EOF

    touch final_ranked.rds
    """
}