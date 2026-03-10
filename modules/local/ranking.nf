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

  output:
    path "final_ranked.tsv", emit: ranked_tsv
    path "final_ranked.rds", optional: true, emit: ranked_rds
    path "id_annot_cache.tsv", optional: true, emit: id_cache

  script:
    def useNeighbour = (params.ranking_network ?: "neighbours") == "neighbours"
    def hhnet_metrics = useNeighbour ? neighbour_metrics : subnet_metrics

    def geneMapArg = (params.gene_map && params.gene_map.toString().trim()) \
      ? "--gene_map '${params.gene_map}'" \
      : ""

    def featArg = (params.ranking_features && params.ranking_features.toString().trim()) \
      ? "--ranking_features '${params.ranking_features}'" \
      : ""

    def cachePath = (params.id_annot_cache && params.id_annot_cache.toString().trim()) \
      ? params.id_annot_cache.toString() \
      : "${params.outdir}/ranking/id_annot_cache.tsv"

    """
    mkdir -p "\$(dirname '${cachePath}')"

    Rscript ${projectDir}/bin/run_ranking.R \
      --hhnet_metrics ${hhnet_metrics} \
      --druggability ${druggability} \
      --ml_scores ${ml_scores} \
      --citations ${citations} \
      ${geneMapArg} \
      ${featArg} \
      --id_annot_cache '${cachePath}' \
      --out_tsv final_ranked.tsv \
      --out_rds final_ranked.rds

    test -s final_ranked.tsv

    if [ -s '${cachePath}' ]; then
      cp '${cachePath}' id_annot_cache.tsv
    fi
    """

  stub:
    """
    cat > final_ranked.tsv <<'EOF'
gene_symbol	final_score
GENE1	0.91
GENE2	0.73
EOF

    touch final_ranked.rds
    touch id_annot_cache.tsv
    """
}