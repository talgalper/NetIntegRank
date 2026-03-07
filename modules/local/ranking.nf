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

    def geneMapPath = params.gene_map ? file(params.gene_map).toString() : null
    def geneMapArg  = geneMapPath ? "--gene_map '${geneMapPath}'" : ""

    def featArg = params.ranking_features
      ? "--ranking_features '${params.ranking_features}'"
      : ""

    def cachePath = params.id_annot_cache
      ? file(params.id_annot_cache).toString()
      : file("${params.outdir}/ranking/id_annot_cache.tsv").toString()

    def cacheDir = new File(cachePath).getParent()

    """
    mkdir -p '${cacheDir}'

    Rscript run_ranking.R \
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
gene_symbol\tfinal_score
GENE1\t0.91
GENE2\t0.73
EOF

    touch final_ranked.rds
    touch id_annot_cache.tsv
    """
}