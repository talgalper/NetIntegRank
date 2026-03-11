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
    path "final_ranked_incomplete.tsv", emit: incomplete_tsv
    path "final_ranked_incomplete.rds", optional: true, emit: incomplete_rds
    path "final_ranked_missing_external_gene_name.csv", emit: missing_gene_name_csv
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

    def negFeatArg = (params.negative_features && params.negative_features.toString().trim()) \
      ? "--negative_features '${params.negative_features}'" \
      : ""

    def stepArg = (params.ranking_step != null && params.ranking_step.toString().trim()) \
      ? "--step ${params.ranking_step}" \
      : ""

    def cachePath = (params.id_annot_cache && params.id_annot_cache.toString().trim()) \
      ? params.id_annot_cache.toString() \
      : "${params.outdir}/ranking/id_annot_cache.tsv"

    def writeRdsArg = (params.write_rds != null && params.write_rds.toString().trim())
      ? "--write_rds '${params.write_rds}'"
      : ""

    """
    mkdir -p "\$(dirname '${cachePath}')"

    Rscript ${projectDir}/bin/run_ranking.R \
      --hhnet_metrics ${hhnet_metrics} \
      --druggability ${druggability} \
      --ml_scores ${ml_scores} \
      --citations ${citations} \
      ${geneMapArg} \
      ${featArg} \
      ${negFeatArg} \
      ${stepArg} \
      ${writeRdsArg} \
      --id_annot_cache '${cachePath}' \
      --out_tsv final_ranked.tsv \
      --out_rds final_ranked.rds \
      --out_incomplete_tsv final_ranked_incomplete.tsv \
      --out_incomplete_rds final_ranked_incomplete.rds \
      --out_missing_gene_name_csv final_ranked_missing_external_gene_name.csv
    
    test -s final_ranked.tsv
    test -e final_ranked_incomplete.tsv
    test -e final_ranked_missing_external_gene_name.csv
    
    if [ "${params.write_rds ?: 'true'}" = "true" ] || [ "${params.write_rds ?: 'TRUE'}" = "TRUE" ]; then
      test -s final_ranked.rds
      test -e final_ranked_incomplete.rds
    fi

    if [ -s '${cachePath}' ]; then
      cp '${cachePath}' id_annot_cache.tsv
    fi
    """

  stub:
    """
    cat > final_ranked.tsv <<'EOF'
external_gene_name	ensembl_gene_id	uniprot_gn_id	avg_rank	rank_variance	counts	Prediction_Score_rf
GENE1	ENSG000001	P12345	1.25	0.30	100	0.91
GENE2	ENSG000002	Q99999	2.10	0.55	250	0.73
EOF

    cat > final_ranked_incomplete.tsv <<'EOF'
ensembl_gene_id	external_gene_name	missing_fields
ENSG000003	GENE3	drug_score
EOF

    cat > final_ranked_missing_external_gene_name.csv <<'EOF'
"ensembl_gene_id","external_gene_name","removal_reason"
"ENSG000004","","missing_external_gene_name"
EOF

    touch final_ranked.rds
    touch final_ranked_incomplete.rds
    touch id_annot_cache.tsv
    """
}