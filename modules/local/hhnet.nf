// modules/local/hhnet.nf
process HHNET {
  tag "${params.hhnet_network_name}_${params.hhnet_score_name}"
  publishDir "${params.outdir}/hhnet", mode: 'copy', overwrite: true

  input:
    path de_scores
    path ppi_edge_list

  output:
    path "clusters_${params.hhnet_network_name}_${params.hhnet_score_name}.tsv", emit: clusters

  script:
    def hhnetDirArg = params.hhnet_dir ? "--hhnet_dir ${params.hhnet_dir}" : ""
    def runIdArg    = params.hhnet_run_id ? "--run_id ${params.hhnet_run_id}" : ""

    """
    run_hhnet.sh \
      --de ${de_scores} \
      --ppi ${ppi_edge_list} \
      --outdir . \
      ${runIdArg} \
      --network_name ${params.hhnet_network_name} \
      --score_name ${params.hhnet_score_name} \
      --num_permutations ${params.hhnet_num_permutations} \
      --num_cores ${task.cpus} \
      --compile_fortran ${params.hhnet_compile_fortran} \
      ${hhnetDirArg}

    test -s clusters_${params.hhnet_network_name}_${params.hhnet_score_name}.tsv
    """

  stub:
    """
    cat > clusters_${params.hhnet_network_name}_${params.hhnet_score_name}.tsv <<'EOF'
cluster_id	node_id
1	ENSG000001
1	ENSG000002
EOF
    """
}