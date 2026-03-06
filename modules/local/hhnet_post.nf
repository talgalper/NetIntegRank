process HHNET_POST {
  tag "hhnet_post"
  publishDir "${params.outdir}/hhnet_post", mode: 'copy', overwrite: true

  input:
    path hhnet_clusters
    path ppi_edge_list

  output:
    path "hhnet_processing/hhnet_subnet_metrics.tsv", emit: subnet_metrics
    path "hhnet_processing/hhnet_neighbour_metrics.tsv", emit: neighbour_metrics
    path "hhnet_processing", emit: artefacts_dir

  script:
    """
    mkdir -p hhnet_processing

    Rscript process_hhnet_results.R \
      --hhnet_clusters ${hhnet_clusters} \
      --ppi ${ppi_edge_list} \
      --outdir hhnet_processing \
      --min_cluster_size ${params.hhnet_min_cluster_size ?: 2} \
      --n_clusters ${params.hhnet_n_clusters ?: 0} \
      --seed ${params.hhnet_seed ?: 1234} \
      --node_id_type ${params.hhnet_node_id_type ?: 'ensembl_gene_id'} \
      --id_annot_cache hhnet_processing/id_annot_cache.tsv \
      --extra_colours "${params.hhnet_extra_colours ?: ''}" \
      --plot_max_nodes ${params.hhnet_plot_max_nodes ?: 5000}

    test -s hhnet_processing/hhnet_subnet_metrics.tsv
    test -s hhnet_processing/hhnet_neighbour_metrics.tsv
    """
}