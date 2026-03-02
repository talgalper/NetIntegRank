process HHNET {
  tag "hhnet"

  publishDir "${params.outdir}/hhnet", mode: 'copy'

  input:
    path de_results
    path ppi_network

  output:
    path "hhnet/metrics.tsv", emit: metrics
    path "hhnet/subnetwork.tsv", optional: true, emit: subnetwork

  script:
    """
    mkdir -p hhnet

    run_hhnet.sh \
      --de ${de_results} \
      --ppi ${ppi_network} \
      --num_cores ${task.cpus} \
      --num_permutations ${params.hhnet_num_permutations} \
      --network_name ${params.hhnet_network_name} \
      --score_name ${params.hhnet_score_name} \
      --outdir hhnet

    test -s hhnet/metrics.tsv
    """
}
