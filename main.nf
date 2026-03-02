nextflow.enable.dsl=2

include { HHNET }   from './modules/local/hhnet'
include { RANKING } from './modules/local/ranking'

workflow {

  // DE results (input) should be a file path, e.g. CSV/TSV
  def de_ch  = Channel.fromPath(params.de_results, checkIfExists: true)

  // PPI defaults to params.ppi_network and can be overridden with --ppi_network
  def ppi_ch = Channel.fromPath(params.ppi_network, checkIfExists: true)

  hhnet_out = HHNET(de_ch, ppi_ch)

  // Ranking consumes HHnet outputs + your annotation/druggability/ML/citation tables
  ranked_out = RANKING(
    hhnet_out.metrics,
    Channel.fromPath(params.druggability, checkIfExists: true),
    Channel.fromPath(params.ml_scores,   checkIfExists: true),
    Channel.fromPath(params.citations,   checkIfExists: true)
  )
}
