nextflow.enable.dsl=2

include { HHNET }      from './modules/local/hhnet'
include { HHNET_POST } from './modules/local/hhnet_post'
include { RANKING }    from './modules/local/ranking'

workflow {

  if( !params.de_results )   error "Missing required --de_results"
  if( !params.ppi_network )  error "Missing required --ppi_network"
  if( !params.druggability ) error "Missing required --druggability"
  if( !params.ml_scores )    error "Missing required --ml_scores"
  if( !params.citations )    error "Missing required --citations"

  def de_ch  = Channel.fromPath(params.de_results,  checkIfExists: true)
  def ppi_ch = Channel.fromPath(params.ppi_network, checkIfExists: true)

  // ppi_ch is used by HHNET and HHNET_POST, so split it
  ppi_ch.into { ppi_for_hhnet; ppi_for_post }

  def drug_ch = Channel.fromPath(params.druggability, checkIfExists: true)
  def ml_ch   = Channel.fromPath(params.ml_scores,    checkIfExists: true)
  def cit_ch  = Channel.fromPath(params.citations,    checkIfExists: true)

  // IMPORTANT: for an optional path input, emit null (not an empty channel),
  // otherwise downstream processes won't fire.
  def gene_map_ch = (params.gene_map && params.gene_map.toString().trim())
    ? Channel.fromPath(params.gene_map, checkIfExists: true)
    : Channel.value(null)

  // 1) HHNet
  def hhnet_out = HHNET(de_ch, ppi_for_hhnet)

  // 2) Post-process HHNet clusters into networks + metrics
  def post_out = HHNET_POST(hhnet_out.clusters, ppi_for_post)

  // 3) Ranking (defaults to neighbour metrics)
  RANKING(
    post_out.subnet_metrics,
    post_out.neighbour_metrics,
    drug_ch,
    ml_ch,
    cit_ch,
    gene_map_ch
  )
}