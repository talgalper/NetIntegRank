nextflow.enable.dsl=2

include { HHNET }      from './modules/local/hhnet'
include { HHNET_POST } from './modules/local/hhnet_post'
include { RANKING }    from './modules/local/ranking'

workflow {

  if( !params.scores )        error "Missing required --scores"
  if( !params.ppi_network )   error "Missing required --ppi_network"
  if( !params.druggability )  error "Missing required --druggability"
  if( !params.ml_scores )     error "Missing required --ml_scores"
  if( !params.citations )     error "Missing required --citations"

  def scores_ch = Channel.fromPath(params.scores, checkIfExists: true)
  def ppi_ch    = Channel.fromPath(params.ppi_network, checkIfExists: true)
  def drug_ch   = Channel.fromPath(params.druggability, checkIfExists: true)
  def ml_ch     = Channel.fromPath(params.ml_scores, checkIfExists: true)
  def cit_ch    = Channel.fromPath(params.citations, checkIfExists: true)

  def hhnet_out = HHNET(scores_ch, ppi_ch)

  // 2) Post-process HHNet clusters into networks + metrics
  def post_out = HHNET_POST(hhnet_out.clusters, ppi_ch)

  // 3) Ranking (defaults to neighbour metrics)
  RANKING(
    post_out.subnet_metrics,
    post_out.neighbour_metrics,
    drug_ch,
    ml_ch,
    cit_ch
  )
}