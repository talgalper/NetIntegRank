nextflow.enable.dsl=2

include { HHNET }      from './modules/local/hhnet'
include { HHNET_POST } from './modules/local/hhnet_post'
include { RANKING }    from './modules/local/ranking'

workflow {

  if( !params.scores )        error "Missing required --scores"

  // Reference data: a file given by its own flag, else the file of that name in --reference_dir,
  // which defaults to the bundled files in assets/
  def ref_dir   = params.reference_dir ?: "${projectDir}/assets"
  def ppi_file  = params.ppi_network  ?: "${ref_dir}/STRING_physical_ENSG.csv"
  def drug_file = params.druggability ?: "${ref_dir}/druggability_scores_annot.csv"
  def ml_file   = params.ml_scores    ?: "${ref_dir}/ML_data.csv"
  def cit_file  = params.citations    ?: "${ref_dir}/PubTator3_counts.csv"

  log.info "Reference data: ppi_network  = ${ppi_file}"
  log.info "                druggability = ${drug_file}"
  log.info "                ml_scores    = ${ml_file}"
  log.info "                citations    = ${cit_file}"

  def scores_ch = channel.fromPath(params.scores, checkIfExists: true)
  def ppi_ch    = channel.fromPath(ppi_file, checkIfExists: true)
  def drug_ch   = channel.fromPath(drug_file, checkIfExists: true)
  def ml_ch     = channel.fromPath(ml_file, checkIfExists: true)
  def cit_ch    = channel.fromPath(cit_file, checkIfExists: true)

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