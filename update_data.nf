nextflow.enable.dsl=2

// Builds new versions of NetIntegRank's reference data into a dated folder outside assets/:
//
//   nextflow run update_data.nf -profile conda --datasets string,citations
//   nextflow run update_data.nf -profile conda --datasets all --data_outdir reference_data/2027-03
//
// The folder then holds the files under the names the pipeline reads, plus provenance/ and build/, and is used
// with: nextflow run main.nf --reference_dir <folder> ... Nothing in assets/ is changed. See "Updating the
// reference data" in the README.

include { STRING_NETWORK; CITATION_COUNTS; POCKETMINER_ENV; GET_ALPHAFOLD; GET_SWISSMODEL; USER_STRUCTURES;
          CHECK_STRUCTURES; FILTER_STRUCTURES; FPOCKET; POCKETMINER; COMBINE_SCORES; ANNOTATE_SCORES;
          PACKAGE_TABLE; PACK_WORKING_FILES; ASSEMBLE_SET } from './modules/local/update_data'

workflow {

  // ---- which datasets
  def all_datasets = ['string', 'citations', 'druggability']
  if( !params.datasets ) error "Give --datasets: string, citations, druggability (comma-separated), or all"
  def requested = params.datasets.toString().trim() == 'all' ? all_datasets
                : params.datasets.toString().tokenize(',').collect { d -> d.trim() }
  def unknown = requested.findAll { d -> !(d in all_datasets) }
  if( unknown ) error "Unknown --datasets ${unknown.join(', ')}: choose from ${all_datasets.join(', ')}, or all"

  // ---- where the new set goes, and where the files that are not rebuilt come from
  def today   = java.time.LocalDate.now().toString()
  def set_dir = file(params.data_outdir ?: "${launchDir}/reference_data/${today}").toAbsolutePath().normalize()
  def ref_dir = file(params.reference_dir ?: "${projectDir}/assets").toAbsolutePath().normalize()
  def assets  = file("${projectDir}/assets").toAbsolutePath().normalize()
  if( set_dir.startsWith(assets) )
    error "--data_outdir must be outside assets/: update_data.nf never changes the bundled files"
  if( !ref_dir.isDirectory() ) error "--reference_dir ${ref_dir} is not a folder"

  log.info "update_data.nf: building ${requested.join(', ')}"
  log.info "  new set:        ${set_dir}"
  log.info "  other files:    copied from ${ref_dir}"
  log.info "  when finished:  nextflow run main.nf -profile conda --reference_dir ${set_dir} --scores <your scores>"

  // ---- checks before anything starts: druggability needs Linux, disk, CPUs and GNU tools
  if( 'druggability' in requested && !workflow.stubRun ) {
    if( !System.getProperty('os.name').toLowerCase().contains('linux') )
      error "The druggability rebuild runs on Linux only (its scripts use GNU tools); STRING and citations run anywhere"
    def gnu = ['bash', '-c', 'tar --version | grep -q "GNU tar" && find --version | grep -q GNU && command -v sha256sum'].execute()
    gnu.waitFor()
    if( gnu.exitValue() != 0 ) error "The druggability rebuild needs GNU tar, GNU find and sha256sum"
    // the nearest existing folder on the way to set_dir: the disk the new set will be written to
    def existing = (0..<set_dir.nameCount).collect { i -> set_dir.root.resolve(set_dir.subpath(0, set_dir.nameCount - i)) }
      .find { p -> p.exists() } ?: set_dir.root
    def free_gb = (existing.toFile().getUsableSpace() / 1e9) as int
    log.info "  free disk:      ${free_gb} GB at ${existing}"
    if( params.min_free_disk_gb && free_gb < (params.min_free_disk_gb as int) )
      error "The druggability rebuild needs about ${params.min_free_disk_gb} GB free at ${existing} (structures, fpocket and PocketMiner working files); ${free_gb} GB is free. Free some space, choose another --data_outdir, or lower --min_free_disk_gb if part of the build is already there"
    def cpus = params.data_cpus ?: Runtime.runtime.availableProcessors()
    log.info "  CPUs:           ${cpus}"
    if( (cpus as int) < 8 ) log.warn "Only ${cpus} CPUs: the druggability rebuild will take many hours (about 3 h of fpocket and PocketMiner on 16 CPUs)"
  }

  def built = channel.empty()

  // ---- STRING network
  if( 'string' in requested ) {
    STRING_NETWORK(set_dir.toString())
    built = built.mix(STRING_NETWORK.out.ensg, STRING_NETWORK.out.symbol)
  }

  // ---- citation counts
  if( 'citations' in requested ) {
    CITATION_COUNTS(set_dir.toString())
    built = built.mix(CITATION_COUNTS.out.counts)
  }

  // ---- druggability
  if( 'druggability' in requested ) {
    def python = POCKETMINER_ENV().python

    // structures: downloaded, or a folder given for that source
    def supplied = channel.empty()
    def downloaded = channel.empty()
    if( params.alphafold_structures )
      supplied = supplied.mix(channel.of(['alphafold', file(params.alphafold_structures, checkIfExists: true).toString()]))
    else
      downloaded = downloaded.mix(GET_ALPHAFOLD(set_dir.toString()).structures)
    if( params.swissmodel_structures )
      supplied = supplied.mix(channel.of(['swissmodel', file(params.swissmodel_structures, checkIfExists: true).toString()]))
    else
      downloaded = downloaded.mix(GET_SWISSMODEL(set_dir.toString()).structures)
    def structures = downloaded.mix(USER_STRUCTURES(supplied).structures)

    CHECK_STRUCTURES(structures)
    FILTER_STRUCTURES(CHECK_STRUCTURES.out.checked, set_dir.toString())
    FPOCKET(FILTER_STRUCTURES.out.filtered, set_dir.toString())
    POCKETMINER(FILTER_STRUCTURES.out.filtered, set_dir.toString(), python)

    def tables = FPOCKET.out.tables.map { _source, full, small -> [full, small] }
      .mix(POCKETMINER.out.tables.map { _source, pm -> pm })
      .flatten().collect()
    COMBINE_SCORES(tables, set_dir.toString())
    ANNOTATE_SCORES(COMBINE_SCORES.out.combined, set_dir.toString())

    def records = CHECK_STRUCTURES.out.checked.map { _source, _raw, check, download -> [check, download] }.flatten().collect()
    PACKAGE_TABLE(ANNOTATE_SCORES.out.annotated, records, set_dir.toString())
    built = built.mix(PACKAGE_TABLE.out.table)

    if( params.pack_working_files ) PACK_WORKING_FILES(PACKAGE_TABLE.out.table, set_dir.toString())
  }

  // ---- the finished set: tables under the pipeline's names, plus a manifest
  ASSEMBLE_SET(built.collect(), set_dir.toString(), ref_dir.toString())
}
