// Processes for update_data.nf. Each calls the scripts in data_scripts/ (params.data_scripts) (the same scripts and settings that
// built the bundled reference data). Downloads and working files go to <data_outdir>/build/, which the scripts
// reuse when the run is repeated; only the finished tables and their provenance records come back into the
// Nextflow work folder.

// ---------------------------------------------------------------------------------------- STRING network ---
process STRING_NETWORK {
  tag "STRING v${params.string_version}"
  label 'data'
  publishDir path: { "${set_dir}/provenance" }, mode: 'copy', pattern: '*.provenance.txt'

  input:
    val set_dir

  output:
    path 'STRING_physical_ENSG.csv',       emit: ensg
    path 'STRING_physical_geneSymbol.csv', emit: symbol
    path 'STRING.provenance.txt',          emit: provenance

  script:
    def work = "${set_dir}/build/string"
    """
    Rscript '${params.data_scripts}/update_STRING_data.R' \\
      --links-version '${params.string_version}' --info-version '${params.string_info_version}' --out-dir '${work}'
    cp '${work}/STRING_physical_ENSG.csv' '${work}/STRING_physical_geneSymbol.csv' '${work}/STRING.provenance.txt' .
    """

  stub:
    """
    printf '"protein1_ENSG","protein2_ENSG","combined_score"\\n"ENSG01","ENSG02",0.9\\n' > STRING_physical_ENSG.csv
    printf '"protein1_symbol","protein2_symbol","combined_score"\\n"A","B",0.9\\n' > STRING_physical_geneSymbol.csv
    printf 'STRING (stub)\\ndata_version: STRING v${params.string_version}\\n' > STRING.provenance.txt
    """
}


// --------------------------------------------------------------------------------------- citation counts ---
process CITATION_COUNTS {
  tag "PubTator3"
  label 'data'
  label 'ensembl'
  publishDir path: { "${set_dir}/provenance" }, mode: 'copy', pattern: '*.provenance.txt'

  input:
    val set_dir

  output:
    path 'PubTator3_counts.csv',     emit: counts
    path 'citations.provenance.txt', emit: provenance

  script:
    def work = "${set_dir}/build/citations"
    def local_copy = params.gene2pubtator3 ? "--gene2pubtator3 '${file(params.gene2pubtator3)}'" : ""
    """
    Rscript '${params.data_scripts}/update_citation_data.R' --out-dir '${work}' \\
      --out '${work}/PubTator3_counts.csv' ${local_copy}
    cp '${work}/PubTator3_counts.csv' '${work}/citations.provenance.txt' .
    """

  stub:
    """
    printf '"symbol","counts"\\n"TP53",100\\n' > PubTator3_counts.csv
    printf 'citations (stub)\\ndata_version: PubTator3 gene2pubtator3 of stub\\n' > citations.provenance.txt
    """
}


// ------------------------------------------------------------------------------- druggability: structures ---
// A Python for PocketMiner: the pocketminer conda environment's, or --pocketminer_python. Checked by loading
// PocketMiner's packages before any structure is downloaded or scored.
process POCKETMINER_ENV {
  tag "PocketMiner Python"
  label 'pocketminer_env'

  output:
    path 'pocketminer_python.txt', emit: python

  script:
    def python = params.pocketminer_python ?: '$(command -v python)'
    """
    echo "${python}" > pocketminer_python.txt
    export PYTHONDONTWRITEBYTECODE=1   # leave no cache files in the authors' code folder
    "\$(cat pocketminer_python.txt)" '${params.data_scripts}/druggability/pocketminer_src/pocketminer_predict.py' versions
    """

  stub:
    """
    echo python > pocketminer_python.txt
    """
}

// Downloads the AlphaFold DB human proteome and unpacks it (download_alphafold_structures.sh).
process GET_ALPHAFOLD {
  tag "AlphaFold DB ${params.alphafold_version}"
  label 'data'

  input:
    val set_dir

  output:
    tuple val('alphafold'), env('RAW_DIR'), path('alphafold_download.provenance.txt'), emit: structures

  script:
    def out = "${set_dir}/build/druggability/alphafold"
    """
    OUT_DIR='${out}' AF_VERSION='${params.alphafold_version}' bash '${params.data_scripts}/druggability/download_alphafold_structures.sh'
    cp '${out}'/*_v*.provenance.txt alphafold_download.provenance.txt
    RAW_DIR='${out}/structures_raw'
    """

  stub:
    """
    printf 'AlphaFold DB (stub)\\n' > alphafold_download.provenance.txt
    RAW_DIR='${set_dir}/build/druggability/alphafold/structures_raw'
    """
}

// Downloads the SWISS-MODEL Repository's human models and unpacks the one model per protein that is scored
// (download_swissmodel_structures.sh).
process GET_SWISSMODEL {
  tag "SWISS-MODEL"
  label 'data'

  input:
    val set_dir

  output:
    tuple val('swissmodel'), env('RAW_DIR'), path('swissmodel_download.provenance.txt'), emit: structures

  script:
    def out = "${set_dir}/build/druggability/swissmodel"
    """
    OUT_DIR='${out}' bash '${params.data_scripts}/druggability/download_swissmodel_structures.sh'
    cp '${out}'/swissmodel_*.provenance.txt swissmodel_download.provenance.txt
    RAW_DIR='${out}/structures_raw'
    """

  stub:
    """
    printf 'SWISS-MODEL (stub)\\n' > swissmodel_download.provenance.txt
    RAW_DIR='${set_dir}/build/druggability/swissmodel/structures_raw'
    """
}

// Structures supplied with --alphafold_structures or --swissmodel_structures: a record of where they came from.
process USER_STRUCTURES {
  tag "${source}"

  input:
    tuple val(source), val(raw_dir)

  output:
    tuple val(source), val(raw_dir), path("${source}_supplied.provenance.txt"), emit: structures

  script:
    """
    {
      echo "${source} structures supplied with --${source}_structures (not downloaded by update_data.nf)"
      echo "data_version: ${source} structures supplied by the user"
      echo "folder: ${raw_dir}"
      echo "pdb_files: \$(find '${raw_dir}' -maxdepth 1 -name '*.pdb' | wc -l)"
      echo "fingerprint (SHA-256 of the sorted per-file SHA-256 list): \$(cd '${raw_dir}' && find . -maxdepth 1 -name '*.pdb' -printf '%f\\n' | LC_ALL=C sort | xargs sha256sum | sha256sum | cut -d' ' -f1)"
    } > ${source}_supplied.provenance.txt
    """

  stub:
    """
    echo "data_version: ${source} structures supplied by the user (stub)" > ${source}_supplied.provenance.txt
    """
}

// Checks a structure set before scoring (file names, ATOM records, confidence scale): check_structures.sh.
process CHECK_STRUCTURES {
  tag "${source}"
  label 'data_heavy'

  input:
    tuple val(source), val(raw_dir), path(download_record)

  output:
    tuple val(source), val(raw_dir), path("${source}_structures.check.txt"), path(download_record), emit: checked

  script:
    """
    bash '${params.data_scripts}/druggability/check_structures.sh' ${source} '${raw_dir}' -t ${task.cpus} | tee ${source}_structures.check.txt
    """

  stub:
    """
    echo "result: all checks passed (stub)" > ${source}_structures.check.txt
    """
}


// --------------------------------------------------------------------------------- druggability: scoring ---
// Overall confidence of each structure and the prepared structure set (filter_structures.sh).
process FILTER_STRUCTURES {
  tag "${source}"
  label 'data_heavy'

  input:
    tuple val(source), val(raw_dir), path(check), path(download_record)
    val set_dir

  output:
    tuple val(source), path("${source}_struct_score.csv"), emit: filtered

  script:
    def min_conf = source == 'alphafold' ? params.min_struct_conf_af : params.min_struct_conf_sm
    def out = "${set_dir}/build/druggability/${source}/structures_filtered"
    """
    bash '${params.data_scripts}/druggability/filter_structures.sh' ${source} -i '${raw_dir}' -o '${out}' -c ${min_conf} -t ${task.cpus}
    cp '${out}/struct_score.csv' ${source}_struct_score.csv
    """

  stub:
    """
    echo "filename,struct_score" > ${source}_struct_score.csv
    """
}

// fpocket on every prepared structure, then one table per source (Fpocket_master.sh). fpocket comes from the
// conda environment (checked to be release 4.2.3), or is built from source when not running with conda.
process FPOCKET {
  tag "${source}"
  label 'data_heavy'

  input:
    tuple val(source), path(struct_score)
    val set_dir

  output:
    tuple val(source), path("${source}_fpocket_full.csv"), path("${source}_fpocket_small.csv"), emit: tables

  script:
    def base = "${set_dir}/build/druggability/${source}"
    """
    if [ -n "\${CONDA_PREFIX:-}" ] && [ -x "\$CONDA_PREFIX/bin/fpocket" ]; then export FPOCKET_BIN="\$CONDA_PREFIX/bin/fpocket"; fi
    bash '${params.data_scripts}/druggability/Fpocket_master.sh' ${source} -i '${base}/structures_filtered/pdb' -o '${base}/fpocket_results' \\
      -t ${task.cpus} -m ${params.min_residues} -c ${params.min_conf}
    cp '${base}/fpocket_results/fpocket_druggability_full.csv' ${source}_fpocket_full.csv
    cp '${base}/fpocket_results/small_structures.csv' ${source}_fpocket_small.csv
    """

  stub:
    """
    touch ${source}_fpocket_full.csv ${source}_fpocket_small.csv
    """
}

// PocketMiner on every prepared structure, then one table per source (pocketminer_master.sh).
process POCKETMINER {
  tag "${source}"
  label 'data_heavy'

  input:
    tuple val(source), path(struct_score)
    val set_dir
    path python

  output:
    tuple val(source), path("${source}_pocketminer.csv"), emit: tables

  script:
    def base = "${set_dir}/build/druggability/${source}"
    """
    export POCKETMINER_PYTHON="\$(cat ${python})"
    bash '${params.data_scripts}/druggability/pocketminer_master.sh' ${source} -i '${base}/structures_filtered/pdb' -o '${base}/pocketminer_results' \\
      -t ${task.cpus} -m ${params.min_residues} -w ${params.pocketminer_window} -s ${params.pocketminer_threshold} \\
      -c ${params.min_conf}
    cp '${base}/pocketminer_results/pocketminer_results.csv' ${source}_pocketminer.csv
    """

  stub:
    """
    touch ${source}_pocketminer.csv
    """
}


// ------------------------------------------------------------------ druggability: one score per protein ---
// combine_scores.R: each protein's best confident fpocket and PocketMiner score; drug_score is the higher.
process COMBINE_SCORES {
  tag "combine"
  label 'data'
  publishDir path: { "${set_dir}/build/druggability/combined" }, mode: 'copy'

  input:
    path tables     // the six per-source tables
    val set_dir

  output:
    path 'drug_scores.csv', emit: combined

  script:
    """
    Rscript '${params.data_scripts}/druggability/combine_scores.R' \\
      --af-fpocket alphafold_fpocket_full.csv,alphafold_fpocket_small.csv \\
      --sm-fpocket swissmodel_fpocket_full.csv,swissmodel_fpocket_small.csv \\
      --af-pocketminer alphafold_pocketminer.csv \\
      --sm-pocketminer swissmodel_pocketminer.csv \\
      --min-residues ${params.min_residues} --fpocket-score druggability_conf --pocketminer-score max_hit_conf \\
      --out drug_scores.csv
    """

  stub:
    """
    echo "uniprot_id,drug_score" > drug_scores.csv
    """
}

// annotate_scores.R: Ensembl 116 gene names, descriptions and biotypes.
process ANNOTATE_SCORES {
  tag "Ensembl 116"
  label 'data'
  label 'ensembl'
  publishDir path: { "${set_dir}/build/druggability/annotated" }, mode: 'copy'

  input:
    path combined
    val set_dir

  output:
    path 'drug_scores_annotated.csv', emit: annotated

  script:
    """
    Rscript '${params.data_scripts}/druggability/annotate_scores.R' --in ${combined} --out drug_scores_annotated.csv
    """

  stub:
    """
    echo "uniprot_id,external_gene_name,drug_score" > drug_scores_annotated.csv
    """
}

// package_table.R: final column names and order, then the provenance record of the whole build.
process PACKAGE_TABLE {
  tag "package"
  label 'data'
  publishDir path: { "${set_dir}/provenance" }, mode: 'copy', pattern: '*.provenance.txt'

  input:
    path annotated
    path records    // download records and structure checks
    val set_dir

  output:
    path 'druggability_scores_annot.csv',     emit: table
    path 'druggability.provenance.txt',       emit: provenance

  script:
    def build = "${set_dir}/build/druggability"
    """
    Rscript '${params.data_scripts}/druggability/package_table.R' --in ${annotated} --out druggability_scores_annot.csv
    bash '${params.data_scripts}/druggability/src/build_provenance.sh' '${build}' druggability_scores_annot.csv \\
      ${params.min_residues} ${records} > druggability.provenance.txt
    """

  stub:
    """
    printf '"uniprot_id","external_gene_name","drug_score"\\n"P04637","TP53",0.9\\n' > druggability_scores_annot.csv
    printf 'druggability (stub)\\ndata_version: stub\\n' > druggability.provenance.txt
    """
}

// pack_results.sh: packs the large working folders once the table is made.
process PACK_WORKING_FILES {
  tag "pack"
  label 'data_heavy'

  input:
    path table      // waits for the finished table
    val set_dir

  output:
    path 'packed.txt', emit: record

  script:
    def build = "${set_dir}/build/druggability"
    """
    bash '${params.data_scripts}/druggability/pack_results.sh' '${build}' -t ${task.cpus}
    cp '${build}/packed.txt' packed.txt
    """

  stub:
    """
    echo "packed (stub)" > packed.txt
    """
}


// ---------------------------------------------------------------------------------- the finished data set ---
// Puts the tables in the set folder under the names the pipeline reads, so the folder works as
// --reference_dir. A file not rebuilt in this run is kept if the folder already has it (an earlier
// update_data.nf run into the same folder), and copied from --reference_dir otherwise. MANIFEST.tsv records
// each file's size, SHA-256 and origin.
process ASSEMBLE_SET {
  tag "${set_dir}"
  label 'data'
  publishDir path: { "${set_dir}" }, mode: 'copy', pattern: '*.csv'
  publishDir path: { "${set_dir}/provenance" }, mode: 'copy', pattern: 'MANIFEST.tsv'

  input:
    path built, stageAs: 'built/*'
    val set_dir
    val ref_dir

  output:
    path '*.csv', optional: true, emit: tables
    path 'MANIFEST.tsv',          emit: manifest

  script:
    """
    bash '${params.data_scripts}/assemble_set.sh' '${set_dir}' '${ref_dir}' built
    """

  stub:
    """
    bash '${params.data_scripts}/assemble_set.sh' '${set_dir}' '${ref_dir}' built
    """
}
