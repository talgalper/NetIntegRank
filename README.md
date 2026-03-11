# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

NetIntegRank is a modular Nextflow pipeline for identifying dysregulated subnetworks from a user-supplied gene-level score file and prioritising candidate targets by integrating network topology, structural druggability, and optional features.

This repository reflects the workflow developed during my PhD. The broader project data and related analyses can be found in the companion repository:  
`https://github.com/talgalper/PhD-MOC.git`

---

## Overview

NetIntegRank currently runs three main stages:

1. **HHNet**
   - infers dysregulated clusters from a gene-level score file and a protein-protein interaction (PPI) network

2. **HHNet post-processing**
   - reconstructs cluster subnetworks
   - creates neighbour-enriched subnetworks
   - exports network artefacts
   - calculates node-level network metrics

3. **Ranking**
   - integrates network metrics with druggability scores
   - appends ML scores and citation counts
   - outputs ranked, incomplete, and filtered result tables

By default, **ML scores and citation counts are appended to the final outputs but are not used in the rank calculation unless explicitly included in `--ranking_features`**.

---

## Current execution model

The pipeline is currently configured primarily for **Conda-based execution**.

`nextflow.config` also contains `docker` and `apptainer` profiles, but this repository is most clearly set up and documented for `-profile conda`. Unless you have also configured suitable container images separately, Conda should be treated as the supported path.

Provided Conda environments:

- `envs/hhnet.yaml` — Python environment for HHNet
- `envs/r.yaml` — R environment for HHNet post-processing and ranking

---

## Requirements

### Core software

- **Java 17 or newer**
- **Nextflow**
- **Conda / Miniconda**
- **nf-test** (for testing)

### Platform

NetIntegRank is intended for POSIX-style environments such as:

- Linux
- macOS
- Windows via WSL

---

## Installation

### 1) Clone the repository

```bash
git clone https://github.com/talgalper/netintegrank.git
cd netintegrank
````

### 2) Install Java 17+

Ensure a Java 17+ runtime is available:

```bash
java -version
```

If the reported version is older than 17, install a newer JDK before continuing.

### 3) Install Nextflow

A simple installation approach is:

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
chmod +x ~/bin/nextflow
```

Make sure `~/bin` is on your `PATH`, then verify:

```bash
nextflow -version
```

### 4) Install Miniconda

Install Miniconda (if you dont already have it) using the official installer for your platform, then initialise your shell:

```bash
conda init
```

Restart your shell and verify:

```bash
conda --version
```

### 5) Install nf-test

```bash
curl -fsSL https://get.nf-test.com | bash
mv nf-test ~/bin/
chmod +x ~/bin/nf-test
```

Verify installation:

```bash
nf-test version
```

### 6) Ensure scripts in `bin/` are executable

This repository relies on helper scripts in `bin/`. In particular, `run_hhnet.sh` is called directly from the workflow, so missing executable permissions can cause task failures.

A safe approach is to mark all scripts in `bin/` as executable:

```bash
find bin -type f \( -name "*.sh" -o -name "*.R" -o -name "*.py" \) -exec chmod +x {} +
```

At minimum, ensure:

```bash
chmod +x bin/run_hhnet.sh
```

---

## Running the pipeline
> **Note -**
> If your run is interupted, include `-resume` in the command and it will continue from the last cahced save:
> `nextflow run main.nf -profile conda -resume \`
>   ...

> If the run is repeatedly being stopped when trying to annotate IDs, try again later or restart the working environment. Sometimes it's just issues on the Ensembl server side.


### Basic example

```bash
nextflow run main.nf -profile conda \
  --outdir results \
  --scores path/to/scores.tsv \
  --ppi_network path/to/ppi.tsv \
  --druggability path/to/druggability.tsv \
  --ml_scores path/to/ml_scores.tsv \
  --citations path/to/citations.tsv
```

### Example starter command

A more explicit example (using DE scores) is:

```bash
nextflow run main.nf -profile conda \
  --outdir results_example \
  --scores example/logFC_scores_abs.tsv \
  --ppi_network assets/STRING_physical_ENSG.csv \
  --druggability assets/druggability_scores_annot2.0.csv \
  --ml_scores assets/ML_data.csv \
  --citations assets/PubTator3_counts.csv \
  --hhnet_run_id example_run \
  --hhnet_num_permutations 100 \
  --hhnet_min_cluster_size 2 \
  --hhnet_n_clusters 0 \
  --ranking_network neighbours \
  --ranking_features "degree,betweenness,closeness,eigen_centrality,page_rank,highest_score"
```

---

## Required inputs

The pipeline currently expects the following required inputs:

* `--scores`
* `--ppi_network`
* `--druggability`
* `--ml_scores`
* `--citations`

Optional inputs include:

* `--gene_map`
* `--id_annot_cache`
* `--hhnet_run_id`
* `--hhnet_dir`

---

## Input formats

### 1) Gene-level score file (`--scores`)

Expected as a two-column table:

* **column 1:** `gene_id`
* **column 2:** `score`

Notes:

* `gene_id` must match the identifier namespace used in the PPI network
* a header row is allowed
* the second column must be numeric
* the pipeline does **not** convert signed values to absolute values internally
* users should provide scores in the exact form they want HHNet to use
* in the testing setup and in my PhD thesis, the supplied scores were typically **absolute differential expression logFC values**

---

### 2) PPI edge list (`--ppi_network`)

Expected format:

* **column 1:** gene/protein ID
* **column 2:** gene/protein ID

Optional:

* header row
* extra columns (these are ignored by the HHNet wrapper if only the first two columns are required)

The HHNet wrapper is expected to:

* remove self-loops
* collapse duplicated undirected edges
* generate internal indexed network files for HHNet

> **Similarity matrix caching-**
> HHNet can construct a `similarity_matrix.h5` file from the processed PPI network. For large networks this can take substantial time. To reduce repeat work, the wrapper caches the similarity matrix and associated `beta.txt` using the processed/canonicalised PPI network as the cache key (i.e. same PPI file = same similarity matrix). 

---

### 3) Druggability (`--druggability`)

Required columns:

* `uniprot_gn_id`
* `highest_score`

These scores represent structural druggability estimates derived from pocket-detection workflows.

---

### 4) ML scores (`--ml_scores`)

Required columns:

* `uniprot_gn_id`
* `Prediction_Score_rf`

These scores are appended to the final output and are **not used for ranking unless included in `--ranking_features`**.

---

### 5) Citations (`--citations`)

Required columns:

* `external_gene_name`
* `counts`

These counts are appended to the final output and are **not used for ranking unless included in `--ranking_features`**.

---

### 6) Optional gene mapping (`--gene_map`)

Optional TSV/CSV used to map Ensembl IDs to UniProt IDs and gene symbols, reducing or avoiding repeated annotation lookups (via biomaRt).

Recommended columns:

* `ensembl_gene_id`
* `uniprot_gn_id`

Useful additional column:

* `external_gene_name`

If `--gene_map` is not provided, the ranking stage attempts to reuse an annotation cache before falling back to online annotation lookup.

---

## Key parameters

### General

* `--outdir`
  Output directory.
  Default: `results`

---

### HHNet

* `--hhnet_num_permutations`
  Default: `100`

* `--hhnet_network_name`
  Default: `STRING`

* `--hhnet_score_name`
  Default: `DE`

* `--hhnet_compile_fortran`
  Default: `auto`

* `--hhnet_run_id`
  Default: `null`

  Behaviour:

  * if set, HHNet writes run-specific outputs using that ID
  * if not set, behaviour depends on the HHNet wrapper configuration

* `--hhnet_dir`
  Optional explicit path to an HHNet installation/directory

---

### HHNet post-processing

* `--hhnet_min_cluster_size`
  Default: `2`

* `--hhnet_n_clusters`
  Default: `0`
  Meaning: include all clusters meeting the minimum size threshold

* `--hhnet_seed`
  Default: `1234`

* `--hhnet_node_id_type`
  Default: `ensembl_gene_id`

* `--hhnet_plot_max_nodes`
  Default: `5000`

* `--hhnet_extra_colours`
  Default: `""`

  Optional comma-separated colours for clusters beyond the default palette.

---

### Ranking

* `--ranking_network`
  Default: `neighbours`

  Options:

  * `neighbours` — use neighbour-enriched network metrics
  * `subnet` — use cluster-only subnet metrics

* `--ranking_features`
  Default:
  `degree,betweenness,closeness,eigen_centrality,page_rank,highest_score`

  Comma-separated numeric columns used in the averaged ranking.

* `--negative_features`
  Optional.
  Comma-separated features that should be treated as negatively oriented during ranking.

* `--ranking_step`
  Optional.
  Passed through to the ranking script when provided.

* `--gene_map`
  Optional mapping table used instead of online lookup where possible.

* `--id_annot_cache`
  Optional explicit cache path for Ensembl ↔ UniProt ↔ gene symbol annotation.

* `--write_rds`
  Optional.
  Controls whether RDS outputs are written. In practice, the ranking module currently behaves as though this is `true` unless explicitly disabled.

---

## Outputs

### 1) HHNet outputs

Published to:

```text
${outdir}/hhnet/
```

Main output:

```text
clusters_<network_name>_<score_name>.tsv
```

For example:

```text
results/hhnet/clusters_STRING_DE.tsv
```

---

### 2) HHNet post-processing outputs

Published directly under:

```text
${outdir}/hhnet_processing/
```

Key outputs include:

* `hhnet_subnet_metrics.tsv`
* `hhnet_neighbour_metrics.tsv`
* network artefacts generated by the post-processing script
* `id_annot_cache.tsv`

This is a change from earlier README text that referred to `hhnet_post/`.

---

### 3) Ranking outputs

Published to:

```text
${outdir}/ranking/
```

Main outputs:

* `final_ranked.tsv`
* `final_ranked.rds` *(optional if RDS writing is enabled)*
* `final_ranked_incomplete.tsv`
* `final_ranked_incomplete.rds` *(optional if RDS writing is enabled)*
* `final_ranked_missing_external_gene_name.csv`
* `id_annot_cache.tsv` *(optional; copied out if a cache file was written/reused)*

### Output interpretation

* `final_ranked.tsv`
  Final ranked table for complete cases included in ranking.

* `final_ranked_incomplete.tsv`
  Rows retained for reporting but excluded from complete ranking because one or more required ranking fields were missing.

* `final_ranked_missing_external_gene_name.csv`
  Rows removed because `external_gene_name` could not be resolved.

---

## Testing

This repository uses **nf-test**.

### Full pipeline smoke test

If your repository includes `tests/main.nf.test`, run:

```bash
nf-test test tests/main.nf.test --profile conda --verbose
```

This is useful for checking that the main workflow can execute end-to-end in the expected testing configuration.

### Module-level tests

You can also run individual module tests, for example:

```bash
nf-test test tests/modules/local/hhnet.nf.test --profile conda --verbose
nf-test test tests/modules/local/hhnet_post.nf.test --profile conda --verbose
nf-test test tests/modules/local/ranking.nf.test --profile conda --verbose
```

### Running all tests

If your repository is organised with the standard `tests/` layout, you can run all tests with:

```bash
nf-test test tests --profile conda --verbose
```

### Notes on testing

* the tests assume the Conda profile is available
* the tests may rely on files under `tests/fixtures/`, `tests/results/`, and `assets/`
* failures involving `command not found` or `permission denied` often indicate that scripts in `bin/` are not executable

---

## Practical notes

### Annotation cache behaviour

Two different stages may write annotation cache files:

1. **HHNet post-processing** writes:

   * `hhnet_processing/id_annot_cache.tsv`

2. **Ranking** writes or reuses:

   * `${outdir}/ranking/id_annot_cache.tsv` by default
   * or a user-specified path via `--id_annot_cache`

This is expected and can be useful when re-running the pipeline.

### Internet access

If `--gene_map` is not provided and no suitable annotation cache exists, the ranking stage may need internet access for annotation retrieval.

### Reproducibility

For more reproducible behaviour across runs, it is preferable to provide one or both of:

* `--gene_map`
* `--id_annot_cache`

### Similarity-matrix cache

Similarity matrices generated by HHNet are cached under the HHNet directory. Reusing the same processed PPI network should speed up future runs, provided the cache location is writable.

---

## Acknowledgement

NetIntegRank is based on work developed during my PhD, with a focus on dysregulated subnetwork discovery and integrated evidence ranking for target prioritisation.


