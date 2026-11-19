# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

NetIntegRank is a modular Nextflow pipeline for identifying dysregulated subnetworks from a user-supplied gene-level score file and prioritising candidate targets by integrating network topology, structural druggability, and optional features.

This repository reflects the workflow developed during my PhD. The broader project data and related analyses can be found in the repository: `https://github.com/talgalper/PhD-MOC.git`.

---

**Current version: 2.0.** Rankings from 2.0 are not byte identical to rankings from the original version, 1.0 (March 2026), and 2.0 ships updated reference data. See [Version history](#version-history) at the end of this README for what changed.

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

By default, **ML scores and citation counts are appended to the final outputs but are not used in the rank calculation**. Citation counts can be added to the ranking through `--ranking_features` (as `counts` or `counts_norm`); ML scores cannot, because they are joined after the ranking is calculated.

---

## Current execution model

The pipeline is currently configured primarily for **Conda-based execution**.

`nextflow.config` also contains `docker` and `apptainer` profiles, but this repository is most clearly set up and documented for `-profile conda`. Unless you have also configured suitable container images separately, Conda should be treated as the supported path. May add proper docker support at a later date.

Provided Conda environments:

- `envs/hhnet.yaml` — Python environment for HHNet
- `envs/r.yaml` — R environment for HHNet post-processing and ranking

---

## Requirements

### Core software

- **Java 17 or newer**
- **Nextflow** (version 2.0 of the pipeline is tested on Nextflow 26.10.5)
- **Conda / Miniconda**
- **nf-test** (for testing; version 2.0 of the pipeline is tested with nf-test 0.9.5)

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
mv nextflow /bin/
chmod +x /bin/nextflow
```

Make sure `/bin` is on your `PATH`, then verify:

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

> Recent Miniconda releases will not build any environment until Anaconda's Terms of Service for its default channels have been accepted: the first run stops with `CondaToSNonInteractiveError`, and the message gives the `conda tos accept` commands. Either accept them, or install [Miniforge](https://conda-forge.org/download/) instead of Miniconda; Miniforge uses only the conda-forge channel and has no such step. The pipeline's environments build with either.

### 5) Install nf-test

```bash
curl -fsSL https://get.nf-test.com | bash
mv nf-test /bin/
chmod +x /bin/nf-test
```

Verify installation:

```bash
nf-test version
```

### 6) Ensure scripts in `bin/` are executable

This repository relies on helper scripts in `bin/`. All three (`run_hhnet.sh`, `process_hhnet_results.R` and `run_ranking.R`) are called directly by name from the workflow, so missing executable permissions cause task failures.

A safe approach is to mark all scripts in `bin/` as executable:

```bash
find bin -type f \( -name "*.sh" -o -name "*.R" -o -name "*.py" \) -exec chmod +x {} +
```

At minimum, ensure:

```bash
chmod +x bin/run_hhnet.sh bin/process_hhnet_results.R bin/run_ranking.R
```

---

## Running the pipeline
> **Note -**
> If your run is interupted, include `-resume` in the command and it will continue from the last cahced save:
> `nextflow run main.nf -profile conda -resume \`
>   ...

> If the run is repeatedly being stopped when trying to annotate IDs, try again later or restart the working environment. Sometimes it's just issues on the Ensembl server side.


### Example run (default params)

It is suggested to defaults unless you are comfortable with this tool which you can run easily by calling the default [parameter file](params/example_default.yaml) which you will need to update at the indicated fields. I would suggested making a copy of this file and changing the name to something more appropriate as well.

```bash
nextflow run main.nf -profile conda -params-file params/example_default.yaml
```

The params can be run in command line as well if you want to do your own scripting:
<details>
  <summary>show example commands</summary>

  ### Basic example

  ```bash
  nextflow run main.nf -profile conda \
    --outdir results \
    --scores path/to/scores.tsv
  ```

  This uses the bundled reference data. To use other reference files, add `--reference_dir` or any of `--ppi_network`, `--druggability`, `--ml_scores` and `--citations` (see [Required inputs](#required-inputs)).

  ### Example starter command

  A more explicit example (using DE scores) is:

  ```bash
  nextflow run main.nf -profile conda \
    --outdir results_example \
    --scores example/logFC_scores_abs.tsv \
    --hhnet_run_id example_run \
    --hhnet_num_permutations 100 \
    --hhnet_min_cluster_size 2 \
    --hhnet_n_clusters 0 \
    --ranking_network neighbours \
    --ranking_features "degree,betweenness,closeness,eigen_centrality,page_rank,drug_score"
  ```
</details>

<br/>

Be sure to test the pipeline before your first run. Jump to the [testing section](#testing) below for more info.

---

## Required inputs

The only required input is:

* `--scores`

The reference data (PPI network, druggability scores, ML scores and citation counts) defaults to the files bundled in `assets/`, described under [Bundled reference data](#bundled-reference-data). To use other files:

* `--reference_dir`: a folder holding all four files under the bundled names (`STRING_physical_ENSG.csv`, `druggability_scores_annot.csv`, `ML_data.csv`, `PubTator3_counts.csv`)
* `--ppi_network`, `--druggability`, `--ml_scores`, `--citations`: one file each, used for that input in place of the file from `--reference_dir` or `assets/`

The reference files a run uses are printed when it starts.

Optional inputs include:

* `--gene_map`
* `--id_annot_cache`
* `--hhnet_run_id`
* `--hhnet_dir`

> To build newer versions of the reference data, see [Updating the reference data](#updating-the-reference-data).

---

## Bundled reference data

The files in `assets/` are the human reference datasets, and the pipeline uses them by default (see [Required inputs](#required-inputs)). Version 2.0 updated the PPI network, the druggability scores and the citation counts. The files that shipped with 1.0 (used in my PhD) are kept in `assets/legacy_data/` (see [Legacy data](#legacy-data)). To build newer versions, see [Updating the reference data](#updating-the-reference-data).

| File | Contents | Source |
|---|---|---|
| `STRING_physical_ENSG.csv`, `STRING_physical_geneSymbol.csv` | Physical protein-protein interactions: 430,638 interactions between 17,828 genes, as Ensembl gene IDs or as gene symbols | STRING v12.5 |
| `druggability_scores_annot.csv` | One structure-based druggability score per protein: 20,516 proteins, 19,433 gene names | AlphaFold DB v6 and SWISS-MODEL Repository structures, scored with fpocket and PocketMiner; gene names from Ensembl 116 |
| `PubTator3_counts.csv` | Number of PubMed articles that mention each gene: 24,507 genes | PubTator3; gene names from Ensembl 116 |
| `ML_data.csv` | Machine-learning target scores (`Prediction_Score_rf`) and the features behind them | Unchanged from 1.0 |

### STRING PPI network

* Built from STRING v12.5 `9606.protein.physical.links.full.v12.5.txt.gz` (physical interactions only), downloaded in November 2026 from https://string-db.org/cgi/download.
* Interactions with a `combined_score` of at least 400 are kept, and the score is stored divided by 1000 (0.4 to 0.999), as in 1.0.
* STRING protein IDs are mapped to Ensembl gene IDs through STRING's protein aliases (source `Ensembl_gene`) and to gene symbols through STRING's protein info file. STRING 12.5 reuses the 12.0 versions of both files.
* Each interaction is listed once.
* The gene-symbol file uses STRING's own gene names, which date from 2023 and differ from current names for some genes (for example `G6PC`, now `G6PC1`). The druggability and citation tables use Ensembl 116 names, so a renamed gene finds no match when the symbol file is the network. The Ensembl ID file avoids this, because the pipeline names its genes from Ensembl 116.

### Druggability scores

`druggability_scores_annot.csv` gives each human protein a structure-based druggability score, `drug_score`, between 0 and 1.

* **Structures.** AlphaFold DB v6 models of the human reference proteome (UP000005640, https://alphafold.ebi.ac.uk), and for each protein the SWISS-MODEL Repository model (https://swissmodel.expasy.org/repository, models of release 2026-09.1) with a sequence identity of at least 60% and the best QMEANDisCo global score. AlphaFold DB provides proteins longer than 2,700 residues as overlapping fragments, and each fragment is scored on its own.
* **Filters.** Models shorter than 80 residues are not scored. Models are not filtered on their overall confidence: the confident-region rule below leaves out their poorly modelled parts instead, so a mostly disordered protein is still scored on its well-modelled regions.
* **Confident regions only.** Per-residue confidence (AlphaFold pLDDT, or the SWISS-MODEL QMEANDisCo local score multiplied by 100, so both are on a 0 to 100 scale) is used to ignore pockets in poorly modelled regions. A pocket or window counts only if its mean confidence is at least 70.
* **`FP_score`:** fpocket 4.2.3 (default settings) druggability of the best pocket whose lining residues pass the confidence rule.
* **`CP_score`:** PocketMiner cryptic-pocket score (Meller et al. 2023, *Nature Communications* 14:1177; code and model from the authors' repository, https://github.com/Mickdub/gvp, commit `187062d`). Each residue's prediction is averaged with the 10 residues on either side, and the highest average among windows that pass the confidence rule is kept.
* **`drug_score`** is the higher of the protein's best `FP_score` and best `CP_score` across its AlphaFold fragments and its SWISS-MODEL model.
* Proteins with no usable score have `drug_score` 0, so their genes are still ranked (928 proteins: 531 with no pocket or window in a confident region, 394 whose only structures are under 80 residues, and 3 that score 0).
* Gene names, descriptions and biotypes come from Ensembl 116, the release the pipeline names its genes with. There is one row per protein and gene name, and the pipeline keeps the highest `drug_score` per gene name.
* Proteins that Ensembl 116 does not link to a gene have no gene name, and their scores cannot reach the pipeline's genes. This includes a few well-known genes whose main UniProt protein Ensembl 116 does not record, for example POLR2A (P24928) and LIN52 (Q52LA3). Such genes have no `drug_score`, so the pipeline lists them in `final_ranked_incomplete.tsv` rather than ranking them.

| Column | Contents |
|---|---|
| `uniprot_id` | UniProt accession |
| `external_gene_name` | Ensembl 116 gene name; empty or `NA` where Ensembl has no match (such rows cannot join to the pipeline's genes) |
| `description`, `gene_biotype` | Ensembl 116 gene description and biotype |
| `FP_score` | Best fpocket druggability among confident pockets; 0 if none |
| `FP_source` | Structure that gave `FP_score`: `AF` (AlphaFold DB) or `SM` (SWISS-MODEL) |
| `FP_struct_score` | Overall confidence of that structure, 0 to 100 |
| `CP_score` | Best PocketMiner score among confident windows; 0 if none |
| `CP_source` | Structure that gave `CP_score`: `AF` or `SM` |
| `drug_score` | The higher of `FP_score` and `CP_score` |
| `drug_score_tool` | Tool that gave `drug_score`: `fpocket` or `PocketMiner` (fpocket on a tie); empty if no score |
| `drug_score_source` | Structure that gave `drug_score`: `AF` or `SM` (AlphaFold on a tie) |

### Citation counts

* `counts` is the number of distinct PubMed articles in which PubTator3 found a mention of the gene. Every article counts, whatever its topic.
* Built from PubTator3's gene annotations (`gene2pubtator3.gz`, downloaded on 13 November 2026 from https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/). Mentions that PubTator links to more than one gene are left out.
* PubTator3 identifies genes by NCBI gene ID. These are named through Ensembl 116, so only human genes are kept, under the same gene names the pipeline uses. Where Ensembl links one NCBI gene to several genes (mostly near-identical paralogs), the gene that carries NCBI's own symbol is used.
* Columns: `symbol`, `counts`.

### ML scores

`ML_data.csv` is unchanged from 1.0. It was built in 2024 from the data available then (including STRING 12.0 network metrics) and has not been rebuilt for 2.0.

### Legacy data

`assets/legacy_data/` holds the reference files that 1.0 shipped with, so that 1.0 results can be reproduced. Pass one explicitly to use it, for example `--druggability assets/legacy_data/druggability_scores_annot2.0.csv`.

| File | Contents |
|---|---|
| `STRING_physical_ENSG_v12.0.csv`, `STRING_physical_geneSymbol_v12.0.csv` | STRING v12.0 physical network, combined score of at least 400 |
| `druggability_scores_annot2.0.csv` | The 1.0 druggability table, from AlphaFold DB v4 and 2024 SWISS-MODEL structures |
| `PubTator3_counts_2024-12.csv` | The 1.0 citation counts, from PubTator3 as downloaded in December 2024 |

---

## Updating the reference data

The bundled files are fixed at their release. To build newer versions, run `update_data.nf`, a separate workflow in the same repository. It downloads the current data, builds the tables with the same scripts and settings that built the bundled files, and writes them to a new folder. It never changes `assets/`.

```bash
# STRING and citation counts: minutes, on any machine
nextflow run update_data.nf -profile conda --datasets string,citations

# everything, druggability included: hours, on Linux (see below)
nextflow run update_data.nf -profile conda --datasets all --data_outdir reference_data/2027-03
```

Then point the pipeline at the new folder:

```bash
nextflow run main.nf -profile conda --scores my_scores.csv --reference_dir reference_data/2027-03
```

The folder (by default `reference_data/<date>` in the folder you run from) holds:

* the four reference files under the names the pipeline reads. A file not rebuilt in the run is copied from `--reference_dir` (default `assets/`), or kept if the folder already has it from an earlier run, so the folder is always a complete set. `ML_data.csv` is always copied: it cannot be rebuilt.
* `provenance/`: a record per dataset (source addresses, versions, download dates, sizes, SHA-256 checksums, software versions and settings) and `MANIFEST.tsv`, which gives each file's checksum and where it came from.
* `build/`: the downloads and working files. They are reused when the run is repeated, so keep the folder if you may rebuild or rescore.

**Datasets** (`--datasets`, comma-separated, or `all`):

| Dataset | Builds | Needs |
|---|---|---|
| `string` | `STRING_physical_ENSG.csv`, `STRING_physical_geneSymbol.csv` | about 30 MB of downloads, a minute |
| `citations` | `PubTator3_counts.csv` | a 760 MB download, about 4 GB of memory, a few minutes, Ensembl 116 BioMart |
| `druggability` | `druggability_scores_annot.csv` | Linux, many CPUs (fpocket and PocketMiner take about 2 hours on 32 CPUs, 4 on 16), 11 GB of downloads, about 80 GB of free disk while it runs, Ensembl 116 BioMart |

The druggability rebuild checks the operating system, the free disk and the tools it needs before it starts, and stops with a message if one falls short. At the end it packs its large working files into compressed archives (`--pack_working_files false` keeps them unpacked), which leaves about 16 GB in `build/`, most of it the two download archives. The SWISS-MODEL archive is replaced with every SWISS-MODEL release and cannot be downloaded again later, so keep it if you may need to rebuild this version.

**Options:**

| Option | Default | Meaning |
|---|---|---|
| `--data_outdir` | `reference_data/<date>` | Folder for the new set |
| `--reference_dir` | `assets/` | Where the files that are not rebuilt are copied from |
| `--data_cpus` | all CPUs | CPUs for the druggability steps |
| `--string_version`, `--string_info_version` | `12.5`, `12.0` | STRING release of the links file, and of the protein info and aliases files (STRING 12.5 reuses the 12.0 files; for a release with its own, give the same version twice) |
| `--gene2pubtator3` | download | A kept copy of `gene2pubtator3.gz`, to rebuild an earlier citation table exactly (NCBI replaces the file with each PubTator3 release) |
| `--alphafold_version` | `latest` | AlphaFold DB version to download, or a number such as `6` |
| `--alphafold_structures`, `--swissmodel_structures` | download | Score a folder of your own structures instead of downloading that source (named `AF-<UniProt ID>-F<n>-model_v<n>.pdb` and `<UniProt ID>.pdb`; checked before scoring) |
| `--min_struct_conf_af`, `--min_struct_conf_sm` | `0`, `0` | Whole-model confidence filter (mean pLDDT 0-100, mean QMEANDisCo 0-1); 0 scores every structure |
| `--min_residues` | `80` | Structures shorter than this are not scored |
| `--min_conf` | `70` | Pockets and windows count only where their mean confidence is at least this |
| `--pocketminer_window`, `--pocketminer_threshold` | `10`, `0.7` | PocketMiner window (residues on either side) and hit threshold |
| `--pack_working_files` | `true` | Pack the druggability working files at the end |
| `--min_free_disk_gb` | `80` | Free disk the druggability rebuild needs before it starts (0 skips the check) |
| `--pocketminer_python` | conda environment | A Python with PocketMiner's packages, for runs without `-profile conda` |

The defaults are the settings that built the bundled files. Changing a scoring setting gives a table built by a different method; the setting is recorded in its provenance.

**Stopping and resuming.** Run the same command again with `-resume` and the same `--data_outdir`. Finished steps are skipped, and fpocket and PocketMiner carry on from the last structure they finished. To rescore with a different `--min_residues`, `--min_conf` or PocketMiner setting, rerun into the same folder with `-resume`: the structures and the per-structure fpocket and PocketMiner outputs are reused and only the tables are rebuilt. A different whole-model filter or a new data release changes the set of structures, so it goes into a new folder: the scoring steps refuse to reuse working files made from different structures.

**Ensembl.** The citation counts and the druggability gene names are looked up in the Ensembl 116 BioMart archive. When it is unreachable, these steps wait and try again (after 10, 20 and 40 minutes) before the run stops.

**Environments.** With `-profile conda`, `update_data.nf` uses `envs/data.yaml` (R, fpocket 4.2.3, aria2, zstd) and `envs/pocketminer.yaml` (Python 3.7 with TensorFlow 2.1.0 for PocketMiner). Building the PocketMiner environment compiles mdtraj, which needs a C compiler on the system (Linux: `build-essential`; macOS: the Xcode command line tools).

Each step can also be run by hand, without Nextflow: see [`data_scripts/`](data_scripts/README.md).

> **Maintainers:** to make a new set the bundled data, run `data_scripts/promote_data.sh <set folder> --label <dataset>=<version of the bundled file>` (for example `--label string=v12.5`). It checks the new files against the bundled ones, moves the bundled files to `assets/legacy_data/` under their version names, copies the new files in and lists what is left to do (README, test snapshot, version history). The provenance records stay in the set folder, outside git.

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

> **Caching-**
> HHNet performs two levels of caching to avoid redundant computation:
>
> - **Similarity matrix** — `similarity_matrix.h5` and `beta.txt` are cached globally under `<hhnet_dir>/cache/similarity_matrices/` keyed by a hash of the indexed PPI edge list. The same PPI network always reuses the cached matrix regardless of run ID or score file.
>
> - **Intermediate files** — if `--hhnet_intermediate_cache_dir` is set, the permutation bins, permuted scores, and hierarchy files written to `intermediate/<network>_<score>/` are cached under that directory, keyed by a combined hash of the indexed PPI and score file. A cache hit skips `find_permutation_bins`, `permute_scores`, and `construct_hierarchy` entirely. This is most valuable when re-running with the same inputs but different downstream parameters (e.g. `--hhnet_min_cluster_size`).

---

### 3) Druggability (`--druggability`)

Required columns:

* `external_gene_name`
* `drug_score` (numeric, no missing values)

Scores are joined to genes by name, keeping the highest `drug_score` per gene name. Genes with no `drug_score` are not ranked; they are reported in `final_ranked_incomplete.tsv`. Other columns are ignored. The bundled table is described under [Druggability scores](#druggability-scores).

---

### 4) ML scores (`--ml_scores`)

Required columns:

* `Protein` (UniProt accession)
* `Prediction_Score_rf`

These scores are joined by UniProt accession after the ranking is calculated, so they are appended to the final output only and **cannot be used in `--ranking_features`**.

---

### 5) Citations (`--citations`)

Required columns:

* `symbol` (gene name)
* `counts`

Counts are summed per symbol and joined to genes by name. They are appended to the final output and are **not used for ranking unless included in `--ranking_features`** (as `counts`, or as `counts_norm`, which is log10 of the count). The bundled counts are described under [Citation counts](#citation-counts).

---

### 6) Optional gene mapping (`--gene_map`)

Optional TSV/CSV used to map Ensembl IDs to UniProt IDs and gene symbols, reducing or avoiding repeated annotation lookups (via Ensembl BioMart).

Recommended columns:

* `ensembl_gene_id`
* `uniprot_gn_id`

Useful additional column:

* `external_gene_name`

If `--gene_map` is not provided, the ranking stage attempts to reuse an annotation cache before falling back to online annotation lookup via Ensembl BioMart.

**Note on the BioMart lookup.** Ensembl 116 (June 2026) is the last Ensembl release that includes BioMart, and it is served from the Ensembl archive site (`jun2026.archive.ensembl.org`). The pipeline queries that BioMart service directly over HTTPS rather than through the `biomaRt` R package, because the `biomaRt` versions available through conda can no longer connect to it. Gene names, descriptions, biotypes and UniProt IDs therefore reflect Ensembl 116. Both the HHNet post-processing stage and the ranking stage use this lookup, and it needs internet access unless every gene is already in the annotation cache (or, for the ranking stage, `--gene_map` is supplied).

---

## Key parameters

### General

* `--outdir`
  Output directory.
  Default: `results`

* `--reference_dir`
  Folder holding the four reference files under the bundled names. A per-file flag (`--ppi_network`, `--druggability`, `--ml_scores`, `--citations`) replaces that one file.
  Default: `assets/` in the pipeline folder

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

* `--ranking_network` \
  Default: `neighbours`

  Options:

  * `neighbours` — use neighbour-enriched network metrics
  * `subnet` — use cluster-only subnet metrics

* `--ranking_features` \
  Default:
  `degree,betweenness,closeness,eigen_centrality,page_rank,drug_score`

  Comma-separated numeric columns used in the averaged ranking.

* `--negative_features` 
  Optional - Comma-separated features that should be treated as negatively oriented during ranking.

* `--ranking_step` 
  Optional - Passed through to the ranking script when provided.

* `--gene_map` 
  Optional - Mapping table used instead of online lookup where possible.

* `--id_annot_cache` 
  Optional - Explicit cache path for Ensembl ↔ UniProt ↔ gene symbol annotation.

* `--write_rds`
  Optional - Controls whether RDS outputs are written. In practice, the ranking module currently behaves as though this is `true` unless explicitly disabled.

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
* network artefacts generated by the post-processing script (images + graphml files)
* `id_annot_cache.tsv`

In both metrics tables, `betweenness`, `closeness`, `eigen_centrality` and `page_rank` are rounded to 10 decimal places. Digits beyond that are floating-point noise that varies with the version of the numerical libraries Conda installs. Rounding removes it, so the tables and the ranking built on them are not affected (checked across OpenBLAS 0.3.30 and 0.3.34).

---

### 3) Ranking outputs

Published to:

```text
${outdir}/ranking/
```

Main outputs:

* `final_ranked.tsv`
* `final_ranked.rds` *(optional: can disable rds writing)*
* `final_ranked_incomplete.tsv`
* `final_ranked_incomplete.rds` *(optional: can disable rds writing)*
* `final_ranked_missing_external_gene_name.csv`
* `final_ranked_missing_annotations.tsv`
* `id_annot_cache.tsv` *(optional; copied out if a cache file was written/reused)*

### Output interpretation

* `final_ranked.tsv`
  Final ranked table for complete cases included in ranking. Genes are ranked once for every combination of feature weights (rank 1 is best); `avg_rank` is a gene's mean rank across those combinations and `rank_variance` is the variance of its ranks. Under any one combination, genes with equal scores share the average of the ranks they span, so genes with identical ranking features receive the same `avg_rank`. Rows are sorted by `avg_rank`, then by gene name.

* `final_ranked_incomplete.tsv`
  Rows retained for reporting but excluded from complete ranking because one or more required ranking fields were missing (e.g. `drug_score`).

* `final_ranked_missing_external_gene_name.csv`
  Rows removed because `external_gene_name` could not be resolved.

* `final_ranked_missing_annotations.tsv`
  Ranked genes that could not be matched to one or more annotation sources (ML scores and/or citation counts). These genes are still present in `final_ranked.tsv` with `NA` in the unmatched columns. The `missing_annotation_fields` column lists which source(s) were absent (e.g. `counts`, `Prediction_Score_rf`, or `counts;Prediction_Score_rf`).

---

## Testing

This repository uses **nf-test**.

### Full pipeline smoke test

If your repository includes `tests/main.nf.test`, run:
> Use this for a quick check before committing: it runs a stub smoke test and one real end-to-end run, and compares the outputs against the stored snapshot.

```bash
nf-test test tests/main.nf.test --profile conda --verbose
```

This is useful for checking that the main workflow can execute end-to-end in the expected testing configuration.

### Module-level tests

You can also run individual module tests, for example:
> Don't worry about these tests unless your only using one component of the pipeline.

```bash
nf-test test tests/modules/local/hhnet.nf.test --profile conda --verbose
nf-test test tests/modules/local/hhnet_post.nf.test --profile conda --verbose
nf-test test tests/modules/local/ranking.nf.test --profile conda --verbose
```

### Running all tests

If your repository is organised with the standard `tests/` layout, you can run all tests with:
> Use this after changing a single stage: it runs the full pipeline test above plus the three module tests, so it takes noticeably longer.

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

This is expected and can be useful when re-running the pipeline. You may experience issues connecting to the Ensembl servers here in which case retry or try again later.

### Internet access

Two stages look up gene annotation online (see the note under "Optional gene mapping"):

* **HHNet post-processing** queries Ensembl BioMart on every run that is not resumed from Nextflow's cache, so it needs internet access.
* **Ranking** needs internet access unless `--gene_map` is provided or its annotation cache already covers every gene.

### Caching

**Similarity matrix**

Similarity matrices generated by HHNet are cached under `<hhnet_dir>/cache/similarity_matrices/` keyed by a hash of the indexed PPI edge list. Reusing the same processed PPI network skips the eigendecomposition step on future runs, provided the cache location is writable.

**Intermediate files**

If `--hhnet_intermediate_cache_dir` is provided, the permutation bins, permuted scores, and hierarchy files are cached at that path, keyed by a combined SHA-256 hash of the indexed PPI and score file. On a cache hit the three most expensive parallel steps are skipped entirely and their outputs are linked into the current run directory before `process_hierarchies` executes as normal.

A stable path outside the Nextflow work directory is recommended so the cache persists across `-resume` runs and `nextflow clean`:

```groovy
// nextflow.config
params.hhnet_intermediate_cache_dir = "${projectDir}/cache/hhnet_intermediates"
```

Add `cache/` to `.gitignore` to keep it out of version control.

---

## Version history

Key changes between versions, newest first.

### 2.0 (November 2026)

The first update to the original version. **Rankings from 2.0 are not byte identical to rankings from 1.0.**

* **Gene annotation queries Ensembl BioMart directly.**
  Ensembl 116 (June 2026) is the last Ensembl release that includes BioMart, and the `biomaRt` R package versions available through Conda can no longer connect to it. Both R stages now send their queries straight to the Ensembl 116 archive over HTTPS, and the R environment installs `httr` instead of `biomaRt`. Gene names, descriptions, biotypes and UniProt IDs are therefore fixed at Ensembl 116. See the note under "Optional gene mapping" above.

* **Network metrics are rounded to 10 decimal places.**
  In 1.0, `betweenness`, `closeness`, `eigen_centrality` and `page_rank` were written at full floating-point precision. Their last digits changed with the version of the linear algebra library (OpenBLAS) that Conda happened to install, and that was enough to change the ranking from one installation to the next. With rounding, the OpenBLAS versions installed in March and in November 2026 (0.3.30 and 0.3.34) give identical metrics tables and an identical ranking from the same inputs.

* **Tied genes share a rank.**
  Under each combination of feature weights, genes with equal scores now receive the average of the ranks they span. In 1.0 they were ranked in the order they appeared in the metrics table, which follows the row order of the PPI file. `final_ranked.tsv` is now sorted by `avg_rank`, then by gene name.

* **Updated for Nextflow 26.10.**
  The workflow now uses the current `channel.fromPath` spelling in place of the deprecated `Channel.fromPath`, and the two R scripts are called by name from `bin/` in the same way as `run_hhnet.sh`, so they must be executable (see Installation, step 6). `nextflow lint` reports no errors or warnings. The full test suite passes on Nextflow 26.10.5 with nf-test 0.9.5, with outputs identical to those from Nextflow 26.04.6.

* **Updated reference data.** The new files are described under [Bundled reference data](#bundled-reference-data); the files that shipped with 1.0 are kept in `assets/legacy_data/` (see [Legacy data](#legacy-data)).
  * *PPI network:* STRING 12.0 to 12.5, at the same cut-off (combined score of at least 400). Each interaction is now listed once.
  * *Druggability:* rebuilt from new structures (AlphaFold DB v4 to v6, SWISS-MODEL Repository 2024 to 2026), with fpocket and PocketMiner run again on every structure. Pockets and windows now count only in confidently modelled regions, which replaces the filter on whole-model confidence, and the 80-residue minimum applies to both tools. Gene names come from Ensembl 116, the release the pipeline names its genes with. The file is renamed `druggability_scores_annot.csv` (was `druggability_scores_annot2.0.csv`) and columns are renamed: `uniprot_gn_id` to `uniprot_id`, `Fpocket_method` to `FP_source`, `struct_score` to `FP_struct_score`, `PocketMiner_method` to `CP_source` and `source` to `drug_score_source`; `drug_score_tool` is new. The pipeline reads only `external_gene_name` and `drug_score`.
  * *Citation counts:* regenerated from PubTator3 as downloaded on 13 November 2026. Every article is counted, for human genes only, named from Ensembl 116.
  * *ML scores:* unchanged.

* **Bundled reference data used by default.** `--ppi_network`, `--druggability`, `--ml_scores` and `--citations` default to the files in `assets/`, so `--scores` is the only required input. `--reference_dir` points a run at another folder of reference files, and each run prints the reference files it uses.

* **Updating the reference data.** `update_data.nf` builds new versions of the PPI network, druggability and citation files into a new folder, used with `--reference_dir` (see [Updating the reference data](#updating-the-reference-data)). The scripts it runs are in `data_scripts/`, with instructions for running them by hand in its README.

* **Input formats corrected.** The input-format section now lists the columns the pipeline actually requires, and says that ML scores cannot be used as a ranking feature.

* **The test snapshot was re-recorded** to match the new outputs.

### 1.0 (March 2026)

The original version: HHNet clustering, post-processing into network metrics, and integrated evidence ranking, run through Nextflow with Conda. Gene annotation used the `biomaRt` R package against the live Ensembl BioMart service.

---

## Acknowledgement

NetIntegRank is based on work developed during my PhD, with a focus on dysregulated subnetwork discovery and integrated evidence ranking for small molecule target prioritisation in cancer.