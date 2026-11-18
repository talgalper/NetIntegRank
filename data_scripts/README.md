# Data scripts

Scripts that download and build the reference data in `assets/`. The pipeline itself does not use this folder. Run these only to rebuild or update the bundled files.

| Script | Builds | Runs on |
|---|---|---|
| `update_STRING_data.R` | `STRING_physical_ENSG.csv`, `STRING_physical_geneSymbol.csv` | any machine with R (about 30 MB download, under a minute) |
| `update_citation_data.R` | `PubTator3_counts.csv` | any machine with R and about 4 GB of free memory (760 MB download, a few minutes) |
| `druggability/` | `druggability_scores_annot.csv` | Linux, many CPUs and about 60 GB of free disk (structures and pocket-detection runs) |

Each script writes into its own working folder. Check the result, then copy it into `assets/` under the name the pipeline uses. The methods and data versions behind the bundled files are described under "Bundled reference data" in the main [README](../README.md).

---

## STRING network

```bash
Rscript update_STRING_data.R
```

* Downloads the human (9606) STRING files into `STRING_data/` in the current folder: `protein.physical.links.full.v12.5`, plus `protein.info.v12.0` and `protein.aliases.v12.0` (STRING 12.5 reuses the 12.0 versions of these).
* Keeps one row per interaction (STRING lists each one in both directions) and interactions with a `combined_score` of at least 400, stored divided by 1000.
* Maps STRING protein IDs to Ensembl gene IDs (aliases with source `Ensembl_gene`) and to gene symbols (protein info).
* Writes `STRING_data/STRING_physical_ENSG.csv` and `STRING_data/STRING_physical_geneSymbol.csv`.
* R package: `data.table`.

## Citation counts

```bash
Rscript update_citation_data.R
```

* Downloads PubTator3's gene annotations (`gene2pubtator3.gz`) into `citation_data/` in the current folder. NCBI replaces this file with every PubTator3 release, so a later download gives somewhat different counts.
* Counts the distinct PubMed articles that mention each NCBI gene ID, whatever the topic, and leaves out mentions linked to more than one gene.
* Names the gene IDs through Ensembl 116 (the release the pipeline names its genes with), which keeps human genes only, and sums the counts per gene name.
* Writes `citation_data/PubTator3_counts_noFilt_<date>.csv` (columns `symbol`, `counts`); copy it to `assets/PubTator3_counts.csv`.
* R packages: `data.table`, `httr`. Needs internet access for Ensembl 116.

---

## Druggability scores

Everything below runs from inside `druggability/`. With default paths each step reads and writes per-source folders next to the scripts (`alphafold/`, `swissmodel/`), which are ignored by git. The three drivers (`filter_structures.sh`, `Fpocket_master.sh`, `pocketminer_master.sh`) also take `-i` and `-o` folders and `-t` threads; run one with no arguments to see its usage. The settings each step uses (confidence cut-offs, minimum length, PocketMiner window) are at the top of the drivers and are written into each results folder's `run_info.txt`, together with the tool versions and the checksums of the inputs.

The shell drivers need Linux (GNU `find`, `xargs`, `sha256sum`, `nproc`) and R with `bio3d`, `readr` and `httr`.

### 1. Get the structures

**AlphaFold DB** (human reference proteome; about 5 GB to download, 10 GB unpacked):

```bash
./download_alphafold_structures.sh            # newest AlphaFold DB version; AF_VERSION=6 for a specific one
```

It needs GNU `tar`, `curl`, `python3` and `gzip`, and `aria2c` for a fast download. It checks the archive, unpacks the PDB files into `alphafold/structures_raw/` and writes a provenance record next to the archive.

**SWISS-MODEL Repository** (human models, https://swissmodel.expasy.org/repository; about 6 GB). These steps are not scripted yet:

```bash
mkdir -p swissmodel && cd swissmodel
curl -O https://swissmodel.expasy.org/repository/download/core_species/9606_coords.tar.gz
curl -O https://swissmodel.expasy.org/repository/download/core_species/9606_meta.tar.gz
tar -xzf 9606_coords.tar.gz && tar -xzf 9606_meta.tar.gz     # gives SWISS-MODEL_Repository/
mkdir -p structures_raw && Rscript ../swissmodel_sort.R      # R packages: tidyverse, bio3d, progress
cd ..
```

`swissmodel_sort.R` keeps, for each protein, the model with a sequence identity of at least 60% and the best QMEANDisCo global score, and copies it to `swissmodel/structures_raw/<UniProt ID>.pdb`. SWISS-MODEL replaces these archives with each release, so record their size and SHA-256.

### 2. Filter by overall confidence

```bash
./filter_structures.sh alphafold
./filter_structures.sh swissmodel
```

Keeps structures with a mean pLDDT of at least 50 (AlphaFold) or a mean QMEANDisCo of at least 0.5 (SWISS-MODEL), and trims SWISS-MODEL models to their first chain. Writes `<source>/structures_filtered/` (`pdb/`, `struct_score.csv`, `low_conf_struct.csv`).

### 3. fpocket

```bash
./Fpocket_master.sh alphafold
./Fpocket_master.sh swissmodel
```

On the first run it builds fpocket 4.2.3 into `fpocket-4.2.3/` (needs `git`, `make`, a C compiler and `libnetcdf-dev`). Writes `<source>/fpocket_results/`: `fpocket_druggability_full.csv` (one row per structure, with the confident-pocket columns), `small_structures.csv` (structures under 80 residues) and `fpocket_pockets.csv` (every pocket with its confidence and lining residues). A stopped run carries on where it left off.

### 4. PocketMiner

```bash
./pocketminer_master.sh alphafold
./pocketminer_master.sh swissmodel
```

Needs a `pocketminer` conda environment with Python 3.7, TensorFlow 2.1.0, numpy 1.18.1 and mdtraj 1.9.9; set `POCKETMINER_PYTHON` if its Python is not at `~/miniconda3/envs/pocketminer/bin/python`. (An environment file will be added to this folder.) The authors' code and trained model are in `pocketminer_src/upstream/` (Meller et al. 2023, https://github.com/Mickdub/gvp, commit `187062d`, MIT licence); the driver checks that every model weight has loaded before it scores anything. Writes `<source>/pocketminer_results/`: `pocketminer_results.csv`, `small_structures.csv` and per-structure predictions. About 20 minutes for the AlphaFold set on 32 CPUs.

### 5. Combine, annotate and package

```bash
mkdir -p combined annotated package

Rscript combine_scores.R \
  --af-fpocket alphafold/fpocket_results/fpocket_druggability_full.csv,alphafold/fpocket_results/small_structures.csv \
  --sm-fpocket swissmodel/fpocket_results/fpocket_druggability_full.csv,swissmodel/fpocket_results/small_structures.csv \
  --af-pocketminer alphafold/pocketminer_results/pocketminer_results.csv \
  --sm-pocketminer swissmodel/pocketminer_results/pocketminer_results.csv \
  --min-residues 80 --fpocket-score druggability_conf --pocketminer-score max_hit_conf \
  --out combined/drug_scores.csv

Rscript annotate_scores.R --in combined/drug_scores.csv --out annotated/drug_scores_annotated.csv

Rscript package_table.R --in annotated/drug_scores_annotated.csv --out package/druggability_scores_annot.csv
```

* `combine_scores.R` takes each protein's best confident fpocket score and best confident PocketMiner score over its structures, ignoring structures under 80 residues; `drug_score` is the higher of the two, and a protein with no usable score gets 0.
* `annotate_scores.R` adds Ensembl 116 gene names, descriptions and biotypes (needs internet access).
* `package_table.R` gives the columns their final names and order. Copy `package/druggability_scores_annot.csv` to `assets/`.
