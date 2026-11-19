# Data scripts

Scripts that download and build the reference data in `assets/`. The pipeline itself does not use this folder. `update_data.nf` runs these scripts for you, with the settings that built the bundled files (see "Updating the reference data" in the main [README](../README.md)); this page is for running them by hand, one step at a time.

| Script | Builds | Runs on |
|---|---|---|
| `update_STRING_data.R` | `STRING_physical_ENSG.csv`, `STRING_physical_geneSymbol.csv` | any machine with R (about 30 MB download, under a minute) |
| `update_citation_data.R` | `PubTator3_counts.csv` | any machine with R and about 4 GB of free memory (760 MB download, a few minutes) |
| `druggability/` | `druggability_scores_annot.csv` | Linux, many CPUs and about 80 GB of free disk while it runs (structures and pocket-detection runs) |

Each script writes into its own working folder, with a provenance record (sources, versions, dates, sizes and SHA-256 checksums) next to its output. To use a result, put the files in a folder under the names the pipeline reads and run the pipeline with `--reference_dir <folder>` (copy `ML_data.csv` from `assets/` into it as well), or give a single file with its own option, such as `--ppi_network`. The methods and data versions behind the bundled files are described under "Bundled reference data" in the main README.

`promote_data.sh` is for the maintainer: it makes a set built by `update_data.nf` the bundled data (see the comments at its top).

---

## STRING network

```bash
Rscript update_STRING_data.R                                  # STRING 12.5, as bundled
Rscript update_STRING_data.R --links-version 13.0 --info-version 13.0 --out-dir STRING_13
```

* Downloads the human (9606) STRING files into `STRING_data/` (or `--out-dir`) in the current folder: the physical links file of `--links-version` (default 12.5), and the protein info and aliases files of `--info-version` (default 12.0: STRING 12.5 reuses the 12.0 versions of these; for a release that ships its own, give the same version for both).
* Keeps one row per interaction (STRING lists each one in both directions) and interactions with a `combined_score` of at least 400, stored divided by 1000.
* Maps STRING protein IDs to Ensembl gene IDs (aliases with source `Ensembl_gene`) and to gene symbols (protein info).
* Writes `STRING_physical_ENSG.csv`, `STRING_physical_geneSymbol.csv` and `STRING.provenance.txt`.
* R package: `data.table`.

## Citation counts

```bash
Rscript update_citation_data.R
Rscript update_citation_data.R --gene2pubtator3 kept/gene2pubtator3.gz --out PubTator3_counts.csv
```

* Downloads PubTator3's gene annotations (`gene2pubtator3.gz`) into `citation_data/` (or `--out-dir`) in the current folder. NCBI replaces this file with every PubTator3 release, so a later download gives somewhat different counts; keep the download, or give a kept copy with `--gene2pubtator3`, to rebuild a table exactly.
* Counts the distinct PubMed articles that mention each NCBI gene ID, whatever the topic, and leaves out mentions linked to more than one gene.
* Names the gene IDs through Ensembl 116 (the release the pipeline names its genes with), which keeps human genes only, and sums the counts per gene name.
* Writes `PubTator3_counts_noFilt_<date>.csv` (or `--out`; columns `symbol`, `counts`) and `citations.provenance.txt`. The pipeline reads the file as `PubTator3_counts.csv`.
* R packages: `data.table`, `httr`. Needs internet access for Ensembl 116.

---

## Druggability scores

Everything below runs from inside `druggability/`. With default paths each step reads and writes per-source folders next to the scripts (`alphafold/`, `swissmodel/`), which are ignored by git. The three drivers (`filter_structures.sh`, `Fpocket_master.sh`, `pocketminer_master.sh`) also take `-i` and `-o` folders and `-t` threads; run one with no arguments to see its usage. The settings each step uses (confidence cut-offs, minimum length, PocketMiner window) are at the top of the drivers, can be changed for one run with the options below, and are written into each results folder's `run_info.txt`, together with the tool versions and the checksums of the scripts and inputs.

The shell drivers need Linux (GNU `find`, `xargs`, `tar`, `sha256sum`, `nproc`) and R with `bio3d`, `readr`, `httr` and `progress`. The conda environment `../../envs/data.yaml` has all of them, plus fpocket 4.2.3, aria2 and zstd:

```bash
conda env create -n netintegrank-data -f ../../envs/data.yaml && conda activate netintegrank-data
```

### 1. Get the structures

**AlphaFold DB** (human reference proteome; about 5 GB to download, 10 GB unpacked):

```bash
./download_alphafold_structures.sh            # newest AlphaFold DB version; AF_VERSION=6 for a specific one
```

It needs GNU `tar`, `curl`, `python3` and `gzip`, and `aria2c` for a fast download. It checks the archive, unpacks the PDB files into `alphafold/structures_raw/` and writes a provenance record next to the archive.

**SWISS-MODEL Repository** (human models, https://swissmodel.expasy.org/repository; about 6 GB to download):

```bash
./download_swissmodel_structures.sh
```

It downloads `9606_coords.tar.gz` and `9606_meta.tar.gz` into `swissmodel/`, checks them and records their size and SHA-256. It then picks, for each protein, the model with a sequence identity of at least 60% and the best QMEANDisCo global score, unpacks only those models (about 10,000 of the 91,000 in the archive) and copies each to `swissmodel/structures_raw/<UniProt ID>.pdb` with `swissmodel_sort.R`. SWISS-MODEL replaces these archives with each release and an old one cannot be downloaded again, so keep them if you may need to rebuild this version. An archive already in `swissmodel/` is used as it is.

**Your own structures** can be scored instead: put them in a folder, named `AF-<UniProt ID>-F<fragment>-model_v<version>.pdb` (AlphaFold) or `<UniProt ID>.pdb` (SWISS-MODEL), with per-residue confidence in the B-factor column (pLDDT 0-100, or QMEANDisCo 0-1), and give the folder to step 2 with `-i`. Check a folder first with:

```bash
./check_structures.sh alphafold path/to/structures/
./check_structures.sh swissmodel path/to/structures/
```

### 2. Score overall confidence and prepare the structures

```bash
./filter_structures.sh alphafold
./filter_structures.sh swissmodel
```

Gives each structure its overall confidence (mean pLDDT, or mean QMEANDisCo for SWISS-MODEL) and trims SWISS-MODEL models to their first chain. Every readable structure is kept, because the scoring steps leave out pockets and windows in low-confidence regions; `-c` sets a minimum overall confidence if you want one (for example `-c 50` for AlphaFold, `-c 0.5` for SWISS-MODEL). Writes `<source>/structures_filtered/` (`pdb/`, `struct_score.csv`, and `low_conf_struct.csv` for any structures below `-c`).

### 3. fpocket

```bash
./Fpocket_master.sh alphafold
./Fpocket_master.sh swissmodel
```

Uses the fpocket given by `FPOCKET_BIN` if it is fpocket 4.2.3 installed by conda (`FPOCKET_BIN=$CONDA_PREFIX/bin/fpocket` in the data environment); otherwise it builds fpocket 4.2.3 into `fpocket-4.2.3/` on the first run (needs `git`, `make`, a C compiler and `libnetcdf-dev`). The two give the same results; only fpocket's pocket volume estimate, which it makes by random sampling and which no table uses, differs between any two runs. Writes `<source>/fpocket_results/`: `fpocket_druggability_full.csv` (one row per structure, with the confident-pocket columns), `small_structures.csv` (structures under 80 residues) and `fpocket_pockets.csv` (every pocket with its confidence and lining residues). `-m` sets the minimum length (default 80) and `-c` the confidence cut-off (default 70). A stopped run carries on where it left off; a results folder is only reused for exactly the same structure files, so use a new one for new structures.

### 4. PocketMiner

```bash
conda env create -n pocketminer -f ../../envs/pocketminer.yaml   # once; builds mdtraj, so it needs a C compiler
./pocketminer_master.sh alphafold
./pocketminer_master.sh swissmodel
```

The `pocketminer` environment has Python 3.7, TensorFlow 2.1.0, numpy 1.18.1 and mdtraj 1.9.9. The driver looks for its Python at `~/miniconda3/envs/pocketminer/bin/python`; set `POCKETMINER_PYTHON` if it is elsewhere. The authors' code and trained model are in `pocketminer_src/upstream/` (Meller et al. 2023, https://github.com/Mickdub/gvp, commit `187062d`, MIT licence); the driver checks that every model weight has loaded before it scores anything. Writes `<source>/pocketminer_results/`: `pocketminer_results.csv`, `small_structures.csv` and per-structure predictions. `-m`, `-w`, `-s` and `-c` set the minimum length (80), the window on either side (10), the hit threshold (0.7) and the confidence cut-off (70). About 20 minutes for the AlphaFold set on 32 CPUs. Like fpocket, a stopped run carries on, and a results folder is only reused for the same structures.

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

* `combine_scores.R` takes each protein's best confident fpocket score and best confident PocketMiner score over its structures, ignoring structures under 80 residues; `drug_score` is the higher of the two, and a protein with no usable score gets 0. These three settings are its defaults, written out above so they are visible; `--fpocket-score druggability`, `--pocketminer-score max_hit` and `--min-residues 0` use all pockets, all windows and every structure instead.
* `annotate_scores.R` adds Ensembl 116 gene names, descriptions and biotypes (needs internet access).
* `package_table.R` gives the columns their final names and order. The pipeline reads the result as `druggability_scores_annot.csv`.

### 6. Pack the working files

```bash
./pack_results.sh
```

Packs the large working folders, which only a rerun needs, into one archive each (zstd, or gzip if zstd is not installed): fpocket's `structural/`, PocketMiner's `txt_results/` and `pocket_results/`, and `structures_filtered/pdb/`. Each archive is checked against its folder before the folder is removed. It also removes the unpacked `structures_raw/` folders, which the download archives beside them still hold. This takes a full run from about 80 GB to about 16 GB, most of it the two download archives. The tables stay as they are. To run the steps again, run the download scripts first (they unpack the structures from the archives again); `Fpocket_master.sh` and `pocketminer_master.sh` unpack their own folders.
