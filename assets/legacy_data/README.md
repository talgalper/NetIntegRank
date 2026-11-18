# Legacy data

The reference files that NetIntegRank 1.0 shipped with in `assets/`, kept so that 1.0 results can be reproduced. Version 2.0 uses the updated files in `assets/`, described under "Bundled reference data" in the main [README](../../README.md).

| File | Was in 1.0 as | Contents |
|---|---|---|
| `STRING_physical_ENSG_v12.0.csv` | `assets/STRING_physical_ENSG.csv` | STRING v12.0 physical network, combined score of at least 400, Ensembl gene IDs |
| `STRING_physical_geneSymbol_v12.0.csv` | `assets/STRING_physical_geneSymbol.csv` | The same network with gene symbols |
| `druggability_scores_annot2.0.csv` | `assets/druggability_scores_annot2.0.csv` | Druggability table from AlphaFold DB v4 and 2024 SWISS-MODEL structures |
| `PubTator3_counts_2024-12.csv` | `assets/PubTator3_counts.csv` | PubTator3 citation counts, downloaded in December 2024 |

To use one, pass it explicitly, for example:

```bash
nextflow run main.nf -profile conda -params-file params/example_default.yaml \
  --druggability assets/legacy_data/druggability_scores_annot2.0.csv
```
