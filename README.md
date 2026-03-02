# NetIntegRank

Dysregulated subnetworks → integrated evidence ranking.

## Overview
NetIntegRank runs two stages:
1) HHnet subnetwork inference from user-provided differential expression results and a PPI network
2) Integrated ranking using network metrics + druggability + ML scores + citation evidence

## Inputs
- DE results (user-supplied)
- PPI network (default provided; optional user override)
- Druggability table
- ML score table
- Citation score table
- Gene ID mapping table (precomputed; no biomaRt dependency)
