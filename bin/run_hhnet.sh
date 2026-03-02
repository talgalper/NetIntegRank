#!/usr/bin/env bash
set -euo pipefail

while [[ $# -gt 0 ]]; do
  case "$1" in
    --de) DE="$2"; shift 2 ;;
    --ppi) PPI="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# IMPORTANT: Make sure to double check and update paths and filenames
# Set the paths for your data, intermediate, and results directories
data=$PWD/MOC/data
intermediate=$PWD/MOC/intermediate
results=$PWD/MOC/results

# Create directories
mkdir -p $data
mkdir -p $intermediate
mkdir -p $results

# Set the name of your network and score file
network=STRING
score=MOC_logFC

# Filenames for your input files
index_gene_file=gene_index.tsv
edge_list_file=edge_list_index.tsv
score_file=logFC_scores_abs.tsv

num_permutations=100

# This script uses gnu parallel to parallelize some of the scripts. If gnu
# parallel is not already installed on your system, then you can install it from
# https://www.gnu.org/software/parallel/. You can change the num_cores variable
# to specify the number of cores for your system.

##### MAKE SURE THIS IS CORRECT #####
num_cores=36

mkdir -p $intermediate/"$network"
mkdir -p $intermediate/"$network"_"$score"

# Compile Fortran module (if needed)
#cd ../src
#f2py -c fortran_module.f95 -m fortran_module > /dev/null
#cd ..

# Construct similarity matrix
echo "Construct similarity matrix..."
python src/construct_similarity_matrix.py \
    -i   $data/"$edge_list_file" \
    -o   $intermediate/"$network"/similarity_matrix.h5 \
    -bof $intermediate/"$network"/beta.txt

# Permute scores
echo "Permuting scores..."
cp $data/"$score_file" $intermediate/"$network"_"$score"/scores_0.tsv

python src/find_permutation_bins.py \
    -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
    -igf $data/"$index_gene_file" \
    -elf $data/"$edge_list_file" \
    -ms  1000 \
    -o   $intermediate/"$network"_"$score"/score_bins.tsv

parallel -u -j $num_cores --bar \
    python src/permute_scores.py \
        -i  $intermediate/"$network"_"$score"/scores_0.tsv \
        -bf $intermediate/"$network"_"$score"/score_bins.tsv \
        -s  {} \
        -o  $intermediate/"$network"_"$score"/scores_{}.tsv \
    ::: `seq $num_permutations`

# Construct hierarchies
echo "Constructing hierarchies..."
parallel -u -j $num_cores --bar \
    python src/construct_hierarchy.py \
        -smf  $intermediate/"$network"/similarity_matrix.h5 \
        -igf  $data/"$index_gene_file" \
        -gsf  $intermediate/"$network"_"$score"/scores_{}.tsv \
        -helf $intermediate/"$network"_"$score"/hierarchy_edge_list_{}.tsv \
        -higf $intermediate/"$network"_"$score"/hierarchy_index_gene_{}.tsv \
    ::: `seq 0 $num_permutations`

# Process hierarchies
echo "Processing hierarchies..."
python src/process_hierarchies.py \
    -oelf $intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv \
    -oigf $intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv \
    -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv "; done) \
    -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv "; done) \
    -lsb  10 \
    -cf   $results/clusters_"$network"_"$score".tsv \
    -pl   $network $score \
    -pf   $results/sizes_"$network"_"$score".pdf \
    -nc   $num_cores

# Placeholder: replace with your real HHnet invocation.
# Example outputs expected by Nextflow:
#   ${OUTDIR}/metrics.tsv
#   ${OUTDIR}/subnetwork.tsv (optional)

echo -e "node\tdegree\tbetweenness\tcloseness\teigenvector" > "${OUTDIR}/metrics.tsv"
