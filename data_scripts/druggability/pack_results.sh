#!/bin/bash


# example commands: ./pack_results.sh
#                   ./pack_results.sh path/to/druggability_build -t 16
#
# arguments: the folder holding alphafold/ and swissmodel/ (default: this script's folder, the default output
#            location of the other scripts); -t threads for zstd (default: all CPUs)
#
# Last step of the druggability rebuild, after package_table.R. Packs the large working folders of each source,
# which only a rerun needs, into one archive each: fpocket's structural/, PocketMiner's txt_results/ and
# pocket_results/, and the prepared structures in structures_filtered/pdb/. Each archive is checked against its
# folder before the folder is removed (src/packed_dirs.sh). Then structures_raw/ is removed where the download
# archive and provenance record it was unpacked from are next to it: download_alphafold_structures.sh and
# download_swissmodel_structures.sh unpack it again when rerun. Structures you supplied yourself (anywhere
# else) are never touched, and the download archives are never touched.
#
# The tables (fpocket_druggability_full.csv, pocketminer_results.csv, small_structures.csv, the run_info.txt
# records) stay as they are. Fpocket_master.sh and pocketminer_master.sh unpack their folders when rerun, so
# resuming and rescoring still work. Uses zstd if it is installed, otherwise gzip (larger archives).



### pack the druggability working folders ###
if [ -z "${BASH_VERSION:-}" ]; then echo "Run with bash, not sh"; exit 1; fi
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$script_dir/src/packed_dirs.sh"

build_dir="$script_dir"
if [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; then build_dir="$1"; shift; fi
threads="$(nproc)"
while getopts "t:" opt; do
  case "$opt" in
    t) threads="$OPTARG" ;;
    *) echo "Usage: $0 [build_dir] [-t threads]"; exit 1 ;;
  esac
done
[ -d "$build_dir" ] || { echo "$build_dir does not exist"; exit 1; }
build_dir="$(cd "$build_dir" && pwd)"
if command -v zstd > /dev/null; then format="zstd -3 ($(zstd -V | head -n 1))"; else format="gzip (zstd not installed)"; fi

record="$build_dir/packed.txt"
{ echo "packed by pack_results.sh on $(date -u '+%Y-%m-%d %H:%M UTC'), format: $format"; } >> "$record"
total_before=0; total_after=0
for source in alphafold swissmodel; do
  src_dir="$build_dir/$source"
  [ -d "$src_dir" ] || continue
  for item in fpocket_results/structural pocketminer_results/txt_results pocketminer_results/pocket_results \
              structures_filtered/pdb; do
    parent="$src_dir/$(dirname "$item")"; name="$(basename "$item")"
    [ -d "$parent/$name" ] || continue
    before=$(du -sb "$parent/$name" | cut -f1)
    archive="$(pack_dir "$parent" "$name" "$threads")"
    after=$(stat -c %s "$parent/$archive")
    total_before=$(( total_before + before )); total_after=$(( total_after + after ))
    printf '  %-45s %8d MB -> %6d MB  %s\n' "$source/$item/" $(( before / 1048576 )) $(( after / 1048576 )) "$archive" \
      | tee -a "$record"
  done

  # unpacked structures: removed only where the download archive and its provenance record are beside them
  raw="$src_dir/structures_raw"
  [ -d "$raw" ] || continue
  if [ "$source" = alphafold ]; then
    archive_glob=("$src_dir"/*_v*.tar); prov_glob=("$src_dir"/*_v*.provenance.txt)
  else
    archive_glob=("$src_dir"/*_coords.tar.gz); prov_glob=("$src_dir"/swissmodel_*.provenance.txt)
  fi
  if [ -f "${archive_glob[0]}" ] && [ -f "${prov_glob[0]}" ]; then
    before=$(du -sb "$raw" | cut -f1)
    rm -rf "${raw:?}"
    total_before=$(( total_before + before ))
    printf '  %-45s %8d MB removed (kept in %s)\n' "$source/structures_raw/" $(( before / 1048576 )) \
      "$(basename "${archive_glob[0]}")" | tee -a "$record"
  else
    echo "  $source/structures_raw/ kept: no download archive and provenance record beside it" | tee -a "$record"
  fi
done
echo "  total: $(( total_before / 1048576 )) MB -> $(( total_after / 1048576 )) MB" | tee -a "$record"
