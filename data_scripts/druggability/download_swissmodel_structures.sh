#!/usr/bin/env bash
#
# download_swissmodel_structures.sh
#
# Downloads the SWISS-MODEL Repository's models for one species (human by
# default), checks them, and unpacks only the model each protein is scored on.
#
# Usage
#   ./download_swissmodel_structures.sh
#   OUT_DIR=/some/folder ./download_swissmodel_structures.sh
#
# What it does
#   1. Downloads <taxid>_coords.tar.gz (all models, about 6 GB for human) and
#      <taxid>_meta.tar.gz (the INDEX of models) into OUT_DIR. Safe to re-run:
#      an interrupted download carries on, a finished one is not fetched again.
#      An archive already in OUT_DIR is always used as it is, even if the
#      server now has a newer release: SWISS-MODEL replaces these files with
#      each release and an old one cannot be downloaded again, so keep them.
#   2. Checks both archives (gzip's built-in error check; SWISS-MODEL
#      publishes no checksum) and records their size and SHA-256.
#   3. Chooses one model per protein from the INDEX, by swissmodel_sort.R's
#      rule (sequence identity >= 60, best QMEANDisCo global score), unpacks
#      only those models, and runs swissmodel_sort.R on them, which copies
#      each to OUT_DIR/structures_raw/<UniProt ID>.pdb. The full repository
#      (about 91,000 models, 26 GB unpacked for human) is never unpacked.
#   4. Writes a provenance record next to the archives.
#
# It never deletes or changes the archives. If OUT_DIR/structures_raw holds
# structures made from other archives, it stops and says so: use a separate
# OUT_DIR for each release.
#
# Needs: Linux (GNU tar, find, mv), curl, gzip, sha256sum, Rscript (R packages bio3d
# and progress). aria2c is optional and makes the download faster.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------- settings ---
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/swissmodel}"   # where the archives and structures_raw go
SM_TAXID="${SM_TAXID:-9606}"                   # NCBI taxonomy ID; 9606 is human
CONNECTIONS="${CONNECTIONS:-8}"                # parallel connections for aria2c

BASE_URL="https://swissmodel.expasy.org/repository/download/core_species"

# ----------------------------------------------------------------- helpers ---
log() { printf '[%s] %s\n' "$(date -u '+%H:%M:%S')" "$*"; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

for tool in curl gzip tar Rscript; do
  command -v "$tool" >/dev/null || die "'$tool' is needed but was not found."
done
tar --version 2>/dev/null | grep -q 'GNU tar' \
  || die "GNU tar is needed (this is not GNU tar). On a Mac: brew install gnu-tar, then put it first on the PATH."
if command -v sha256sum >/dev/null; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
[[ "$SM_TAXID" =~ ^[0-9]+$ ]] || die "SM_TAXID must be a number, got '$SM_TAXID'."

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
COORDS="$OUT_DIR/${SM_TAXID}_coords.tar.gz"
META="$OUT_DIR/${SM_TAXID}_meta.tar.gz"
STRUCT_DIR="$OUT_DIR/structures_raw"
PROVENANCE="$OUT_DIR/swissmodel_${SM_TAXID}.provenance.txt"
TMP="$OUT_DIR/.unpack_tmp"   # this script's own scratch folder, removed at the end

# ------------------------------------------------------------ 1. download ---
# A finished download has its final name; a part-download is <name>.part.
server_info() {   # prints "<size> <last-modified>" for a URL, or nothing if the server cannot be reached
  curl -sSIL --retry 3 --max-time 120 "$1" 2>/dev/null | awk '
    toupper($1) ~ /^HTTP/ {size = ""; date = ""}
    tolower($1) == "content-length:" {gsub("\r", "", $2); size = $2}
    tolower($1) == "last-modified:" {sub(/^[^:]*: */, ""); gsub("\r", ""); date = $0}
    END {if (size != "") print size, date}'
}

fetch() {   # fetch <url> <file>
  local url="$1" file="$2" info size
  if [ -f "$file" ]; then
    info="$(server_info "$url" || true)"
    size="${info%% *}"
    if [ -n "$size" ] && [ "$size" != "$(wc -c < "$file" | tr -d ' ')" ]; then
      log "NOTE: the server's $(basename "$file") ($size bytes) differs from the copy here; using the copy here"
    else
      log "$(basename "$file") is already here; not downloading again"
    fi
    return 0
  fi
  info="$(server_info "$url")" || true
  [ -n "$info" ] || die "Could not reach $url"
  size="${info%% *}"
  log "Downloading $url ($size bytes)"
  if command -v aria2c >/dev/null; then
    aria2c -c -x "$CONNECTIONS" -s "$CONNECTIONS" --file-allocation=none \
           --summary-interval=0 --console-log-level=warn \
           -d "$(dirname "$file")" -o "$(basename "$file").part" "$url" \
      || die "The download did not finish. Run the script again to carry on from where it stopped."
  else
    curl -fL --retry 5 -C - -o "$file.part" "$url" \
      || die "The download did not finish. Run the script again to carry on from where it stopped."
  fi
  [ "$(wc -c < "$file.part" | tr -d ' ')" = "$size" ] && [ ! -e "$file.part.aria2" ] \
    || die "The download of $(basename "$file") looks unfinished. Run the script again."
  mv "$file.part" "$file"
  echo "$info" > "$file.server_info"
}

fetch "$BASE_URL/${SM_TAXID}_coords.tar.gz" "$COORDS"
fetch "$BASE_URL/${SM_TAXID}_meta.tar.gz" "$META"

# -------------------------------------------------------------- 2. checks ---
log "Checking both archives (a few minutes for human)"
gzip -t "$META" || die "$META is damaged. Move it aside and run the script again."
gzip -t "$COORDS" || die "$COORDS is damaged. Move it aside and run the script again."
log "Computing SHA-256"
COORDS_SHA="$($SHA "$COORDS" | awk '{print $1}')"
META_SHA="$($SHA "$META" | awk '{print $1}')"

# The provenance record is this copy's fingerprint. A different fingerprint means
# the archives changed, and that needs a human.
if [ -f "$PROVENANCE" ] && ! { grep -q "$COORDS_SHA" "$PROVENANCE" && grep -q "$META_SHA" "$PROVENANCE"; }; then
  die "$PROVENANCE records different archives from the ones here.
       Nothing was unpacked. Check which copy is the right one before going on."
fi

# ---------------------------------------------- 3. choose and unpack models ---
count_pdb() { find "$1" -maxdepth 1 -name '*.pdb' | wc -l | tr -d ' '; }

rm -rf "$TMP"
mkdir -p "$TMP"
tar -xzf "$META" -C "$TMP" SWISS-MODEL_Repository/INDEX \
  || die "$META holds no SWISS-MODEL_Repository/INDEX."
INDEX_HEADER="$(grep '^#' "$TMP/SWISS-MODEL_Repository/INDEX" || true)"

log "Choosing one model per protein from the INDEX"
Rscript "$SCRIPT_DIR/src/swissmodel_select.R" "$TMP/SWISS-MODEL_Repository/INDEX" "$TMP/selected_models.tsv"
N_SELECTED="$(( $(wc -l < "$TMP/selected_models.tsv") - 1 ))"

if [ -d "$STRUCT_DIR" ] && [ "$(count_pdb "$STRUCT_DIR")" -gt 0 ]; then
  [ -f "$PROVENANCE" ] \
    || die "$STRUCT_DIR already holds structures, with no record of which archives they came from.
       Nothing was changed. Move that folder aside, or give a different OUT_DIR."
  log "Structures are already unpacked in $STRUCT_DIR (made from these archives); not unpacking again"
else
  # every archive member in a chosen folder whose file name holds the chosen model's coordinate ID
  # (the models are named <from>_<to>_<template>_<coordinate ID>.pdb)
  log "Listing the coordinate archive"
  tar -tzf "$COORDS" > "$TMP/members.txt"
  awk -F'\t' 'NR == FNR {if (FNR > 1) want[$2] = want[$2] " " $3; next}
              {n = split($0, p, "/"); dir = substr($0, 1, length($0) - length(p[n]) - 1)
               if (dir in want) {m = split(want[dir], c, " ")
                 for (i = 1; i <= m; i++) if (index(p[n], c[i]) > 0) {print; break}}}' \
    "$TMP/selected_models.tsv" "$TMP/members.txt" > "$TMP/to_unpack.txt"
  log "Unpacking $(wc -l < "$TMP/to_unpack.txt" | tr -d ' ') of $(wc -l < "$TMP/members.txt" | tr -d ' ') archive members"
  tar -xzf "$COORDS" -C "$TMP" -T "$TMP/to_unpack.txt"

  # swissmodel_sort.R makes the final choice from what was unpacked, into a folder renamed at the end,
  # so an interrupted run never leaves a half-filled structures_raw behind
  log "Copying the chosen models with swissmodel_sort.R"
  mkdir -p "$TMP/structures_raw"
  (cd "$TMP" && Rscript "$SCRIPT_DIR/swissmodel_sort.R" > sort.log 2>&1) \
    || die "swissmodel_sort.R failed, see $TMP/sort.log"
  N_COPIED="$(count_pdb "$TMP/structures_raw")"
  N_UNPACKED="$(wc -l < "$TMP/to_unpack.txt" | tr -d ' ')"
  [ "$N_COPIED" = "$N_UNPACKED" ] \
    || die "swissmodel_sort.R copied $N_COPIED models but $N_UNPACKED were unpacked for it. Left in $TMP."
  N_BAD="$(find "$TMP/structures_raw" -maxdepth 1 -name '*.pdb' -exec sh -c '
    for f; do grep -q -m1 "^ATOM" "$f" || echo "$f"; done' sh {} + | wc -l | tr -d ' ')"
  [ "$N_BAD" -eq 0 ] || die "$N_BAD copied models have no atoms. Left in $TMP."
  if [ -d "$STRUCT_DIR" ]; then
    find "$TMP/structures_raw" -maxdepth 1 -name '*.pdb' -exec mv -t "$STRUCT_DIR" {} +
  else
    mv "$TMP/structures_raw" "$STRUCT_DIR"
  fi
fi

# one fingerprint for the whole set: the SHA-256 of the sorted list of each file's SHA-256
N_FILES="$(count_pdb "$STRUCT_DIR")"
FINGERPRINT="$(cd "$STRUCT_DIR" && find . -maxdepth 1 -name '*.pdb' -printf '%f\n' \
  | LC_ALL=C sort | xargs $SHA | $SHA | awk '{print $1}')"
rm -rf "$TMP"

# ---------------------------------------------------------- 4. provenance ---
if [ -f "$PROVENANCE" ]; then
  log "Provenance record already exists and matches these archives; left as it is"
else
  server_line() { [ -f "$1.server_info" ] && echo "  server: $(cat "$1.server_info") (size, file date)" || echo "  server: not recorded (archive was already here)"; }
  cat > "$PROVENANCE" <<EOF
SWISS-MODEL Repository models: provenance record
Written by $(basename "${BASH_SOURCE[0]}") on $(date -u '+%Y-%m-%d %H:%M UTC')

data_version: SWISS-MODEL Repository $(echo "$INDEX_HEADER" | awk -F'\t' '$1 == "#Created:" {print $2}')

Archives (taxonomy ID $SM_TAXID)
  $BASE_URL/${SM_TAXID}_coords.tar.gz
    $(wc -c < "$COORDS" | tr -d ' ') bytes, SHA-256 $COORDS_SHA, file date here $(date -u -r "$COORDS" '+%Y-%m-%d %H:%M UTC')
  $(server_line "$COORDS")
  $BASE_URL/${SM_TAXID}_meta.tar.gz
    $(wc -c < "$META" | tr -d ' ') bytes, SHA-256 $META_SHA, file date here $(date -u -r "$META" '+%Y-%m-%d %H:%M UTC')
  $(server_line "$META")

INDEX header
$(echo "$INDEX_HEADER" | sed 's/^/  /')

Models chosen: one per protein by swissmodel_sort.R (sequence identity >= 60, best QMEANDisCo global score)
  $N_SELECTED proteins have a model of sequence identity >= 60; $N_FILES models copied to structures_raw/
  structures_raw fingerprint (SHA-256 of the sorted per-file SHA-256 list): $FINGERPRINT
EOF
fi
log "Done: $N_FILES structures in $STRUCT_DIR"
