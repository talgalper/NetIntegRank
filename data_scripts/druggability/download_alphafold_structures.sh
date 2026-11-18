#!/usr/bin/env bash
#
# download_alphafold_structures.sh
#
# Downloads the AlphaFold DB bulk archive for one reference proteome (human by
# default) from EBI, checks that it arrived whole, and unpacks the PDB files.
# It repeats, step for step, what was done by hand on 2026-11-16 for AlphaFold
# DB v6.
#
# Usage
#   ./download_alphafold_structures.sh                  # newest version EBI offers
#   AF_VERSION=6 ./download_alphafold_structures.sh     # a specific version
#   OUT_DIR=/some/folder ./download_alphafold_structures.sh
#
# What it does
#   1. Works out which archive to fetch. With AF_VERSION=latest it reads EBI's
#      own list of current archives, so it follows new releases by itself.
#   2. Downloads it over several connections at once (aria2c), because EBI
#      limits the speed of each single connection. Safe to re-run: an
#      interrupted download carries on, a finished one is not fetched again.
#   3. Checks it: the size, then the built-in error check of every file inside.
#      EBI publishes no checksum, so this is what proves the archive is intact.
#   4. Unpacks the .pdb files into OUT_DIR/structures_raw. The .cif copies stay
#      inside the archive, which is kept.
#   5. Writes a provenance record (address, version, dates, size, SHA-256)
#      next to the archive.
#
# It never deletes or overwrites structures that are already there. If
# OUT_DIR/structures_raw holds structures of a different version, it stops and
# says so. Use a separate OUT_DIR for each version.
#
# Needs: Linux with GNU tar, curl, python3, gzip, sha256sum. aria2c is
# optional but strongly recommended (sudo apt install aria2); without it the
# download uses one connection and takes hours.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------- settings ---
# Each can be set here or given on the command line as shown under Usage.

# AlphaFold DB version: "latest", or a number such as 6.
AF_VERSION="${AF_VERSION:-latest}"

# UniProt reference proteome. UP000005640 is human.
AF_PROTEOME="${AF_PROTEOME:-UP000005640}"

# Where the archive and the structures_raw folder go.
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/alphafold}"

# Parallel connections for the download. aria2c allows at most 16.
CONNECTIONS="${CONNECTIONS:-16}"

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/alphafold"

# ----------------------------------------------------------------- helpers ---
log() { printf '[%s] %s\n' "$(date -u '+%H:%M:%S')" "$*"; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

for tool in curl python3 gzip sha256sum tar; do
  command -v "$tool" >/dev/null || die "'$tool' is needed but was not found."
done
tar --version 2>/dev/null | grep -q 'GNU tar' \
  || die "GNU tar is needed (this is not GNU tar). On a Mac: brew install gnu-tar, then put it first on the PATH."
[[ "$CONNECTIONS" =~ ^[0-9]+$ ]] && [ "$CONNECTIONS" -ge 1 ] || die "CONNECTIONS must be a whole number, got '$CONNECTIONS'."
[ "$CONNECTIONS" -le 16 ] || CONNECTIONS=16

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --------------------------------------------- 1. which archive to fetch ---
log "Reading EBI's list of current AlphaFold archives"
curl -fsS --retry 3 --max-time 120 -o "$WORK/metadata.json" "$BASE_URL/download_metadata.json" \
  || die "Could not read $BASE_URL/download_metadata.json"

read -r LATEST_ARCHIVE LATEST_SIZE LATEST_COUNT < <(python3 - "$WORK/metadata.json" "$AF_PROTEOME" <<'PY'
import json, sys
hits = [x for x in json.load(open(sys.argv[1])) if x.get("reference_proteome") == sys.argv[2]]
if len(hits) == 1:
    print(hits[0]["archive_name"], hits[0]["size_bytes"], hits[0]["num_predicted_structures"])
PY
) || true
[ -n "${LATEST_ARCHIVE:-}" ] || die "EBI's list has no single archive for reference proteome '$AF_PROTEOME'."
[[ "$LATEST_ARCHIVE" =~ ^(.+)_v([0-9]+)\.tar$ ]] || die "Unexpected archive name from EBI: $LATEST_ARCHIVE"
ARCHIVE_STEM="${BASH_REMATCH[1]}"
LATEST_VERSION="${BASH_REMATCH[2]}"

EXPECTED_COUNT=""   # number of structures EBI states; only known for the latest version
if [ "$AF_VERSION" = "latest" ]; then
  VERSION="$LATEST_VERSION"
else
  VERSION="${AF_VERSION#v}"
  [[ "$VERSION" =~ ^[0-9]+$ ]] || die "AF_VERSION must be 'latest' or a number such as 6, got '$AF_VERSION'."
fi
[ "$VERSION" = "$LATEST_VERSION" ] && EXPECTED_COUNT="$LATEST_COUNT"

ARCHIVE="${ARCHIVE_STEM}_v${VERSION}.tar"
URL="$BASE_URL/v${VERSION}/$ARCHIVE"

curl -sSIL --retry 3 --max-time 120 "$URL" > "$WORK/head.txt" || die "Could not reach $URL"
HTTP_STATUS="$(awk 'toupper($1) ~ /^HTTP/ {s=$2} END {print s}' "$WORK/head.txt")"
[ "$HTTP_STATUS" = "200" ] || die "EBI answered HTTP $HTTP_STATUS for $URL. Is v$VERSION a real version? The newest is v$LATEST_VERSION."
EXPECTED_SIZE="$(awk 'tolower($1) == "content-length:" {gsub("\r", "", $2); s=$2} END {print s}' "$WORK/head.txt")"
SERVER_DATE="$(awk 'tolower($1) == "last-modified:" {sub(/^[^:]*: */, ""); gsub("\r", ""); d=$0} END {print d}' "$WORK/head.txt")"
[[ "$EXPECTED_SIZE" =~ ^[0-9]+$ ]] || die "EBI did not report a size for $URL"
if [ -n "$EXPECTED_COUNT" ] && [ "$LATEST_SIZE" != "$EXPECTED_SIZE" ]; then
  log "NOTE: EBI's list says $LATEST_SIZE bytes but the file on the server is $EXPECTED_SIZE bytes. Using the server's figure."
fi

# The version exists, so the output folder can be created now.
mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
ARCHIVE_PATH="$OUT_DIR/$ARCHIVE"
STRUCT_DIR="$OUT_DIR/structures_raw"
PROVENANCE="$OUT_DIR/${ARCHIVE%.tar}.provenance.txt"

echo
echo "  AlphaFold DB version : v$VERSION  (asked for: $AF_VERSION; newest available: v$LATEST_VERSION)"
echo "  Archive              : $ARCHIVE"
echo "  Address              : $URL"
echo "  Size                 : $EXPECTED_SIZE bytes"
echo "  Structures (EBI)     : ${EXPECTED_COUNT:-not stated for older versions}"
echo "  Saving into          : $OUT_DIR"
echo

# Structures from another version must never be mixed with this one. Only .pdb
# files count here; a stray log file in the folder is ignored and left alone.
if [ -d "$STRUCT_DIR" ]; then
  OTHER="$(find "$STRUCT_DIR" -maxdepth 1 -name '*.pdb' ! -name "AF-*-model_v${VERSION}.pdb" -print -quit)"
  [ -z "$OTHER" ] || die "$STRUCT_DIR already holds structures that are not v$VERSION (for example $(basename "$OTHER")).
       Nothing was changed. Move that folder aside, or give a different OUT_DIR."
fi

# ------------------------------------------------------------ 2. download ---
have_complete_archive() {
  [ -f "$ARCHIVE_PATH" ] && [ ! -e "$ARCHIVE_PATH.aria2" ] \
    && [ "$(wc -c < "$ARCHIVE_PATH")" = "$EXPECTED_SIZE" ]
}

if have_complete_archive; then
  log "Archive is already here at full size; not downloading again"
else
  FREE_KB="$(df -Pk "$OUT_DIR" | awk 'NR == 2 {print $4}')"
  NEED_KB=$(( EXPECTED_SIZE * 3 / 1024 ))
  [ "$FREE_KB" -ge "$NEED_KB" ] || die "Not enough free disk in $OUT_DIR: about $(( NEED_KB / 1048576 )) GB needed (archive plus unpacked files), $(( FREE_KB / 1048576 )) GB free."
  if command -v aria2c >/dev/null; then
    log "Downloading over $CONNECTIONS connections with aria2c"
    aria2c -c -x "$CONNECTIONS" -s "$CONNECTIONS" --file-allocation=none \
           --summary-interval=0 --console-log-level=warn \
           -d "$OUT_DIR" -o "$ARCHIVE" "$URL" \
      || die "The download did not finish. Run the script again to carry on from where it stopped."
  else
    log "aria2c is not installed: downloading over ONE connection, which takes hours for the human archive."
    log "To make it fast: sudo apt install aria2, then run this script again."
    curl -fL --retry 5 -C - -o "$ARCHIVE_PATH" "$URL" \
      || die "The download did not finish. Run the script again to carry on from where it stopped."
  fi
  have_complete_archive \
    || die "The download looks unfinished (wrong size, or a .aria2 helper file is still there). Run the script again."
fi

# -------------------------------------------------------------- 3. checks ---
# Size alone proves nothing for a download made in pieces, so test every file
# inside. --ignore-zeros (-i) is needed: without it tar stops quietly at a
# blank stretch and reports success.
log "Checking every file inside the archive (a few minutes for the human archive)"
tar -i -xf "$ARCHIVE_PATH" --to-command='gzip -t' \
  || die "The archive is damaged: at least one file inside failed its error check.
       Move $ARCHIVE_PATH aside and run the script again to download a fresh copy."

tar -i -tf "$ARCHIVE_PATH" > "$WORK/members.txt"
N_PDB="$(awk '/\.pdb\.gz$/ {n++} END {print n + 0}' "$WORK/members.txt")"
N_CIF="$(awk '/\.cif\.gz$/ {n++} END {print n + 0}' "$WORK/members.txt")"
N_OTHER="$(awk '!/\.(pdb|cif)\.gz$/ {n++} END {print n + 0}' "$WORK/members.txt")"
N_UNSAFE="$(awk '/\// || /\.\./ || /^-/ {n++} END {print n + 0}' "$WORK/members.txt")"
N_DUP="$(sort "$WORK/members.txt" | uniq -d | wc -l)"
[ "$N_UNSAFE" -eq 0 ] || die "The archive holds $N_UNSAFE file names with a folder or '..' in them, which is not expected. Not unpacking."
[ "$N_DUP" -eq 0 ] || die "The archive holds $N_DUP duplicated file names, which is not expected. Not unpacking."
[ "$N_PDB" -gt 0 ] || die "The archive holds no .pdb.gz files."
if [ -n "$EXPECTED_COUNT" ] && [ "$N_PDB" != "$EXPECTED_COUNT" ]; then
  die "The archive holds $N_PDB structures but EBI states $EXPECTED_COUNT."
fi
log "Archive is intact: $N_PDB .pdb.gz, $N_CIF .cif.gz, $N_OTHER other"

log "Computing SHA-256"
SHA256="$(sha256sum "$ARCHIVE_PATH" | awk '{print $1}')"

# The provenance record is this copy's fingerprint. A different fingerprint for
# the same archive name means the file changed, and that needs a human.
if [ -f "$PROVENANCE" ] && ! grep -q "$SHA256" "$PROVENANCE"; then
  die "$PROVENANCE records a different SHA-256 from this archive ($SHA256).
       Nothing was unpacked. Check which copy is the right one before going on."
fi

# -------------------------------------------------------------- 4. unpack ---
count_pdb() { find "$1" -maxdepth 1 -name '*.pdb' | wc -l; }

if [ -d "$STRUCT_DIR" ] && [ "$(count_pdb "$STRUCT_DIR")" -gt 0 ]; then
  [ "$(count_pdb "$STRUCT_DIR")" = "$N_PDB" ] \
    || die "$STRUCT_DIR holds $(count_pdb "$STRUCT_DIR") v$VERSION structures but the archive holds $N_PDB.
       Nothing was changed. Move that folder aside and run the script again."
  log "Structures are already unpacked in $STRUCT_DIR; not unpacking again"
else
  # Unpack beside the final folder and rename at the end, so an interrupted
  # run can never leave a half-filled structures_raw behind.
  PART="$STRUCT_DIR.incomplete"
  mkdir -p "$PART"
  log "Unpacking $N_PDB PDB files"
  tar -xf "$ARCHIVE_PATH" -C "$PART" --wildcards '*.pdb.gz'
  find "$PART" -maxdepth 1 -name '*.pdb.gz' -exec gunzip -f {} +
  [ "$(count_pdb "$PART")" = "$N_PDB" ] || die "Unpacked $(count_pdb "$PART") files, expected $N_PDB. Left in $PART."
  if [ -d "$STRUCT_DIR" ]; then
    # The folder exists but holds no structures yet: move them in.
    find "$PART" -maxdepth 1 -name '*.pdb' -exec mv -t "$STRUCT_DIR" {} +
    rmdir "$PART" 2>/dev/null || log "NOTE: $PART still holds other files and was left in place"
  else
    mv "$PART" "$STRUCT_DIR"
  fi
fi

# ------------------------------------------------- 5. check what was unpacked ---
log "Checking the unpacked files"
N_FILES="$(count_pdb "$STRUCT_DIR")"
N_LEFT_GZ="$(find "$STRUCT_DIR" -maxdepth 1 -name '*.gz' | wc -l)"
N_EMPTY="$(find "$STRUCT_DIR" -maxdepth 1 -name '*.pdb' -size 0 | wc -l)"
N_BAD="$(find "$STRUCT_DIR" -maxdepth 1 -name '*.pdb' -exec sh -c '
  for f; do
    grep -q -m1 "^ATOM" "$f" && [ "$(tail -n 1 "$f" | cut -c1-3)" = "END" ] || echo "$f"
  done' sh {} + | wc -l)"
[ "$N_FILES" = "$N_PDB" ] || die "$STRUCT_DIR holds $N_FILES .pdb files, expected $N_PDB."
[ "$N_LEFT_GZ" -eq 0 ] || die "$N_LEFT_GZ files in $STRUCT_DIR are still compressed."
[ "$N_EMPTY" -eq 0 ] || die "$N_EMPTY files in $STRUCT_DIR are empty."
[ "$N_BAD" -eq 0 ] || die "$N_BAD files in $STRUCT_DIR have no atoms or do not end with END."

UNPACKED_BYTES="$(find "$STRUCT_DIR" -maxdepth 1 -name '*.pdb' -printf '%s\n' | awk '{s += $1} END {printf "%d", s}')"
find "$STRUCT_DIR" -maxdepth 1 -name '*.pdb' -printf '%f\n' \
  | sed -E 's/-model_v[0-9]+\.pdb$//; s/^AF-//' > "$WORK/ids.txt"
N_PROTEINS="$(sed -E 's/-F[0-9]+$//' "$WORK/ids.txt" | sort -u | wc -l)"
N_FRAGMENTED="$(grep -E -- '-F2$' "$WORK/ids.txt" | wc -l || true)"

# One fingerprint for the whole unpacked set: the SHA-256 of the sorted list of
# each file's SHA-256. Two machines with the same fingerprint have the same files.
FINGERPRINT="$(cd "$STRUCT_DIR" && find . -maxdepth 1 -name '*.pdb' -printf '%f\n' \
  | LC_ALL=C sort | xargs sha256sum | sha256sum | awk '{print $1}')"

# ---------------------------------------------------------- 6. provenance ---
if [ -f "$PROVENANCE" ]; then
  log "Provenance record already exists and matches this archive; left as it is"
else
  cat > "$PROVENANCE" <<EOF
AlphaFold DB bulk archive: provenance record
Written by $(basename "${BASH_SOURCE[0]}") on $(date -u '+%Y-%m-%d %H:%M UTC')

Archive                 $ARCHIVE
Source address          $URL
AlphaFold DB version    v$VERSION (asked for: $AF_VERSION; newest available at the time: v$LATEST_VERSION)
Reference proteome      $AF_PROTEOME
File date on server     $SERVER_DATE
Downloaded              $(date -u -r "$ARCHIVE_PATH" '+%Y-%m-%d %H:%M UTC')
Size                    $EXPECTED_SIZE bytes
SHA-256                 $SHA256
Structures, per EBI     ${EXPECTED_COUNT:-not stated for older versions}
Files in the archive    $N_PDB .pdb.gz, $N_CIF .cif.gz, $N_OTHER other
Integrity check         passed: tar -i -xf $ARCHIVE --to-command='gzip -t'

Unpacked to             $STRUCT_DIR
Unpacked files          $N_FILES .pdb, $UNPACKED_BYTES bytes
Proteins                $N_PROTEINS, of which $N_FRAGMENTED are split into fragments (F1, F2, ...)
Unpacked fingerprint    $FINGERPRINT
  (reproduce inside the folder with:
   find . -maxdepth 1 -name '*.pdb' -printf '%f\n' | LC_ALL=C sort | xargs sha256sum | sha256sum)

Tools                   $(tar --version | head -n 1); $(gzip --version | head -n 1); $(aria2c --version 2>/dev/null | head -n 1 || echo 'aria2c not installed'); $(curl --version | head -n 1 | cut -d' ' -f1-2)
EOF
  log "Provenance record written"
fi

echo
echo "  Done. AlphaFold DB v$VERSION, $AF_PROTEOME"
echo "  Archive      : $ARCHIVE_PATH"
echo "  SHA-256      : $SHA256"
echo "  Structures   : $STRUCT_DIR"
echo "                 $N_FILES files, $N_PROTEINS proteins ($N_FRAGMENTED split into fragments), $UNPACKED_BYTES bytes"
echo "  Fingerprint  : $FINGERPRINT"
echo "  Provenance   : $PROVENANCE"
echo
