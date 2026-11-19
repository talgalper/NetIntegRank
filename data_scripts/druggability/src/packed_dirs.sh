#!/bin/bash
# Packing and unpacking of the large working folders (sourced by the drivers and by pack_results.sh).
#
# A folder <parent>/<name>/ is packed into <parent>/<name>.tar.zst (zstd), or <name>.tar.gz when zstd is not
# installed. Packing checks the archive against the folder (file count, then every file's contents) before the
# folder is removed. Unpacking restores the folder, checks it against the archive, then removes the archive.

# pack_dir <parent> <name> <threads>: pack and check, then remove the folder. Prints the archive's name.
pack_dir() {
  local parent="$1" name="$2" threads="${3:-1}" archive comp n_files n_members
  [ -d "$parent/$name" ] || return 0
  if command -v zstd > /dev/null; then
    archive="$name.tar.zst"; comp=(-I "zstd -3 -q -T$threads")
  else
    archive="$name.tar.gz"; comp=(-z)
  fi
  rm -f "$parent/$name".tar.zst "$parent/$name".tar.gz "$parent/$archive.part"
  tar -C "$parent" "${comp[@]}" -cf "$parent/$archive.part" "$name"
  n_files=$(find "$parent/$name" | wc -l)
  n_members=$(tar -C "$parent" "${comp[@]}" -tf "$parent/$archive.part" | wc -l)
  if [ "$n_files" -ne "$n_members" ]; then
    echo "packing $parent/$name: the archive holds $n_members entries but the folder $n_files; folder kept" >&2
    rm -f "$parent/$archive.part"; return 1
  fi
  if ! tar -C "$parent" "${comp[@]}" -df "$parent/$archive.part" > /dev/null; then
    echo "packing $parent/$name: the archive does not match the folder; folder kept" >&2
    rm -f "$parent/$archive.part"; return 1
  fi
  mv "$parent/$archive.part" "$parent/$archive"
  rm -rf "${parent:?}/${name:?}"
  echo "$archive"
}

# unpack_dir <parent> <name>: if <name>.tar.zst or <name>.tar.gz is there, restore the folder from it.
unpack_dir() {
  local parent="$1" name="$2" archive comp
  if [ -f "$parent/$name.tar.zst" ]; then
    archive="$name.tar.zst"; comp=(-I zstd)
    command -v zstd > /dev/null || { echo "$parent/$archive needs zstd to unpack" >&2; return 1; }
  elif [ -f "$parent/$name.tar.gz" ]; then
    archive="$name.tar.gz"; comp=(-z)
  else
    return 0
  fi
  if [ -d "$parent/$name" ] && [ -n "$(ls -A "$parent/$name")" ]; then
    echo "$parent holds both $name/ and $archive; move one aside" >&2; return 1
  fi
  echo "Unpacking $parent/$archive"
  tar -C "$parent" "${comp[@]}" -xf "$parent/$archive"
  tar -C "$parent" "${comp[@]}" -df "$parent/$archive" > /dev/null \
    || { echo "$parent/$name does not match $archive after unpacking" >&2; return 1; }
  rm -f "$parent/$archive"
}
