#!/usr/bin/env bash
# Downloads the monarch-initiative/phenopacket-store release archive
# (≈17 MB, ≈9.6k real published case-report phenopackets with ground-truth
# OMIM diagnoses) into data/phenopacket_store/.
#
# Re-run safe: skips download if the archive already exists; overwrite by
# deleting data/phenopacket_store/all_phenopackets.zip first.

set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$HERE/data/phenopacket_store"
URL="https://github.com/monarch-initiative/phenopacket-store/releases/latest/download/all_phenopackets.zip"

mkdir -p "$DEST"
if [[ -f "$DEST/all_phenopackets.zip" ]]; then
  echo "[skip] $DEST/all_phenopackets.zip already exists"
else
  echo "[fetch] $URL"
  curl -sSL -o "$DEST/all_phenopackets.zip" "$URL"
fi

# Unpack into versioned subdir (e.g. data/phenopacket_store/0.1.26/<gene>/*.json)
if ! find "$DEST" -mindepth 2 -maxdepth 2 -type d -name '*.*' | grep -q .; then
  echo "[unpack] $DEST/all_phenopackets.zip"
  (cd "$DEST" && unzip -q -o all_phenopackets.zip)
fi

echo "[done] phenopacket-store ready at $DEST"
find "$DEST" -mindepth 2 -maxdepth 2 -type d -name '*.*' | head -1
echo "  $(find "$DEST" -name '*.json' -type f | wc -l | tr -d ' ') phenopackets"
