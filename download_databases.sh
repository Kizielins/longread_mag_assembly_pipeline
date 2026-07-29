#!/usr/bin/env bash
## Download databases required by pipeline
set -eo pipefail

echo "=== Pipeline database download ==="
echo "This script downloads the CheckM2, GUNC, GTDB-Tk (pinned to release226 /"
echo "R226), and hostile (host-read removal) reference databases. Total size:"
echo "~150+ GB. Make sure you have sufficient disk space."
echo ""

## CheckM2 database (~3 GB) - required

echo "[1/4] Downloading CheckM2 database..."
CHECKM2_DB="${1:-$HOME/databases/checkm2}"
mkdir -p "$CHECKM2_DB"
# Check if dmnd file already present- if yes, skip
CHECKM2_DMND=$(find "$CHECKM2_DB" -maxdepth 2 -name 'uniref100.KO.1.dmnd' -print -quit 2>/dev/null)
if [ -n "$CHECKM2_DMND" ]; then
    echo "CheckM2 database already present at $CHECKM2_DMND -- skipping."
else
    conda run -n checkm2 checkm2 database --download --path "$CHECKM2_DB"
    CHECKM2_DMND=$(find "$CHECKM2_DB" -maxdepth 2 -name 'uniref100.KO.1.dmnd' -print -quit 2>/dev/null)
    echo "CheckM2 database downloaded to: $CHECKM2_DB"
fi
echo "config/config.yaml's checkm2_db default already points here unless you passed a custom path above."

## GUNC database (~13 GB) - required unless --skip_gunc

echo ""
echo "[2/4] Downloading GUNC database..."
GUNC_DB="${2:-$HOME/databases/gunc}"
mkdir -p "$GUNC_DB"
GUNC_GZ="$GUNC_DB/gunc_db_progenomes2.1.dmnd.gz"
GUNC_DMND="$GUNC_DB/gunc_db_progenomes2.1.dmnd"
GUNC_GZ_MD5="bc93a855e0760aad5c4e5f2d0e26da46"
GUNC_DMND_MD5="447c9330056b02f29f30fe81fe4af4eb"
if [ -f "$GUNC_DMND" ] && [ "$(md5sum "$GUNC_DMND" | awk '{print $1}')" = "$GUNC_DMND_MD5" ]; then
    echo "GUNC database already present at $GUNC_DB -- skipping."
elif curl -f --progress-bar -o "$GUNC_GZ" "https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz" \
    && [ "$(md5sum "$GUNC_GZ" | awk '{print $1}')" = "$GUNC_GZ_MD5" ]; then
    gunzip -f "$GUNC_GZ"
    if [ "$(md5sum "$GUNC_DMND" | awk '{print $1}')" != "$GUNC_DMND_MD5" ]; then
        echo "GUNC database md5 mismatch after decompression from EMBL mirror -- falling back to gunc's own downloader."
        rm -f "$GUNC_DMND"
        conda run -n pipeline gunc download_db "$GUNC_DB"
    fi
else
    echo "EMBL mirror unreachable or md5 mismatch -- falling back to gunc's own downloader (Zenodo, slower)."
    rm -f "$GUNC_GZ"
    conda run -n pipeline gunc download_db "$GUNC_DB"
fi
echo "GUNC database downloaded to: $GUNC_DB"
echo "config/config.yaml's gunc_db default already points here unless you passed a custom path above"

## GTDB-Tk database (~130 GB) - required for taxonomy
## Pinned to release226 (R226) for reproducibility

echo ""
echo "[3/4] Downloading GTDB-Tk database (release226 / R226)..."
GTDB_DB="${3:-$HOME/databases/gtdbtk}"
mkdir -p "$GTDB_DB"
GTDB_DONE_MARKER="$GTDB_DB/.download_complete"
if [ -f "$GTDB_DONE_MARKER" ]; then
    echo "GTDB-Tk database already present at $GTDB_DB -- skipping."
else
    wget -c --progress=bar:force:noscroll https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz -O "$GTDB_DB/gtdbtk_data.tar.gz"
    tar xzf "$GTDB_DB/gtdbtk_data.tar.gz" -C "$GTDB_DB" --strip-components=1
    rm -f "$GTDB_DB/gtdbtk_data.tar.gz"
    touch "$GTDB_DONE_MARKER"
    echo "GTDB-Tk database (R226) downloaded to: $GTDB_DB"
fi
echo "config/config.yaml's gtdbtk_data_path default already points here unless you passed a custom path above."

## hostile reference index (human host-read removal) - required unless --skip_host_removal

echo ""
echo "[4/4] Pre-fetching hostile's default human reference index (Minimap2 only)..."
conda run -n pipeline hostile index fetch -n human-t2t-hla -m

echo ""
echo "=== All databases downloaded successfully ==="
echo ""
echo "The pipeline reads all three DB paths straight from config/config.yaml"
