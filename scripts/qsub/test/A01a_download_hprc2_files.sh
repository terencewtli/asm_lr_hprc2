#!/bin/bash
#$ -N A01a_download_hprc2_files
#$ -cwd
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 10
#$ -t 1-229:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01a_download_hprc2_files.$JOB_ID.$TASK_ID
#$ -e /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01a_download_hprc2_files.$JOB_ID.$TASK_ID

# One array task per row of the sample manifest (row = SGE_TASK_ID, 1-indexed, header excluded).
# Downloads that donor's hap1/hap2 ONT modbed (+.tbi index) from the public hprc-epigenome S3
# bucket, and hap1/hap2 assembly FASTA from human-pangenomics S3, where resolved.
#
# No `aws` CLI on Hoffman2 as of 2026-09-08 (checked) — both buckets are public/unsigned, so
# plain `curl` against the https://<bucket>.s3.amazonaws.com/<key> path works without credentials.
#
# Re-running is safe: every download is skip-if-exists (see D03_download_hprc_assemblies.sh in
# ../../../asm_lr/scripts/download/ for the precedent this follows). Never deletes anything.
#
# Local test (single sample, no qsub): see test/A01a_download_hprc2_files.sh (ID=1 hardcoded).

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/data/hprc2_sample_manifest.tsv
MODBED_DIR=$PROJDIR/modbed
ASSEMBLY_DIR=$PROJDIR/assemblies
EPIGENOME_BASE=https://hprc-epigenome.s3.amazonaws.com/samples

mkdir -p "$MODBED_DIR" "$ASSEMBLY_DIR" "$PROJDIR/logs"

ID=1  # hardcoded for testing; production uses SGE_TASK_ID
# ID=$SGE_TASK_ID

# Manifest columns (1-indexed, tab-separated):
# 1 sample_id  2 population  3 superpopulation  4 ont_hap1_modbed_mb  5 ont_hap2_modbed_mb
# 6 has_pacbio_methylc  7 has_fiberseq  8 has_hic  9 has_isoseq_expression
# 10 has_harmonized_sup5_basecall  11 had_guppy_raw  12 had_dorado06_raw
# 13 assembly_naming  14 assembly_hap1_url  15 assembly_hap2_url
# 16 in_asm_lr_18donor_cohort  17 in_asm_lr_30donor_plan
ROW=$(awk -F'\t' -v t="$ID" 'NR==t+1 {print}' "$MANIFEST")
if [ -z "$ROW" ]; then
    echo "ERROR: no manifest row for task $ID (manifest has $(($(wc -l < "$MANIFEST") - 1)) data rows)"
    exit 1
fi

SAMPLE=$(echo "$ROW" | cut -f1)
ASM_NAMING=$(echo "$ROW" | cut -f13)
ASM_HAP1_URL=$(echo "$ROW" | cut -f14)
ASM_HAP2_URL=$(echo "$ROW" | cut -f15)

echo "$(date): A01a task $ID — $SAMPLE"

download_if_missing () {
    local url=$1
    local out=$2
    if [ -z "$url" ]; then
        echo "$(date): no URL given, skipping $(basename "$out")"
        return 0
    fi
    if [ -s "$out" ]; then
        echo "$(date): exists, skipping: $out"
        return 0
    fi
    echo "$(date): downloading $(basename "$out")"
    # curl to a temp name first so a killed/failed job never leaves a truncated file
    # sitting at the final path looking like a completed download on the next run.
    time curl -sS --fail -o "${out}.partial" "$url"
    mv "${out}.partial" "$out"
    echo "$(date): done — $(du -sh "$out" | cut -f1)  $out"
}

# --- ONT hap1/hap2 modbed + index ---
for HAP in 1 2; do
    MODBED_OUT=$MODBED_DIR/${SAMPLE}_hap${HAP}.modbed.gz
    download_if_missing "$EPIGENOME_BASE/$SAMPLE/methylation.ONT.hap${HAP}.modbed.gz" "$MODBED_OUT"
    download_if_missing "$EPIGENOME_BASE/$SAMPLE/methylation.ONT.hap${HAP}.modbed.gz.tbi" "${MODBED_OUT}.tbi"
done

# --- hap1/hap2 assembly FASTA (skipped cleanly if this sample's path was never resolved — see
#     docs/data_sources.md §4 for the ~27 samples with no release2 assembly deposited yet) ---
if [ "$ASM_NAMING" = "UNRESOLVED" ] || [ -z "$ASM_HAP1_URL" ]; then
    echo "$(date): $SAMPLE has no resolved assembly URL (assembly_naming=$ASM_NAMING) — skipping assembly download"
else
    download_if_missing "$ASM_HAP1_URL" "$ASSEMBLY_DIR/${SAMPLE}_hap1.fa.gz"
    download_if_missing "$ASM_HAP2_URL" "$ASSEMBLY_DIR/${SAMPLE}_hap2.fa.gz"
    # best-effort md5 sidecar check — HPRC release2 assemblies are usually deposited with a
    # matching .md5 file; not present for every sample/release, so a missing sidecar is a
    # skip, not a failure.
    for HAP in 1 2; do
        FA=$ASSEMBLY_DIR/${SAMPLE}_hap${HAP}.fa.gz
        URL_VAR=ASM_HAP${HAP}_URL
        MD5_URL="${!URL_VAR}.md5"
        if curl -sS --fail -o "${FA}.md5.tmp" "$MD5_URL" 2>/dev/null; then
            EXPECTED=$(awk '{print $1}' "${FA}.md5.tmp")
            ACTUAL=$(md5sum "$FA" | awk '{print $1}')
            if [ "$EXPECTED" = "$ACTUAL" ]; then
                echo "$(date): hap${HAP} md5 OK"
                mv "${FA}.md5.tmp" "${FA}.md5"
            else
                echo "WARNING: hap${HAP} md5 MISMATCH (expected $EXPECTED, got $ACTUAL) — not deleting, flagging only"
                mv "${FA}.md5.tmp" "${FA}.md5.MISMATCH"
            fi
        else
            rm -f "${FA}.md5.tmp"
            echo "$(date): no .md5 sidecar available for hap${HAP}, skipping checksum"
        fi
    done
fi

echo "$(date): $SAMPLE (task $ID) complete"
