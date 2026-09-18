#!/bin/bash
# Sync the curated docs/code from the working directory into the git mirror
# (github/asm_lr_hprc2/), then commit. Run this at the end of a session — it's the single
# place new content should flow from working dir -> mirror; never hand-edit files inside
# github/asm_lr_hprc2/ directly (exception: JOURNAL*.md, which live only there), or the two will silently diverge again (see
# JOURNAL.md 2026-09-08 and asm_lr's md/20260907.results.summary.md incident for why this
# matters).
#
# Usage: bash scripts/github/sync_to_github.sh ["commit message"]
# With no message, stages+shows the diff but does not commit (dry-run-ish; you still push
# yourself, this script never pushes).

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
GH=$PROJDIR/github/asm_lr_hprc2
MSG=${1:-}

# Curated set only — mirrors asm_lr's convention of keeping heavy/binary/generated files
# (logs/, modbed/, assemblies/, tmp/) out of the trimmed git mirror.
#
# `data/` USED to be safe to rsync wholesale back when it only held the small sample manifest
# (see 2026-09-08 entries). It no longer is: a since-reorg now nests modbed/, chains/,
# assemblies/, het_snps/, hmmflagger/ under data/ (~2.8TB combined as of 2026-09-16) -- rsyncing
# that into a git mirror would be a disaster. Dropped from this script for that reason (caught
# 2026-09-16 before ever running it post-reorg). If a genuinely small curated file needs to be
# mirrored again (e.g. tsv/meta/hprc2_sample_manifest.tsv), add it back explicitly by path, not
# as a blanket directory rsync -- the whole point of this script is that a directory can quietly
# start holding something huge without anyone updating the mirror logic to match.
cp "$PROJDIR/README.md" "$GH/README.md"
# JOURNAL.md / JOURNAL.archive.md live ONLY in the mirror ($GH) since 2026-09-17 (user request):
# edit them there directly; they are not copied from the working dir.
mkdir -p "$GH/docs" "$GH/scripts/qsub/test" "$GH/scripts/harmonize_hg38/pilot"
rsync -av --delete "$PROJDIR/docs/" "$GH/docs/"
rsync -av --delete --exclude='__pycache__/' "$PROJDIR/scripts/qsub/" "$GH/scripts/qsub/"
# harmonize_hg38/pilot/ code only — excludes tmp_chr*/ (per-chrom scratch, tens of GB each)
# and any other non-.py output the pilot scripts write alongside themselves.
rsync -av --delete --include='*.py' --exclude='tmp_chr*/' --exclude='__pycache__/' \
    --exclude='test/' --exclude='*' \
    "$PROJDIR/scripts/harmonize_hg38/pilot/" "$GH/scripts/harmonize_hg38/pilot/"

# PMD calling + bigWig export: code only (.py/.sh), recursive.
for d in call_pmds bigwig meth_bins asm ucsc; do
    mkdir -p "$GH/scripts/$d"
    rsync -av --delete --include='*/' --include='*.py' --include='*.sh' \
        --exclude='__pycache__/' --exclude='*' --prune-empty-dirs \
        "$PROJDIR/scripts/$d/" "$GH/scripts/$d/"
done

# notebooks/ — the analysis code behind every figure. Executed notebooks carry their outputs,
# which is the point: a future session can read what a cell produced without rerunning it.
# notebooks/old/ (retired pilots) IS mirrored, deliberately: it is small and it is the provenance
# for results that were superseded. Only checkpoints are excluded.
mkdir -p "$GH/notebooks"
rsync -av --delete --include='*/' --include='*.ipynb' --exclude='.ipynb_checkpoints/' \
    --exclude='*' --prune-empty-dirs \
    "$PROJDIR/notebooks/" "$GH/notebooks/"

# figures/ (small PNGs) are mirrored so RESULTS.md can embed them on GitHub
mkdir -p "$GH/figures"
rsync -av --delete --include='*/' --include='*.png' --exclude='*' --prune-empty-dirs \
    "$PROJDIR/figures/" "$GH/figures/"

# the sync script itself, so the mirror contains the tool that defines it
mkdir -p "$GH/scripts/github"
cp "$PROJDIR/scripts/github/sync_to_github.sh" "$GH/scripts/github/sync_to_github.sh"

cd "$GH"
git add -A
if git diff --cached --quiet; then
    echo "$(date): no changes to sync"
    exit 0
fi

echo "=== staged diff summary ==="
git diff --cached --stat

if [ -z "$MSG" ]; then
    echo ""
    echo "No commit message given — changes staged but not committed."
    echo "Re-run as: bash scripts/sync_to_github.sh \"your message\""
    exit 0
fi

git commit -m "$MSG"
echo "$(date): committed. Push manually when ready (git push) — this script never pushes."
