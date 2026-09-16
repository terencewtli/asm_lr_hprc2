#!/usr/bin/env python3
# Export the relevant slice of hprc2_supp.xlsx (S6: per-sample seq/coverage QC; S3: 1000G
# superpopulation secondary labels) to one flat CSV for the R/ggplot2 QC notebook to load --
# same architecture as this project's figure-notebook convention (Python exports CSV, R plots
# it), and replaces this session's own P04-computed QC now that the paper's authoritative
# Supplementary Table 6 is in hand (bioRxiv's supplementary endpoint was 429-rate-limited
# earlier; the user supplied the xlsx directly instead).

import re
from pathlib import Path

import pandas as pd

XLSX = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/reference/hprc2/hprc2_supp.xlsx')
OUT = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/results/qc/data/supp_seq_qc.csv')
OUT.parent.mkdir(parents=True, exist_ok=True)

xl = pd.ExcelFile(XLSX)

s6 = xl.parse('S6', header=1)
s6 = s6[~s6['project'].isna()]  # drops the CHM13/GRCh38 reference rows (no seq QC)

# S3: one row per superpopulation (AFR, EUR, EAS, SAS, AMR -- skip the AFR/non-AFR aggregate
# rows), with a comma-separated "Accessions" column -- explode into sample_id -> superpop.
s3 = xl.parse('S3', header=2)
s3.columns = ['secondary_label', 'primary_labels', 'accessions']
# the AFR row's label is literally "AFR; African" (not "AFR") -- caught by a first pass that
# silently dropped every AFR sample as unmapped (74/234 NaN, and AFR absent from the value
# counts entirely, which is what actually surfaced it). Match the token before ';', not the
# raw string.
clean_label = s3['secondary_label'].str.split(';').str[0].str.strip()
superpop_rows = s3[clean_label.isin(['AFR', 'EUR', 'EAS', 'SAS', 'AMR'])].assign(
    secondary_label=clean_label[clean_label.isin(['AFR', 'EUR', 'EAS', 'SAS', 'AMR'])])

pop_map = []
for _, row in superpop_rows.iterrows():
    for sample_id in re.split(r',\s*', row['accessions'].strip()):
        pop_map.append({'sample_id': sample_id, 'superpopulation': row['secondary_label']})
pop_df = pd.DataFrame(pop_map)

merged = s6.merge(pop_df, on='sample_id', how='left')
merged.to_csv(OUT, index=False)
print(f'wrote {len(merged)} rows to {OUT}')
print(merged['superpopulation'].value_counts(dropna=False))
