#!/bin/sh
# Rebuild derived/ and figures/ from data/. Run from the repository root.
set -eu

P=.venv/bin/python
F=data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta
S='data/search/*-60min-Phospho.parquet'
M=data/metadata.xlsx
RX='CF[_-](?:[A-Za-z]{1,3}[_-])?(\d{4})'
T=$(mktemp -d "${TMPDIR:-/tmp}/multidia.XXXXXX")
trap 'rm -rf "$T"' EXIT
export PYTHONHASHSEED=0

# Derived tables
$P analysis/prep_counts.py "$S" --sample-regex "$RX" --min-run-count 10 --out derived/counts/
$P analysis/prep_phospho.py scan "$S" --fasta $F --outdir derived/mods
$P analysis/prep_phospho.py export --scan derived/mods --complete --out derived/mods/sites_phospho_complete.tsv
$P analysis/fig2ab_isoform_coverage.py scan "$S" --fasta $F --tsv derived/da_iso/isoform_disc_coverage.tsv
$P analysis/fig2c_supp2_donor_levels.py scan "$S" --fasta $F --genes SNCA TH DDC MAPT --out derived/peptide_levels
$P analysis/prep_sdrf.py "$S" --metadata $M --out derived/multiDIA.sdrf.tsv

# Figure 1
$P analysis/fig1a_workflow.py "$S" --letter a --fasta $F --protein Q96RQ3 --label MCCC1 \
    --acquisition '65 min · 60 windows' --out figures/fig1a.svg
$P analysis/fig1b_depth.py derived/counts/counts.csv --units 'precursors,peptides|protein_isoform_groups' \
    --letter b --width 1215 --height 355 --out figures/fig1b.svg
$P analysis/fig1cd_coverage.py "$S" --fasta $F --rank --by-sample --sample-regex "$RX" \
    --letter c --counts derived/counts/counts.csv --out figures/fig1c.svg
$P analysis/fig1cd_coverage.py "$S" --fasta $F --beeswarm --by-sample --sample-regex "$RX" \
    --letter d --out figures/fig1d.svg

# Figure 2
$P analysis/fig2ab_isoform_coverage.py --tsv derived/da_iso/isoform_disc_coverage.tsv --letter a --width 655 \
    --out figures/fig2a.svg
$P analysis/fig2ab_isoform_coverage.py pd "$S" --fasta $F --letter b --out figures/fig2b.svg
$P analysis/fig2c_supp2_donor_levels.py profile SNCA --fasta $F --metadata $M --letter c --out figures/fig2c.svg

# Figure 3
$P analysis/fig3a_phospho_sites.py --sites derived/mods/sites_phospho_complete.tsv --letter a --out figures/fig3a.svg
$P analysis/fig3bc_phospho_atlas.py heatmap --scan derived/mods --out "$T/heatmap.svg"
$P analysis/lib_compose.py "$T/heatmap.svg" --width 1215 --letter b --out figures/fig3b.svg
$P analysis/fig3bc_phospho_atlas.py protein data/search/GluC-60min-Phospho.parquet \
    data/search/LysC-60min-Phospho.parquet data/search/Trypsin-60min-Phospho.parquet \
    --scan derived/mods --fasta $F --isoform P10636-8 --region N1:45-73 --region N2:74-102 \
    --region "Proline-rich:151-243" --uniprot-cache data/uniprot --rename "Tau/MAP =R" \
    --cell 12 --pad 3 --label-size 10 --out "$T/mapt.svg"
$P analysis/lib_compose.py "$T/mapt.svg" --width 1215 --letter c --out figures/fig3c.svg

# Supplementary Figure 1
$P analysis/supp1_isoform_da.py "$S" --fasta $F --metadata $M \
    --model-plots "MAP4:P27816-3,MRTFB:Q9ULH7-4;Q9ULH7-5" --width 1215 --out "$T/models.svg" \
    --volcano "$T/volcano.svg"
$P analysis/lib_compose.py "$T/volcano.svg" --width 1215 --letter a --out figures/supp1a.svg
$P analysis/lib_compose.py "$T/models.svg" --width 1215 --letter b --out figures/supp1b.svg

# Supplementary Figures 2 and 3
$P analysis/fig2c_supp2_donor_levels.py along TH DDC --fasta $F --metadata $M --out figures/supp2.svg
$P analysis/supp3_tau_proteoforms.py "$S" --fasta $F --metadata $M --protein P10636 --exon 2:45-73 --exon 3:74-102 \
    --exon 4a:125-375 --exon 6:394-459 --exon 10:592-622 --out figures/supp3.svg

# Whole figures and PNGs, then the checks
$P analysis/audit.py --sync
