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
$P analysis/fig2cd_isoform_coverage.py scan "$S" --fasta $F --tsv derived/da_iso/isoform_disc_coverage.tsv
$P analysis/prep_sdrf.py "$S" --metadata $M --out derived/multiDIA.sdrf.tsv

# Figure 1
$P analysis/fig1a_workflow.py "$S" --letter a --fasta $F --protein Q96RQ3 --label MCCC1 \
    --acquisition '65 min · 60 windows' --out figures/fig1a.svg
$P analysis/fig1b_depth.py derived/counts/counts.csv --units 'precursors,peptides|protein_isoform_groups' \
    --letter b --no-xticks --no-outlier-labels --width 1215 --height 422 --out figures/fig1b.svg
$P analysis/fig1cd_coverage.py "$S" --fasta $F --rank --by-sample --sample-regex "$RX" \
    --letter c --counts derived/counts/counts.csv --out figures/fig1c.svg
$P analysis/fig1cd_coverage.py "$S" --fasta $F --beeswarm --by-sample --sample-regex "$RX" \
    --letter d --out figures/fig1d.svg

# Figure 2
$P analysis/fig2a_isoform_strip.py "$S" --fasta $F --protein P27816 --isoform P27816-3 --label MAP4 \
    --min-equal 20 --margin 130 --bar-height 22 --out "$T/strip.svg"
$P analysis/lib_compose.py "$T/strip.svg" --letter a --key Glu-C,Lys-C,Trypsin --out figures/fig2a.svg
$P analysis/fig2b_ed1_isoform_da.py "$S" --fasta $F --metadata $M \
    --model-plots "MAP4:P27816-3,MRTFB:Q9ULH7-4;Q9ULH7-5" --out "$T/models" --volcano "$T/volcano.svg"
$P analysis/lib_compose.py "$T/models/MAP4.svg" "$T/models/MRTFB.svg" --row --out "$T/row1.svg"
$P analysis/lib_compose.py "$T/row1.svg" --width 1215 --pad-top 30 --letter b --out figures/fig2b.svg
$P analysis/fig2cd_isoform_coverage.py pd "$S" --letter c --out figures/fig2c.svg
$P analysis/fig2cd_isoform_coverage.py --letter d --out figures/fig2d.svg

# Figure 3
$P analysis/fig3a_phospho_sites.py --sites derived/mods/sites_phospho_complete.tsv --letter a --out figures/fig3a.svg
$P analysis/fig3bc_phospho_atlas.py heatmap --scan derived/mods --fasta $F --out "$T/heatmap.svg"
$P analysis/lib_compose.py "$T/heatmap.svg" --width 1215 --letter b --out figures/fig3b.svg
$P analysis/fig3bc_phospho_atlas.py heatmap --scan derived/mods --fasta $F --width 350 --height 234 \
    --out "$T/pd.svg" --genes SNCA LRRK2 PINK1 PRKN "Ubiquitin:UBB,UBC,UBA52,RPS27A@P62987" \
    RAB10 RAB8A RAB12 RAB29 RAB35 TH@P07101-3 GSK3B@P49841 PLK2 "CK2:CSNK2A1,CSNK2A2,CSNK2B" \
    PARK7 SYNJ1 EIF4EBP1 "ERM:MSN,EZR,RDX" VPS35
$P analysis/fig3bc_phospho_atlas.py protein data/search/GluC-60min-Phospho.parquet \
    data/search/LysC-60min-Phospho.parquet data/search/Trypsin-60min-Phospho.parquet \
    --scan derived/mods --fasta $F --isoform P10636-8 --region N1:45-73 --region N2:74-102 \
    --region "Proline-rich:151-243" --uniprot-cache data/uniprot --rename "Tau/MAP =R" \
    --width 370 --height 234 --cell 8 --pad 3 --label-size 8 --out "$T/mapt.svg"
$P analysis/lib_compose.py "$T/pd.svg" "$T/mapt.svg" --row --out "$T/row.svg"
$P analysis/lib_compose.py "$T/row.svg" --width 1215 --letter c --out figures/fig3c.svg

# Extended Data Figure 1: the volcano written with Fig. 2b
$P analysis/lib_compose.py "$T/volcano.svg" --out figures/ed1_isoform_volcano.svg

# Whole figures and PNGs, then the checks
$P analysis/audit.py --sync
