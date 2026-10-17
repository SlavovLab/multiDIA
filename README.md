# multiDIA

Analysis and figures for the preprint *Scalable, proteoform-resolved proteomics using multiDIA for neurodegenerative disease*
(Orhun Kok, Corinna Friedrich, Kirsten Altig, Derek Oakley, Bradley T. Hyman, Mahlon Collins, Nikolai Slavov, 2026): Spectronaut directDIA reports from three parallel digests (trypsin,
Lys-C, Glu-C) of one human PD brain cohort. Raw data, search results and settings are
in MassIVE, MSV000103483.

## Layout

```
data/         inputs; nothing here is produced by the code
derived/      tables built from data/ and read by the figures
analysis/     all code, and its tests
figures/      the preprint's figures and supplementary figures: panels, whole figures and PNGs
build.sh      rebuilds derived/ and figures/ from data/
index.md …    the website, https://multidia.slavovlab.net (GitHub Pages: _config.yml, CNAME, *.md)
```

| path | contents |
|---|---|
| `data/search/*-60min-Phospho.parquet` | the Spectronaut reports, one per digest, searched with variable phospho. Not in the repository (too large for GitHub); place them here to rebuild. |
| `data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta` | the search database (UniProt Swiss-Prot human with isoforms, 2024-01-01) |
| `data/metadata.xlsx` | the brain bank's neuropathology sheet; `ADRC #` is the patient id and `NPDX1` gives the condition |
| `data/uniprot/features_P10636.tsv.gz` | UniProt features of tau (Fig. 3c), downloaded once by `fig3bc_phospho_atlas.py` |
| `derived/counts/` | per-sample counts (Fig. 1b–c) |
| `derived/da_iso/` | discriminating-region and junction coverage of every proteoform, by peptides absent from the canonical (Fig. 2a) |
| `derived/peptide_levels/` | run-normalised log2 abundance per peptide and donor for SNCA, TH, DDC and MAPT, and peptides per run (Fig. 2c, Supp. Figs. 2 and 3) |
| `derived/mods/` | the phospho scan and its site table (Fig. 3) |
| `derived/multiDIA.sdrf.tsv` | the SDRF-Proteomics sample table for the PRIDE submission, one row per raw file (`prep_sdrf.py`) |

## Code

Each file in `analysis/` is named for what it produces:

| prefix | role | files |
|---|---|---|
| `fig…_`, `supp…_` | draws that figure's panels | `fig1a_workflow`, `fig1b_depth`, `fig1cd_coverage`, `fig2ab_isoform_coverage`, `fig2c_supp2_donor_levels` (Fig. 2c and Supp. Fig. 2), `fig3a_phospho_sites`, `fig3bc_phospho_atlas`, `supp1_isoform_da` (Supp. Fig. 1, and the peptide quantities Fig. 2c reads), `supp3_tau_proteoforms` |
| `prep_` | builds `derived/` and the reports | `prep_counts`, `prep_phospho`, `prep_sdrf`, `prep_parquet` (Spectronaut TSV export to the parquet reports) |
| `lib_` | shared by the figures | `lib_report` (reports, run names, metadata), `lib_fasta`, `lib_isoform` (canonical against proteoform), `lib_palette`, `lib_svg`, `lib_compose` (stacking panels) |

`audit.py` checks `figures/` and `test_units.py` is the unit tests.

## Rebuild

With the environment from `pyproject.toml` (`uv sync`), from the repository root:

```bash
./build.sh                                  # derived/ and every figure, then the checks
.venv/bin/python analysis/test_units.py     # unit tests, no data needed
.venv/bin/python analysis/audit.py          # whole figures and PNGs current, type >= 10 pt at 6.5 in wide, one font
```
