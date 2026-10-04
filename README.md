# multiDIA

Analysis and figures for the multiDIA manuscript: Spectronaut directDIA reports
from three parallel digests (trypsin, Lys-C, Glu-C) of one human PD brain cohort.

## Layout

```
data/         inputs; nothing here is produced by the code
derived/      tables built from data/ and read by the figures
analysis/     all code, and its tests
figures/      the paper's figures: panels, whole figures and PNGs
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
| `derived/da_iso/` | discriminating-region and junction coverage of every isoform (Fig. 2d) |
| `derived/mods/` | the phospho scan and its site table (Fig. 3) |
| `derived/multiDIA.sdrf.tsv` | the SDRF-Proteomics sample table for the PRIDE submission, one row per raw file (`prep_sdrf.py`) |

## Code

Each file in `analysis/` is named for what it produces:

| prefix | role | files |
|---|---|---|
| `fig…_` | draws that figure's panels | `fig1a_workflow`, `fig1b_depth`, `fig1cd_coverage`, `fig2a_isoform_strip`, `fig2b_ed1_isoform_da` (Fig. 2b and the Extended Data Fig. 1 volcano), `fig2cd_isoform_coverage`, `fig3a_phospho_sites`, `fig3bc_phospho_atlas` |
| `prep_` | builds `derived/` and the reports | `prep_counts`, `prep_phospho`, `prep_sdrf`, `prep_parquet` (Spectronaut TSV export to the parquet reports) |
| `lib_` | shared by the figures | `lib_report` (reports, run names, metadata), `lib_fasta`, `lib_palette`, `lib_svg`, `lib_compose` (stacking panels) |

`audit.py` checks `figures/` and `test_units.py` is the unit tests.

## Rebuild

With the environment from `pyproject.toml` (`uv sync`), from the repository root:

```bash
./build.sh                                  # derived/ and every figure, then the checks
.venv/bin/python analysis/test_units.py     # unit tests, no data needed
.venv/bin/python analysis/audit.py          # whole figures and PNGs current, one font
```
