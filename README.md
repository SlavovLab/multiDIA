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
docs.zip      archived notes, captions, references and tool READMEs
```

| path | contents |
|---|---|
| `data/search/*-60min-Phospho.parquet` | the Spectronaut reports, one per digest, searched with variable phospho. Not in the repository (too large for GitHub); place them here to rebuild. |
| `data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta` | the search database (UniProt Swiss-Prot human with isoforms, 2024-01-01) |
| `data/metadata.xlsx` | the brain bank's neuropathology sheet; `ADRC #` is the patient id and `NPDX1` gives the condition |
| `data/gencode/`, `data/uniprot/` | GENCODE v50 annotation with its splice events (`events.tsv`), and cached UniProt features. `gencode.v50.annotation.gtf.gz` is not in the repository; download it from GENCODE to rerun `prep_splice_events.py`. |
| `data/da_panel.txt` | the marker panel for `extra_differential.py --panel` |
| `derived/counts/` | per-sample counts (Fig. 1a–c) |
| `derived/da_iso/` | discriminating-region coverage (Fig. 2c) |
| `derived/mods/` | the phospho scan and its site table (Fig. 3) |

## Code

Each file in `analysis/` is named for what it produces:

| prefix | role | files |
|---|---|---|
| `fig…_` | draws that figure's panels | `fig1a_workflow`, `fig1b_depth`, `fig1cd_coverage`, `fig2a_supp1_diagnostic_peptides`, `fig2b_isoform_strip`, `fig2cd_isoform_coverage`, `fig3a_phospho_sites`, `fig3bc_phospho_atlas` |
| `prep_` | builds `derived/` or the `data/` caches | `prep_counts`, `prep_phospho`, `prep_splice_events`, `prep_uniprot`, `prep_parquet` |
| `lib_` | shared by the figures | `lib_report` (reports, run names, metadata), `lib_fasta`, `lib_palette`, `lib_svg`, `lib_compose` (stacking panels) |
| `extra_` | exploratory; no final figure | `extra_differential`, `extra_isoform_da`, `extra_isoform_unique`, `extra_changed_proteoforms`, `extra_digest_upset`, `extra_reported_groups`, `extra_perm_null`, `extra_why_multienzyme` |

`audit.py` checks `figures/`, `test_units.py` is the unit tests, `make_synthetic.py`
writes a test report and `preview.py` renders an SVG to PNG.

## Rebuild

With the environment from `pyproject.toml` (`uv sync`), from the repository root:

```bash
./build.sh                                  # derived/ and every figure, then the checks
.venv/bin/python analysis/test_units.py     # unit tests, no data needed
.venv/bin/python analysis/audit.py          # whole figures and PNGs current, one font
```
