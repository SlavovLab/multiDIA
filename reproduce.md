---
layout: default
title: Reproduce
nav_order: 3
permalink: /reproduce/
---

# Reproduce the figures

The [repository](https://github.com/SlavovLab/multiDIA) holds the code, the
derived tables and every figure. `build.sh` rebuilds all of them from the inputs
in `data/`.

## Layout

| path | contents |
|---|---|
| `data/` | inputs: the Spectronaut reports, the FASTA, sample metadata, GENCODE and UniProt annotation |
| `derived/` | tables built from `data/` and read by the figures |
| `analysis/` | the code; each `fig…_` module draws the figure in its name |
| `figures/` | the figures as SVG, with PNGs of the whole figures |
| `build.sh` | rebuilds `derived/` and `figures/` from `data/` |

## Inputs not in the repository

The Spectronaut reports (`data/search/*-60min-Phospho.parquet`, one per
protease) are too large for GitHub; place them in `data/search/` to rebuild.
`data/gencode/gencode.v50.annotation.gtf.gz` comes from
[GENCODE](https://www.gencodegenes.org/human/release_50.html) and is needed only
to recompute the splice events.

## Build

With Python 3.14 and [uv](https://docs.astral.sh/uv/), from the repository root:

```bash
uv sync
./build.sh
```

`.venv/bin/python analysis/test_units.py` runs the unit tests, which need no data.
