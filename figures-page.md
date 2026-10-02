---
layout: default
title: Figures
nav_order: 2
permalink: /figures-page/
---

# Figures
{: .no_toc }

Each figure is available as [SVG](https://github.com/SlavovLab/multiDIA/tree/master/figures) and PNG.

1. TOC
{:toc}

## Figure 1 | Workflow and depth

![Figure 1]({{ site.baseurl }}/figures/fig1.png)

**(a)** Three digests of the same tissue, searched separately and pooled. On
ALDH1A1 no single digest covers more than 58% of the sequence; together they
cover 85%. **(b)** Pooling adds precursors, peptides and protein groups per
sample. **(c)** Protein groups ranked by sequence coverage: 4.0K reach ≥ 50%
pooled, against 2.6K with trypsin alone. **(d)** Coverage of the 3,924 protein
groups detected by all three digests: the same proteins are covered more deeply.

[PNG]({{ site.baseurl }}/figures/fig1.png) · [SVG]({{ site.baseurl }}/figures/fig1.svg)

## Figure 2 | Isoform-resolved quantification

![Figure 2]({{ site.baseurl }}/figures/fig2.png)

**(a)** RUFY3 isoform Q7L099-3, seen only by Lys-C and Glu-C peptides inside its
110-residue insert.
**(b)** The per-peptide model behind each isoform called in Supp. Fig. 1: every
peptide measurement in every patient, minus that peptide's fitted baseline, for the
canonical's and the isoform's peptides in control and LBD patients; black lines, the
model's fitted levels. Δ is the isoform's LBD effect minus the canonical's.
**(c)** Isoforms of PD-implicated genes reached as Lys-C and Glu-C are added to
trypsin.
**(d)** Pooled digests cover more of the sequence that distinguishes isoforms, and
more junctions are spanned by a single peptide.

[PNG]({{ site.baseurl }}/figures/fig2.png) · [SVG]({{ site.baseurl }}/figures/fig2.svg)

## Figure 3 | Phosphosites across digests

![Figure 3]({{ site.baseurl }}/figures/fig3.png)

**(a)** The 3,649 sites by the digests that place them, and by the number of
peptides containing each site. **(b)** All sites by protease, coloured by the best
precursor PEP; about 82% are seen by a single protease. **(c)** Phosphosites on
PD-associated proteins and on tau (2N4R), including the p-tau epitopes T181,
S202/T205, T217 and T231. Sites come from unenriched tissue; Spectronaut
reports only phosphosites localised with probability ≥ 0.75.

[PNG]({{ site.baseurl }}/figures/fig3.png) · [SVG]({{ site.baseurl }}/figures/fig3.svg)

## Supplementary Figure 1 | Isoform-specific differential abundance

![Supplementary Figure 1]({{ site.baseurl }}/figures/supp1_isoform_volcano.png)

Each isoform's change in LBD relative to its canonical form (Δ), from a per-peptide
linear model on every measurement (log2 quantity ~ peptide + LBD + LBD × isoform;
runs normalised to the peptides every run of a digest quantified). Isoform and
canonical peptides each span ≥ 3 LBD and ≥ 3 control patients; red, the 4 of 1,119
at Benjamini–Hochberg q ≤ 0.05 that also hold with any one patient left out
(Fig. 2b).

[PNG]({{ site.baseurl }}/figures/supp1_isoform_volcano.png) · [SVG]({{ site.baseurl }}/figures/supp1_isoform_volcano.svg)
