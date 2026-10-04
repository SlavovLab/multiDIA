---
layout: default
title: Figures
nav_order: 2
permalink: /figures-page/
---
{% assign v = site.time | date: '%s' %}

# Figures
{: .no_toc }

Each figure is available as [SVG](https://github.com/SlavovLab/multiDIA/tree/master/figures) and PNG.

1. TOC
{:toc}

## Figure 1 | multiDIA workflow and sequence coverage

![Figure 1]({{ site.baseurl }}/figures/fig1.png?v={{ v }})

**(a)** Three digests of the same tissue, searched separately and pooled. The
observed coverage of MCCC1, from its gene-specific peptides across all runs:
trypsin alone covers 57% of the sequence and the three digests together 79%, each
adding residues the others miss. **(b)** Pooling adds precursors, peptides and protein groups per
sample. **(c)** Protein groups ranked by sequence coverage: 4.0K reach ≥ 50%
pooled, against 2.6K with trypsin alone. **(d)** Coverage of the 3,924 protein
groups detected by all three digests: the same proteins are covered more deeply.

[PNG]({{ site.baseurl }}/figures/fig1.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig1.svg?v={{ v }})

## Figure 2 | Isoform-resolved quantification

![Figure 2]({{ site.baseurl }}/figures/fig2.png?v={{ v }})

**(a)** MAP4 isoform P27816-3, the strongest isoform called in Extended Data Fig. 1,
against canonical P27816 on one residue scale; a thin line marks sequence the form lacks.
Two changes distinguish them: residues 2–666 replaced by 488, and 703–1152 replaced by 14,
around 36 shared residues. Bars beneath each form are its own peptides, by protease,
placed by sequence.
**(b)** Two isoforms called in Extended Data Fig. 1 (MAP4, MRTFB): log2 LBD / Control of
each peptide in each LBD donor, coloured by protease as in (a).
**(c)** Isoforms of PD-implicated genes reached as Lys-C and Glu-C are added to
trypsin.
**(d)** Pooled digests cover more of the sequence regions and junctions that distinguish
isoforms.

[PNG]({{ site.baseurl }}/figures/fig2.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig2.svg?v={{ v }})

## Figure 3 | Well localised phosphosites detected with multiDIA

![Figure 3]({{ site.baseurl }}/figures/fig3.png?v={{ v }})

**(a)** The 3,649 sites by the digests that place them, and by the number of
distinct peptides containing each site. **(b)** All sites by protease, coloured by the
best precursor PEP; about 82% are seen by a single protease. **(c)** Phosphosites on
PD-implicated proteins and on tau (2N4R), including the p-tau epitopes T181,
S202/T205, T217 and T231.

[PNG]({{ site.baseurl }}/figures/fig3.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig3.svg?v={{ v }})

## Extended Data Figure 1 | Isoform-specific differential abundance

![Extended Data Figure 1]({{ site.baseurl }}/figures/ed1_isoform_volcano.png?v={{ v }})

Each isoform's change in LBD relative to its canonical form: Δ = (isoform LBD −
control) − (canonical LBD − control). Isoforms whose canonical form was not measured
are not tested. The 4 isoforms that diverge from their canonical form at
Benjamini–Hochberg q ≤ 0.05 are highlighted in red.

[PNG]({{ site.baseurl }}/figures/ed1_isoform_volcano.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/ed1_isoform_volcano.svg?v={{ v }})
