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

**(a)** Workflow schematic: aliquots of each sample are digested separately with
Glu-C, Lys-C or trypsin; each digest is analyzed and searched separately, and the
results are pooled (top). Below, the observed coverage of MCCC1 across all runs:
trypsin alone covers 57% of the sequence and the three digests together cover 79%,
with each digest adding residues missed by the others. **(b)** Number of quantified
precursors, peptides and protein groups per sample. **(c)** Protein groups ranked by
sequence coverage: 4.0K reach ≥50% coverage with all digests pooled, versus 2.6K
with trypsin alone. **(d)** Sequence coverage of the 3,924 protein groups detected in
all three digests; the same proteins are covered more deeply by multiDIA.

[PNG]({{ site.baseurl }}/figures/fig1.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig1.svg?v={{ v }})

## Figure 2 | Proteoforms and PD-associated peptides by donor

![Figure 2]({{ site.baseurl }}/figures/fig2.png?v={{ v }})

**(a)** Proteoforms by discriminating-region coverage and with a junction-spanning
peptide. **(b)** Proteoforms of PD-associated genes; grey, not covered. a and b count
only peptides absent from the canonical sequence. **(c)** α-Synuclein (SNCA) peptide
abundance by donor, log2 relative to the control mean (each control against the
other controls). Top: one row per donor, ordered within each group by profile
similarity; each residue shows the mean of its covering peptides; grey, not
quantified. Bottom: each PD donor's peptides as lines over the residues they cover;
vertical lines mark the deletions of P37840-3 (Δ41–54) and P37840-2 (Δ103–130).

[PNG]({{ site.baseurl }}/figures/fig2.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig2.svg?v={{ v }})

## Figure 3 | Phosphosites detected without enrichment

![Figure 3]({{ site.baseurl }}/figures/fig3.png?v={{ v }})

**(a)** The 3,649 phosphosites by the digests in which they were detected (left) and
by the number of distinct peptides containing each site (right). **(b)** All
phosphosites by protease, colored by the posterior error probability (PEP) of the
best precursor; about 82% were detected by a single protease. **(c)** Modified
residues of MAPT on its 2N4R 441-residue isoform (P10636-8), one cell per protease,
colored as in b; N1 and N2 are the inserts encoded by exons 2 and 3, and R1–R4 the
microtubule-binding repeats, R2 encoded by exon 10.

[PNG]({{ site.baseurl }}/figures/fig3.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/fig3.svg?v={{ v }})

## Supplementary Figure 1 | Proteoform-specific differential abundance

![Supplementary Figure 1]({{ site.baseurl }}/figures/supp1.png?v={{ v }})

**(a)** Change in each proteoform's abundance in LBD relative to its canonical form:
Δ = (proteoform log2 LBD/control) − (canonical log2 LBD/control); proteoforms whose
canonical form was not measured were not tested, and the four that diverge from their
canonical form at Benjamini–Hochberg q ≤ 0.05 are highlighted in red.
**(b)** Two of them (MAP4, MRTFB): log2 LBD/control ratio of each peptide in each LBD
donor, colored by protease.

[PNG]({{ site.baseurl }}/figures/supp1.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/supp1.svg?v={{ v }})

## Supplementary Figure 2 | Dopaminergic enzyme peptide abundance by donor

![Supplementary Figure 2]({{ site.baseurl }}/figures/supp2.png?v={{ v }})

Tyrosine hydroxylase (TH) and DDC peptide abundance along each protein, as in the
top of Fig. 2c: one row per donor, ordered within each group by profile similarity;
each residue shows the mean log2 abundance of the peptides covering it, relative to
the control mean (each control against the other controls); grey, not quantified.
Bars above the rows mark trypsin (green), Lys-C (blue) and Glu-C (orange) peptides.

[PNG]({{ site.baseurl }}/figures/supp2.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/supp2.svg?v={{ v }})

## Supplementary Figure 3 | Peptides detected on MAPT proteoforms

![Supplementary Figure 3]({{ site.baseurl }}/figures/supp3.png?v={{ v }})

Big tau (P10636) and the six adult brain proteoforms (0N3R to 2N4R, P10636-2 to -8),
aligned on the big-tau sequence: black, residues each proteoform contains; thin line,
exons it skips (exons 2, 3, 4a, 6 and 10 marked above). Bars above big tau mark the
residues covered by each protease's tau-unique peptides. Below each proteoform, the
peptides unique to tau that occur in its sequence, colored by protease; a peptide
drawn across a gap spans the junction of the skipped exon. Big tau's peptides are
colored instead by log2 PD / control, each peptide against its own control mean as in
Fig. 2c; grey, not computed.

[PNG]({{ site.baseurl }}/figures/supp3.png?v={{ v }}) · [SVG]({{ site.baseurl }}/figures/supp3.svg?v={{ v }})
