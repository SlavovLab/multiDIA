---
layout: default
title: Home
nav_order: 1
description: "Scalable, proteoform-resolved proteomics using multiDIA for neurodegenerative disease"
permalink: /
---
{% assign v = site.time | date: '%s' %}

# **multiDIA**
{: .fs-9 }

Scalable, proteoform-resolved proteomics using multiDIA for neurodegenerative disease
{: .fs-6 .fw-300 }

[Figures]({{ site.baseurl }}/figures-page/){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [Reproduce]({{ site.baseurl }}/reproduce/){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 } [GitHub Repository](https://github.com/SlavovLab/multiDIA){: .btn .fs-5 .mb-4 .mb-md-0 }

![Figure 1]({{ site.baseurl }}/figures/fig1.png?v={{ v }}){: width="90%" }

## Aim

Trypsin, the standard protease of bottom-up proteomics, sees only part of each
protein. multiDIA digests the same tissue with three proteases, **Glu-C, Lys-C
and trypsin**, which cut after different residues and so make different
peptides from one protein. Each digest is acquired by dia-PASEF and searched
separately, and the results are pooled.

## Study

Dorsolateral striatum from 12 donors, 6 with Parkinson's disease (PD) and 6
controls, was digested with each protease, acquired by dia-PASEF on a timsTOF
Ultra (65 min, 60 windows) and searched with Spectronaut 21.1 directDIA at 1% FDR.

## Findings

- **Depth.** Pooling the digests adds identifications at every level and covers
  the same proteins more deeply: 4.0K protein groups reach ≥ 50% sequence
  coverage pooled, against 2.6K with trypsin alone.
- **Proteoforms.** Counting only peptides absent from the canonical sequence,
  pooled digests give 757 proteoforms a junction-spanning peptide, against 584
  with trypsin alone. Glu-C covers all of α-synuclein, so its peptides can be
  compared along the protein in every donor. Four of 1,217 tested proteoforms
  (MAP4, LRRFIP1, MRTFB, SPTB) change differently from their canonical form in PD.
- **Phosphosites.** Without enrichment, each added protease contributes
  phosphosites the others miss: 1,382 of 3,649 sites are placed only by Lys-C or
  Glu-C.

## Preprint

Orhun Kok, Corinna Friedrich, Kirsten Altig, Derek Oakley, Bradley T. Hyman, Mahlon Collins, Nikolai Slavov.  [*Scalable, proteoform-resolved proteomics using multiDIA for neurodegenerative disease*](https://www.biorxiv.org/content/10.64898/2026.10.16.756994v1), *bioRxiv*, doi: [10.64898/2026.10.16.756994](https://doi.org/10.64898/2026.10.16.756994)  2026.

Raw data, search results and Spectronaut settings: MassIVE
[MSV000103483](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000103483).

## Contact

[Slavov Laboratory](https://slavovlab.net), Northeastern University.
