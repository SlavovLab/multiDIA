---
layout: default
title: Home
nav_order: 1
description: "multiDIA: multi-protease DIA proteomics of human Lewy body disease brain"
permalink: /
---

# **multiDIA**
{: .fs-9 }

Multi-protease data-independent acquisition of human Lewy body disease brain
{: .fs-6 .fw-300 }

[Figures]({{ site.baseurl }}/figures-page/){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [Reproduce]({{ site.baseurl }}/reproduce/){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 } [GitHub Repository](https://github.com/SlavovLab/multiDIA){: .btn .fs-5 .mb-4 .mb-md-0 }

![Figure 1]({{ site.baseurl }}/figures/fig1.png){: width="90%" }

## Aim

Trypsin, the standard protease of bottom-up proteomics, sees only part of each
protein. multiDIA digests the same tissue with three proteases, **Glu-C, Lys-C
and trypsin**, which cut after different residues and so make different
peptides from one protein. Each digest is acquired by dia-PASEF and searched
separately, and the results are pooled.

## Study

Brain tissue from 12 donors, 6 with Lewy body disease (LBD) and 6 Lewy-body-free
controls, was digested with each protease, acquired by dia-PASEF (65 min, 60
windows) and searched with Spectronaut 21 directDIA at 1% FDR.

## Findings

- **Depth.** Pooling the digests adds identifications at every level and covers
  the same proteins more deeply: 4.0K protein groups reach ≥ 50% sequence
  coverage pooled, against 2.6K with trypsin alone.
- **Isoforms.** Pooled digests reach more of the sequence that tells isoforms
  apart. Four of 1,217 tested isoforms (MAP4, LRRFIP1, MRTFB, SPTB) change
  differently from their canonical form in LBD.
- **Phosphosites.** Without enrichment, each added protease contributes
  phosphosites the others miss: 1,382 of 3,649 sites are placed only by Lys-C or
  Glu-C.

## Contact

[Slavov Laboratory](https://slavovlab.net), Northeastern University.
