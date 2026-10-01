#!/usr/bin/env python3
"""Check p-value calibration by re-running the screen on shuffled labels.

    python3 extra_perm_null.py
"""
import collections
import glob
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig2a_supp1_diagnostic_peptides as dp
from lib_fasta import gene_map
from lib_fasta import read_fasta
from lib_report import read_metadata
FASTA = os.environ.get('FASTA', 'data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta')
seqs = read_fasta(FASTA)
genes = gene_map(FASTA)
cond = read_metadata(os.environ.get('METADATA', 'data/metadata.xlsx'))
idx, med, bygene = dp.read_reports(
    sorted(glob.glob(os.environ.get('REPORTS', 'data/search/*-60min-Phospho.parquet'))), seqs, genes)
pats=sorted(cond); labs=[cond[s] for s in pats]
def frac(c, unit):
    rows=[r for r in dp.screen(idx,med,bygene,seqs,genes,c,"LBD","Control",3,2,2,unit)
          if r["p"]==r["p"]]
    qs=[r["q"] for r in rows if r["q"]==r["q"]]
    return (100*sum(1 for r in rows if r["p"]<0.05)/len(rows),
            sum(1 for q in qs if q<0.10), min(qs))
rng=random.Random(20260915)
for unit in ("vs-canonical","patients"):
    real=frac(cond, unit)
    print(f"\n{unit}:  REAL labels  p<0.05 {real[0]:.1f}%   q<0.10 {real[1]}   "
          f"best q {real[2]:.1e}", flush=True)
    for i in range(6):
        sh=labs[:]; rng.shuffle(sh)
        c=dict(zip(pats,sh))
        f=frac(c, unit)
        print(f"   perm {i+1}: p<0.05 {f[0]:5.1f}%   q<0.10 {f[1]:3d}   "
              f"best q {f[2]:.1e}   ({collections.Counter(sh)})", flush=True)
