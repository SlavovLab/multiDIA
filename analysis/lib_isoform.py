#!/usr/bin/env python3
"""Canonical against isoform: the differences between their sequences, and the peptides each gene's groups report."""

import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ORDER = ["GluC", "LysC", "Trypsin"]
ISO = re.compile(r"-\d+$")

# Residues of agreement needed to keep two differing blocks apart.
MIN_EQUAL = 3


def align(c, i):
    """-> (common-prefix length, non-equal difflib opcodes of the trimmed cores)."""
    p = 0
    while p < min(len(c), len(i)) and c[p] == i[p]:
        p += 1
    s = 0
    while s < min(len(c), len(i)) - p and c[len(c) - 1 - s] == i[len(i) - 1 - s]:
        s += 1
    cc, ii = c[p:len(c) - s], i[p:len(i) - s]
    if not cc and not ii:
        return p, []
    if not ii:
        return p, [("delete", 0, len(cc), 0, 0)]
    if not cc:
        return p, [("insert", 0, 0, 0, len(ii))]
    import difflib
    ops = [o for o in difflib.SequenceMatcher(None, cc, ii, autojunk=False)
           .get_opcodes() if o[0] != "equal"]
    return p, ops


def merge(p, ops, min_equal=MIN_EQUAL):
    """Pool differing blocks closer than `min_equal`, shifted back by the prefix `p`."""
    merged = []
    for tag, a1, a2, b1, b2 in ops:
        if merged and a1 - merged[-1][2] < min_equal:
            merged[-1][2], merged[-1][4] = a2, b2
            merged[-1][0] = "replace"
        else:
            merged.append([tag, a1, a2, b1, b2])
    out = []
    for tag, a1, a2, b1, b2 in merged:
        if a2 > a1 and b2 > b1:
            tag = "replace"
        elif a2 > a1:
            tag = "delete"
        else:
            tag = "insert"
        out.append((tag, p + a1, p + a2, p + b1, p + b2))
    return out


def events(c, i, min_equal=MIN_EQUAL):
    """Discrete differences between canonical and isoform."""
    return merge(*align(c, i), min_equal=min_equal)


def collect_groups(paths, precursor_q, bases=None, protein_q=0.01):
    """{group: {digest: set(peptides)}} in one pass, optionally filtered."""
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    import lib_report as rp
    for r in rp.open_reports(paths, order=ORDER):
        path, dig = r.path, r.protease
        F = ["protein_groups", "peptide", "precursor_q"]
        has_pg = protein_q is not None and r.has("protein_q_run")
        if has_pg:
            F.append("protein_q_run")
        for b in r.batches(F):
            g = b["protein_groups"]
            s = b["peptide"]
            q = b["precursor_q"]
            pgq = b["protein_q_run"] if has_pg else None
            for k in range(b["_n"]):
                if not g[k] or not s[k] or q[k] is None or q[k] > precursor_q:
                    continue
                if pgq is not None and (pgq[k] is None or pgq[k] > protein_q):
                    continue
                if bases is not None:
                    accs = [a.strip() for a in g[k].split(";") if a.strip()]
                    if not ({ISO.sub("", a) for a in accs} & bases):
                        continue
                out[g[k]][dig].add(s[k])
        print(f"  scanned {os.path.basename(path)[:46]:<46s} {dig}")
    return out
