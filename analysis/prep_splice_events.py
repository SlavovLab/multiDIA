#!/usr/bin/env python3
"""Splice-event types for UniProt isoforms, from GENCODE exon structures.

    python3 prep_splice_events.py --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --cache data/uniprot --gencode data/gencode --out data/gencode/events.tsv
"""

import argparse
import collections
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# report order, most specific first
ORDER = ["SE", "MXE", "A5SS", "A3SS", "RI", "AF", "AL"]

NAMES = {
    "SE": "Skipped exon",
    "MXE": "Mutually exclusive exons",
    "A5SS": "Alternative 5' splice site",
    "A3SS": "Alternative 3' splice site",
    "RI": "Retained intron",
    "AF": "Alternative first exon",
    "AL": "Alternative last exon",
}

# reasons an isoform has no event
UNMAPPED = "no Ensembl transcript"
SAME = "identical exon structure"
CONFLICT = "canonical and isoform share a transcript"


def read_xrefs(path):
    """-> {isoform accession: [transcript ids]}, from UniProt's Ensembl column."""
    out = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if not header.startswith("Entry"):
            sys.exit(f"{path}: unexpected header {header[:60]!r}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[1]:
                continue
            acc = parts[0]
            for item in parts[1].split(";"):
                item = item.strip()
                if not item:
                    continue
                if "[" in item:
                    enst, tag = item.split("[", 1)
                    iso = tag.rstrip("]").strip()
                else:
                    enst, iso = item, acc
                enst = enst.strip().split(".")[0]
                if enst:
                    out[iso].append(enst)
    return dict(out)


def read_refseq_map(path):
    """-> {RefSeq id (versionless): {transcript}} from GENCODE's metadata."""
    out = collections.defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            enst = f[0].split(".")[0]
            for r in f[1:]:
                if r:
                    out[r.split(".")[0]].add(enst)
    return out


def read_translations(path):
    """-> {protein sequence: {transcript}} from GENCODE's pc_translations FASTA."""
    out = collections.defaultdict(set)
    tid, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if tid:
                    out["".join(buf)].add(tid)
                tid, buf = line[1:].split("|")[1].split(".")[0], []
            else:
                buf.append(line.strip())
    if tid:
        out["".join(buf)].add(tid)
    return out


class Router:
    """Accession -> transcripts, by the first route in ORDER that has any."""

    ORDER = ("ensembl", "translation", "refseq")

    def __init__(self, ens, seqs, translations=None, refseq=None, r2e=None):
        self.ens, self.seqs = ens, seqs
        self.translations = translations or {}
        self.refseq, self.r2e = refseq or {}, r2e or {}
        self.used = collections.Counter()

    def __call__(self, acc):
        """-> (transcripts, route). Empty set and None when nothing maps."""
        for route in self.ORDER:
            ts = self._route(acc, route)
            if ts:
                return ts, route
        return set(), None

    def _route(self, acc, route):
        if route == "ensembl":
            return set(self.ens.get(acc, ()))
        if route == "translation":
            return set(self.translations.get(self.seqs.get(acc, "\0"), ()))
        return {e for r in self.refseq.get(acc, ())
                for e in self.r2e.get(r.split(".")[0], ())}


def read_exons(path, wanted=None, feature="CDS"):
    """-> {transcript: (chrom, strand, [(start, end), ...])} sorted, 1-based."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != feature:
                continue
            attrs = f[8]
            i = attrs.find('transcript_id "')
            if i < 0:
                continue
            j = attrs.find('"', i + 15)
            tid = attrs[i + 15:j].split(".")[0]
            if wanted is not None and tid not in wanted:
                continue
            rec = out.get(tid)
            if rec is None:
                rec = out[tid] = (f[0], f[6], [])
            rec[2].append((int(f[3]), int(f[4])))
    for chrom, strand, ex in out.values():
        ex.sort()
    return out


def introns(exons):
    """-> [(start, end)] of the gaps between consecutive exons, genomic order."""
    return [(exons[k][1] + 1, exons[k + 1][0] - 1)
            for k in range(len(exons) - 1)
            if exons[k + 1][0] - exons[k][1] > 1]


def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def classify(can, iso, strand):
    """Pairwise splice events between two genomic-ascending exon lists -> set of codes."""
    found = set()
    if not can or not iso:
        return found
    plus = strand == "+"
    can_i, iso_i = introns(can), introns(iso)
    can_set, iso_set = set(can), set(iso)

    # skipped exon: wholly inside an intron of the other transcript
    def skipped(a, b_introns, b):
        return [e for e in a
                if e not in set(b)
                and any(i[0] <= e[0] and e[1] <= i[1] for i in b_introns)]

    can_only = skipped(can, iso_i, iso)
    iso_only = skipped(iso, can_i, can)
    if can_only or iso_only:
        found.add("SE")
    # mutually exclusive: each side skips a non-overlapping exon of the other
    for a in can_only:
        for b in iso_only:
            if not overlaps(a, b):
                found.add("MXE")
                break

    # retained intron: an intron covered by a single exon of the other
    for i in can_i:
        if any(e[0] <= i[0] and i[1] <= e[1] for e in iso):
            found.add("RI")
            break
    for i in iso_i:
        if any(e[0] <= i[0] and i[1] <= e[1] for e in can):
            found.add("RI")
            break

    # shifted splice site: overlapping exons sharing exactly one boundary
    for a in can:
        if a in iso_set:
            continue
        for b in iso:
            if b in can_set or not overlaps(a, b):
                continue
            same_lo, same_hi = a[0] == b[0], a[1] == b[1]
            if same_lo == same_hi:
                continue
            # genomic-low end is the acceptor on +, the donor on -
            donor_moved = same_lo if plus else same_hi
            found.add("A5SS" if donor_moved else "A3SS")

    # alternative first / last exon
    c_first, i_first = (can[0], iso[0]) if plus else (can[-1], iso[-1])
    c_last, i_last = (can[-1], iso[-1]) if plus else (can[0], iso[0])
    if not overlaps(c_first, i_first):
        found.add("AF")
    if not overlaps(c_last, i_last):
        found.add("AL")
    return found


def type_isoforms(seqs, router, gtf, base_of):
    """-> {isoform: (codes, reason)}; `reason` is set only when codes is empty."""
    pairs_of, wanted = {}, set()
    for iso in seqs:
        base = base_of(iso)
        if base == iso or base not in seqs:
            continue
        iso_ts, iso_r = router(iso)
        can_ts, can_r = router(f"{base}-1")
        if not can_ts:
            can_ts, can_r = router(base)
        pairs_of[iso] = (can_ts, iso_ts)
        if can_ts and iso_ts:
            router.used[f"{can_r}/{iso_r}"] += 1
            wanted |= can_ts | iso_ts

    exons = read_exons(gtf, wanted)
    print(f"  {len(exons):,} transcripts with exons, of {len(wanted):,} mapped")

    out = {}
    for iso, (can_ts, iso_ts) in pairs_of.items():
        best = best_pair_codes(can_ts, iso_ts, exons)
        if best is None:
            shared = can_ts & iso_ts
            out[iso] = (frozenset(),
                        CONFLICT if shared and not (can_ts - shared)
                        and not (iso_ts - shared) else UNMAPPED)
            continue
        out[iso] = (frozenset(best), "" if best else SAME)
    return out


def best_pair_codes(can_ts, iso_ts, exons):
    """-> events of the transcript pair with the fewest, or None if none comparable."""
    pairs = [(c, i) for c in sorted(can_ts) for i in sorted(iso_ts)
             if c != i                              # never self-compare
             and c in exons and i in exons
             and exons[c][0] == exons[i][0]]        # same chromosome
    if not pairs:
        return None
    best = None
    for c, i in pairs:
        _chrom, strand, cex = exons[c]
        codes = classify(cex, exons[i][2], strand)
        if best is None or len(codes) < len(best):
            best = codes
        if not best:
            break
    return best


def load(path):
    """-> {isoform: (code, ...)} from the table this module writes; untyped keep their reason."""
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            codes = tuple(c for c in r["events"].split(",") if c)
            out[r["isoform"]] = codes or ((r["note"],) if r["note"] else ())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--cache", default="data/uniprot",
                    help="where data/uniprot/ensembl.tsv.gz lives")
    ap.add_argument("--gencode", default="data/gencode",
                    help="directory holding the GENCODE GTF")
    ap.add_argument("--gtf", default=None, help="explicit path to the GTF")
    ap.add_argument("--out", default=None, help="write the per-isoform table")
    args = ap.parse_args(argv)

    from lib_fasta import read_fasta
    import prep_uniprot as uv

    seqs = read_fasta(args.fasta)
    ens = read_xrefs(uv.fetch(args.cache, fields="accession,xref_ensembl",
                              name="ensembl.tsv.gz"))
    print(f"  {len(ens):,} accessions with an Ensembl transcript")

    # fallback routes for isoforms without an Ensembl cross-reference
    translations = refseq = r2e = None
    tr = os.path.join(args.gencode, "gencode.v50.pc_translations.fa.gz")
    if os.path.exists(tr):
        translations = read_translations(tr)
        print(f"  {len(translations):,} distinct GENCODE translations")
    rs = os.path.join(args.gencode, "gencode.v50.metadata.RefSeq.gz")
    if os.path.exists(rs):
        r2e = read_refseq_map(rs)
        refseq = read_xrefs(uv.fetch(args.cache, fields="accession,xref_refseq",
                                     name="refseq.tsv.gz"))
        print(f"  {len(r2e):,} RefSeq ids in GENCODE, {len(refseq):,} "
              f"accessions with one in UniProt")
    router = Router(ens, seqs, translations, refseq, r2e)

    gtf = args.gtf
    if not gtf:
        cand = sorted(f for f in os.listdir(args.gencode) if f.endswith(".gtf.gz"))
        if not cand:
            sys.exit(f"no .gtf.gz in {args.gencode}/")
        gtf = os.path.join(args.gencode, cand[-1])
    print(f"  exons from {gtf}")

    import re as _re
    iso_re = _re.compile(r"-\d+$")
    typed = type_isoforms(seqs, router, gtf, lambda a: iso_re.sub("", a))

    tally = collections.Counter()
    for codes, reason in typed.values():
        if codes:
            for c in codes:
                tally[c] += 1
        else:
            tally[reason] += 1
    n = len(typed)
    print(f"\n{n:,} isoforms typed against their canonical transcript")
    for k in ORDER:
        if tally[k]:
            print(f"  {tally[k]:>6,}  {tally[k] / n:6.2%}  {NAMES[k]} ({k})")
    for k in (SAME, CONFLICT, UNMAPPED):
        if tally[k]:
            print(f"  {tally[k]:>6,}  {tally[k] / n:6.2%}  {k}")
    multi = sum(1 for c, _r in typed.values() if len(c) > 1)
    print(f"  {multi:,} isoforms carry more than one event type")
    print("\n  route used (canonical/isoform), for the pairs that mapped:")
    for k, v in router.used.most_common():
        print(f"    {v:>6,}  {v / max(sum(router.used.values()), 1):6.1%}  {k}")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                    exist_ok=True)
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["isoform", "events", "note"])
            for iso in sorted(typed):
                codes, reason = typed[iso]
                w.writerow([iso, ",".join(c for c in ORDER if c in codes),
                            reason])
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
