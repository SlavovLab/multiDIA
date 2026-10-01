#!/usr/bin/env python3
"""UniProt's VAR_SEQ account of how each isoform differs from its canonical form.

    # download once (cached), then report on the cohort's FASTA
    python3 prep_uniprot.py --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --cache data/uniprot

    # just the fetch
    python3 prep_uniprot.py --fetch-only --cache data/uniprot
"""

import argparse
import collections
import gzip
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,ft_var_seq,cc_alternative_products"
QUERY = "(reviewed:true) AND (organism_id:9606)"

# Notes contain "; ", so anchor on the keyword rather than split.
FEATURE = re.compile(r"VAR_SEQ\s+(\d+)(?:\.\.(\d+))?;\s*(.*?)(?=VAR_SEQ\s+\d|$)",
                     re.S)
NOTE = re.compile(r'/note="((?:[^"\\]|\\.)*)"')
VSP_ID = re.compile(r'/id="(VSP_\d+)"')
# "IsoId=P37840-3; Sequence=VSP_006363, VSP_006364;" -- IsoId may list several.
ISO_BLOCK = re.compile(r"IsoId=([^;]+);\s*Sequence=([^;]+);")
EVENT = re.compile(r"Event=([^;]+);")


def fetch(cache_dir, force=False, query=QUERY, page=500, tries=4,
          fields=None, name="varseq.tsv.gz"):
    """Download a human reviewed UniProt table once and keep it."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, name)
    if os.path.exists(path) and not force:
        print(f"  cached {path} ({os.path.getsize(path):,} bytes)")
        return path

    part = path + ".part"
    url = SEARCH + "?" + urllib.parse.urlencode(
        {"query": query, "format": "tsv", "fields": fields or FIELDS,
         "size": page})
    rows = pages = 0
    with gzip.open(part, "wt") as fh:
        while url:
            for attempt in range(1, tries + 1):
                try:
                    with urllib.request.urlopen(url, timeout=120) as r:
                        body = r.read().decode()
                        link = r.headers.get("Link", "")
                    break
                except Exception as e:                       # noqa: BLE001
                    if attempt == tries:
                        os.unlink(part)
                        sys.exit(f"\n  giving up on page {pages + 1}: "
                                 f"{type(e).__name__}: {e}")
                    print(f"\n  page {pages + 1} attempt {attempt} failed "
                          f"({type(e).__name__}), retrying", file=sys.stderr)
                    time.sleep(2 * attempt)
            lines = body.splitlines()
            if pages == 0:
                fh.write(lines[0] + "\n")                    # header once
            for line in lines[1:]:
                fh.write(line + "\n")
                rows += 1
            pages += 1
            print(f"\r  {rows:,} entries in {pages} pages", end="", file=sys.stderr)
            m = re.search(r'<([^>]+)>;\s*rel="next"', link)
            url = m.group(1) if m else None
    os.replace(part, path)
    print(f"\r  wrote {path}: {rows:,} entries, {pages} pages, "
          f"{os.path.getsize(path):,} bytes")
    return path


def parse(path):
    """-> {isoform accession: record}, one record per named isoform."""
    out = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if not header.startswith("Entry"):
            sys.exit(f"{path}: unexpected header {header[:60]!r}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, var_seq, alt = parts[0], parts[1], parts[2]
            if not alt:
                continue

            # VSP id -> (start, end, kind, note)
            feats = {}
            for m in FEATURE.finditer(var_seq):
                start, end, body = m.group(1), m.group(2), m.group(3)
                vid = VSP_ID.search(body)
                note = NOTE.search(body)
                if not vid:
                    continue
                text = note.group(1) if note else ""
                # "Missing (in isoform 2-5)" vs "MAEP... -> MLRA... (in isoform X)"
                kind = "delete" if text.startswith("Missing") else "replace"
                feats[vid.group(1)] = (int(start), int(end or start), kind, text)

            ev = EVENT.search(alt)
            events = ([e.strip() for e in ev.group(1).split(",")] if ev else [])
            for m in ISO_BLOCK.finditer(alt):
                seq = m.group(2).strip()
                for iso in (i.strip() for i in m.group(1).split(",")):
                    if not iso:
                        continue
                    if seq == "Displayed":
                        applied, status = [], "canonical"
                    elif seq.startswith("VSP_"):
                        applied = [v.strip() for v in seq.split(",")]
                        status = "described"
                    else:
                        # "Not described" / "External": no sequence from these features.
                        applied, status = [], seq.lower().replace(" ", "-")
                    out[iso] = {
                        "base": acc,
                        "status": status,
                        "mechanism": events,
                        "vsps": [(v, ) + feats[v] for v in applied if v in feats],
                        "missing_vsps": [v for v in applied if v not in feats],
                    }
    return out


def classify(rec):
    if rec["status"] == "canonical":
        return "canonical"
    if rec["status"] != "described":
        return f"sequence {rec['status']}"
    n = len(rec["vsps"])
    if n == 0:
        return "no features resolved"
    if n > 1:
        return f"{min(n, 6)}{'+' if n >= 6 else ''} separate changes"
    _vid, start, end, kind, _note = rec["vsps"][0]
    return kind, start, end


def place(kind, start, end, length):
    at_n = start <= 3
    at_c = end >= length - 2 if length else False
    where = "N-terminal" if at_n else ("C-terminal" if at_c else "internal")
    verb = "segment absent" if kind == "delete" else "segment replaced"
    return f"{where} {verb}"


def alignments(isos, seqs):
    """-> [(accession, canonical length, offset, opcodes)], difflib per isoform."""
    import fig2b_isoform_strip as ie
    out = []
    for acc, rec in sorted(isos.items()):
        if rec["status"] != "described" or not rec["vsps"]:
            continue
        base = rec["base"]
        if base not in seqs or acc not in seqs:
            continue
        p, ops = ie.align(seqs[base], seqs[acc])
        out.append((acc, len(seqs[base]), len(seqs[acc]), p, ops))
    return out


def compare(isos, aligned, min_equal=None):
    """Score the difflib decomposition against the curation, isoform by isoform."""
    import fig2b_isoform_strip as ie
    if min_equal is None:
        min_equal = ie.MIN_EQUAL

    n_agree = lab_agree = total = 0
    delta = collections.Counter()
    confusion = collections.Counter()
    for acc, clen, ilen, p, ops in aligned:
        rec = isos[acc]
        total += 1
        want_n = len(rec["vsps"])
        got = ie.merge(p, ops, min_equal)
        delta[len(got) - want_n] += 1
        if len(got) == want_n:
            n_agree += 1

        want = classify(rec)
        want = want if isinstance(want, str) else place(*want, clen)
        have = ie.label(got, clen, ilen)
        if want == have:
            lab_agree += 1
        confusion[(want, have)] += 1
    return {"total": total, "n_agree": n_agree, "lab_agree": lab_agree,
            "delta": delta, "confusion": confusion, "min_equal": min_equal}


def report_compare(c):
    t = max(c["total"], 1)
    print(f"\ndifflib vs curation, min_equal={c['min_equal']}, "
          f"{c['total']:,} isoforms with both")
    print(f"  same NUMBER of changes  {c['n_agree']:>7,}  {c['n_agree'] / t:6.1%}")
    print(f"  same LABEL              {c['lab_agree']:>7,}  {c['lab_agree'] / t:6.1%}")
    print("\n  difflib count minus curated count")
    for d, v in sorted(c["delta"].items()):
        print(f"    {d:+3d}  {v:>7,}  {v / t:6.1%}")
    print("\n  most common disagreements (curated -> difflib)")
    wrong = [(k, v) for k, v in c["confusion"].items() if k[0] != k[1]]
    for (want, have), v in sorted(wrong, key=lambda kv: -kv[1])[:12]:
        print(f"    {v:>6,}  {want:<28s} -> {have}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default="data/uniprot")
    ap.add_argument("--fasta", default=None,
                    help="restrict the report to accessions in this FASTA")
    ap.add_argument("--force", action="store_true", help="re-download")
    ap.add_argument("--fetch-only", action="store_true")
    ap.add_argument("--compare", action="store_true",
                    help="score fig2b_isoform_strip.py's difflib call against the curation")
    ap.add_argument("--sweep", default=None,
                    help="comma-separated min_equal values to score, e.g. 0,3,5,10")
    args = ap.parse_args(argv)

    path = fetch(args.cache, args.force)
    if args.fetch_only:
        return 0

    recs = parse(path)
    isos = {a: r for a, r in recs.items() if r["status"] != "canonical"}
    print(f"\n{len(recs):,} isoform records, {len(isos):,} non-canonical")

    if args.fasta:
        from lib_fasta import read_fasta
        seqs = read_fasta(args.fasta)
        have = {a for a in isos if a in seqs}
        print(f"{len(have):,} of them are in {os.path.basename(args.fasta)}"
              f"  ({len(set(isos) - set(seqs)):,} annotated but absent, "
              f"{len([a for a in seqs if '-' in a and a not in recs]):,} in the "
              f"FASTA but unannotated)")
        isos = {a: isos[a] for a in have}
        lengths = {a: len(seqs[r['base']]) for a, r in isos.items()
                   if r["base"] in seqs}
    else:
        lengths = {}

    mech = collections.Counter()
    for r in isos.values():
        mech[", ".join(r["mechanism"]) or "<none>"] += 1
    print("\nmechanism (as annotated per entry)")
    for k, v in mech.most_common():
        print(f"  {v:>7,}  {k}")

    kinds = collections.Counter()
    for a, r in isos.items():
        c = classify(r)
        kinds[c if isinstance(c, str)
              else place(*c, lengths.get(a, 0))] += 1
    print("\nchange type (curated)")
    for k, v in kinds.most_common():
        print(f"  {v:>7,}  {k}")

    unresolved = sum(1 for r in isos.values() if r["missing_vsps"])
    if unresolved:
        print(f"\n{unresolved:,} isoforms reference a VSP with no feature line")

    if args.compare or args.sweep:
        if not args.fasta:
            sys.exit("--compare needs --fasta")
        print("\n  aligning ...")
        aligned = alignments(isos, seqs)
        if args.sweep:
            print(f"\nmin_equal sweep over {len(aligned):,} isoforms")
            print(f"  {'min_equal':>9}  {'same n':>8}  {'same label':>11}  "
                  f"{'undercount':>11}  {'overcount':>10}")
            for me in (int(v) for v in args.sweep.split(",")):
                c = compare(isos, aligned, me)
                t = max(c["total"], 1)
                under = sum(v for d, v in c["delta"].items() if d < 0)
                over = sum(v for d, v in c["delta"].items() if d > 0)
                print(f"  {me:>9}  {c['n_agree'] / t:>7.1%}  {c['lab_agree'] / t:>10.1%}"
                      f"  {under / t:>10.1%}  {over / t:>9.1%}")
        else:
            report_compare(compare(isos, aligned))
    return 0


if __name__ == "__main__":
    sys.exit(main())
