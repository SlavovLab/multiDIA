"""Reading the FASTA."""

import re


def read_fasta(path):
    """-> {accession: sequence}."""
    seqs, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                head = line[1:].split()[0]
                parts = [p for p in head.split("|") if p]
                acc = parts[1] if len(parts) >= 3 else parts[0]
                buf = []
            elif line:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def gene_map(fasta):
    """accession -> gene symbol, from the FASTA headers."""
    out = {}
    for line in open(fasta):
        if line.startswith(">"):
            h = line[1:].rstrip()
            acc = h.split("|")[1] if h.count("|") >= 2 else h.split()[0]
            m = re.search(r"GN=(\S+)", h)
            out[acc] = m.group(1) if m else acc
    return out
