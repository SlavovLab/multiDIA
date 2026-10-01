#!/usr/bin/env python3
"""Rasterise an SVG with MuPDF to preview a figure.

    python3 preview.py figures/fig1b.svg            # -> figures/fig1b.png
    python3 preview.py figures/*.svg --dpi 200 --outdir /tmp/preview
"""

import argparse
import glob
import os
import sys


def render(path, out, dpi):
    import pymupdf
    doc = pymupdf.open(path)
    pdf = pymupdf.open("pdf", doc.convert_to_pdf())
    pdf[0].get_pixmap(dpi=dpi).save(out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("svgs", nargs="+")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--outdir", default=None,
                    help="default: alongside the SVG")
    args = ap.parse_args(argv)

    paths = []
    for pat in args.svgs:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    for p in paths:
        if not os.path.exists(p):
            print(f"  missing: {p}")
            continue
        d = args.outdir or os.path.dirname(os.path.abspath(p))
        os.makedirs(d, exist_ok=True)
        out = os.path.join(d, os.path.splitext(os.path.basename(p))[0] + ".png")
        render(p, out, args.dpi)
        print(f"  {p} -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
