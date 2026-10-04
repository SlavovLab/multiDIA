#!/usr/bin/env python3
"""Check figures/: whole figures match their panels, PNGs are current, type is
readable at print width, one font stack.

    python3 audit.py [--sync]     # --sync rebuilds the whole figures and PNGs
"""

import argparse
import contextlib
import filecmp
import glob
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES = os.path.join(ROOT, "figures")


def rel(p):
    return os.path.relpath(p, ROOT)


PANEL = re.compile(r"fig(\d+)[a-z]\.svg")


def combined_groups():
    """-> {figN.svg: [figNa.svg, figNb.svg, ...]}, panels in letter order."""
    groups = {}
    for p in sorted(glob.glob(os.path.join(FIGURES, "fig*.svg"))):
        m = PANEL.fullmatch(os.path.basename(p))
        if m:
            whole = os.path.join(FIGURES, f"fig{m.group(1)}.svg")
            groups.setdefault(whole, []).append(p)
    return groups


# Whole figures not laid out as one stack: rows top to bottom, each row's
# space-separated columns side by side, each column's letters stacked.
LAYOUT = {"fig1.svg": ["a", "b", "c d"], "fig2.svg": ["abd c"]}
COLUMN_GAP = 24.0


def build_whole(compose, panels, name, out, tmp):
    """Compose one whole figure from its panels, per `LAYOUT`."""
    rows = LAYOUT.get(name)
    if not rows:
        compose.compose(panels, out, natural=True, reletter=True)
        return
    by_letter = {os.path.basename(p)[-5]: p for p in panels}
    if sorted(by_letter) != sorted("".join(rows).replace(" ", "")):
        sys.exit(f"panels {sorted(by_letter)} do not match the layout {rows}")
    row_svgs = []
    for i, row in enumerate(rows):
        cols = []
        for j, letters in enumerate(row.split()):
            col = os.path.join(tmp, f"{name}.r{i}c{j}.svg")
            compose.compose([by_letter[x] for x in letters], col, natural=True,
                            reletter=True)
            cols.append(col)
        if len(cols) == 1:
            row_svgs.append(cols[0])
            continue
        row_svg = os.path.join(tmp, f"{name}.r{i}.svg")
        compose.compose_row(cols, row_svg, gap=COLUMN_GAP, natural=True)
        row_svgs.append(row_svg)
    if len(row_svgs) == 1:
        shutil.copyfile(row_svgs[0], out)
    else:
        compose.compose(row_svgs, out, natural=True)


def combined(sync=False):
    """-> (checked, stale, failed) for the whole figures."""
    import lib_compose as compose
    stale, failed, groups = [], [], combined_groups()
    with tempfile.TemporaryDirectory() as tmp:
        for whole, panels in groups.items():
            built = os.path.join(tmp, os.path.basename(whole))
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    build_whole(compose, panels, os.path.basename(whole), built,
                                tmp)
            except SystemExit as e:           # a panel whose letter is not found
                failed.append(f"{rel(whole)}: {e}")
                continue
            if os.path.exists(whole) and filecmp.cmp(built, whole, shallow=False):
                continue
            if sync:
                shutil.copyfile(built, whole)
            else:
                stale.append(rel(whole))
    return len(groups), stale, failed


PNG_ZOOM = "3"


def render_png(svg, png):
    subprocess.run(["rsvg-convert", "--zoom", PNG_ZOOM, "--background-color",
                    "white", svg, "-o", png], check=True)


def png_sources():
    """-> the SVGs that ship with a PNG: the whole figures and the Extended Data figures."""
    return sorted(set(combined_groups())
                  | set(glob.glob(os.path.join(FIGURES, "ed[0-9]*.svg"))))


def pngs(sync=False):
    """-> (checked, stale, skipped) for the PNGs."""
    wholes = png_sources()
    if not shutil.which("rsvg-convert"):
        return len(wholes), [], True
    stale = []
    with tempfile.TemporaryDirectory() as tmp:
        for svg in wholes:
            png = svg[:-4] + ".png"
            fresh = os.path.join(tmp, os.path.basename(png))
            render_png(svg, fresh)
            if os.path.exists(png) and filecmp.cmp(fresh, png, shallow=False):
                continue
            if sync:
                shutil.copyfile(fresh, png)
            else:
                stale.append(rel(png))
    return len(wholes), stale, False


SVG_NS = "{http://www.w3.org/2000/svg}"


def fonts():
    """-> ({font stack: files}, files with <text> that inherits no font)."""
    stacks, unstyled = {}, []
    for p in sorted(glob.glob(os.path.join(FIGURES, "*.svg"))):
        root = ET.parse(p).getroot()
        parent = {ch: el for el in root.iter() for ch in el}
        found = set()
        for t in root.iter(SVG_NS + "text"):
            el, ff = t, None
            while el is not None:
                ff = el.get("font-family")
                if ff:
                    break
                el = parent.get(el)
            if ff:
                found.add(" ".join(ff.split()))
            else:
                unstyled.append(rel(p))
                break
        for f in found:
            stacks.setdefault(f, []).append(rel(p))
    return stacks, unstyled


def _scale(transform):
    """-> the uniform scale factor a transform attribute applies."""
    k = 1.0
    for name, args in re.findall(r"(\w+)\(([^)]*)\)", transform or ""):
        v = [float(x) for x in re.split(r"[\s,]+", args.strip()) if x]
        if name == "scale":
            k *= v[0]
        elif name == "matrix":
            k *= abs(v[0] * v[3] - v[1] * v[2]) ** 0.5
    return k


def printed_type(svg):
    """-> (printed height in pt, [(size in pt, text)]) at the print width."""
    from lib_palette import PRINT_W_PT
    root = ET.parse(svg).getroot()
    _x, _y, w, h = (float(v) for v in root.get("viewBox").split())
    unit = PRINT_W_PT / w
    out = []

    def walk(el, k):
        if el.tag == SVG_NS + "text":
            txt = "".join(el.itertext()).strip()
            if txt and el.get("font-size"):
                out.append((float(el.get("font-size")) * k * unit, txt))
            return
        k *= _scale(el.get("transform"))
        for ch in el:
            walk(ch, k)
    walk(root, 1.0)
    return h * unit, out


def type_sizes():
    """-> [(file, printed height in pt, smallest size in pt, [texts below MIN_PT])]."""
    from lib_palette import MIN_PT
    rows = []
    for svg in png_sources():
        h, sizes = printed_type(svg)
        small = sorted((round(s, 1), t) for s, t in sizes if s < MIN_PT - 0.05)
        rows.append((rel(svg), h, min((s for s, _ in sizes), default=0.0), small))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sync", action="store_true",
                    help="rebuild the whole figures and PNGs")
    args = ap.parse_args(argv)

    bad = 0
    n, stale, failed = combined(args.sync)
    if stale:
        bad += len(stale)
        print(f"STALE         {len(stale)} whole figure(s) missing or behind their panels:")
        for s in stale:
            print(f"    {s}   (re-run with --sync)")
    elif args.sync:
        print(f"ok  {n} whole figure(s) rebuilt from their panels")
    else:
        print(f"ok  {n} whole figure(s) match their panels")
    if failed:
        bad += len(failed)
        print(f"UNBUILDABLE   {len(failed)} whole figure(s):")
        for f in failed:
            print(f"    {f}")

    n, stale, skipped = pngs(args.sync)
    if skipped:
        print(f"skip {n} whole-figure figure PNG(s): rsvg-convert is not installed")
    elif stale:
        bad += len(stale)
        print(f"STALE         {len(stale)} PNG(s) missing or behind their figure:")
        for s_ in stale:
            print(f"    {s_}   (re-run with --sync)")
    elif args.sync:
        print(f"ok  {n} figure PNG(s) rendered (whole figures and Extended Data)")
    else:
        print(f"ok  {n} figure PNG(s) match their SVG (whole figures and Extended Data)")

    from lib_palette import MIN_PT, PRINT_W_PT
    for f, h, smallest, small in type_sizes():
        size = (f"{PRINT_W_PT / 72:.1f} x {h / 72:.1f} in, smallest type "
                f"{smallest:.1f} pt")
        if small:
            bad += 1
            print(f"SMALL TYPE    {f}: {size}; {len(small)} text(s) under {MIN_PT:g} pt:")
            for s_, t in small[:6]:
                print(f"    {s_:>5.1f} pt  {t[:40]}")
        else:
            print(f"ok  {f}: {size}")

    stacks, unstyled = fonts()
    if unstyled:
        bad += len(unstyled)
        print(f"NO FONT       {len(unstyled)} figures/ file(s) with <text> that "
              f"inherits no font-family:")
        for u in sorted(set(unstyled)):
            print(f"    {u}")
    if len(stacks) > 1:
        bad += 1
        print(f"MIXED FONTS   {len(stacks)} different stacks across figures/:")
        for f, files in sorted(stacks.items(), key=lambda kv: -len(kv[1])):
            print(f"    {', '.join(sorted(set(files)))}\n        {f}")
    elif stacks and not unstyled:
        only = next(iter(stacks))
        print(f"ok  one font stack across figures/  ({only.split(',')[0]}, …)")

    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
