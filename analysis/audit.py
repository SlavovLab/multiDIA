#!/usr/bin/env python3
"""Check figures/: whole figures match their panels, PNGs are current, one font stack.

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


# Whole figures laid out in columns rather than one stack
COLUMNS = {"fig2.svg": ["abd", "c"]}
COLUMN_GAP = 24.0


def build_whole(compose, panels, name, out, tmp):
    """Compose one whole figure from its panels, per `COLUMNS`."""
    cols = COLUMNS.get(name)
    if not cols:
        compose.compose(panels, out, natural=True, reletter=True)
        return
    by_letter = {os.path.basename(p)[-5]: p for p in panels}
    if sorted(by_letter) != sorted("".join(cols)):
        sys.exit(f"panels {sorted(by_letter)} do not match the columns {cols}")
    parts = []
    for i, letters in enumerate(cols):
        part = os.path.join(tmp, f"{name}.col{i}.svg")
        compose.compose([by_letter[x] for x in letters], part, natural=True,
                        reletter=True)
        parts.append(part)
    compose.compose_row(parts, out, gap=COLUMN_GAP, natural=True)


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
