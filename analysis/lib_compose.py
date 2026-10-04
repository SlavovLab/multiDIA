#!/usr/bin/env python3
"""Stack SVG panels into one figure.

    python3 lib_compose.py strip.svg --letter a --key Glu-C,Lys-C,Trypsin --out figures/fig2a.svg
"""

import argparse
import os
import re
import sys


LETTER = re.compile(
    r'<text x="2[0-9](?:\.\d+)?" y="3[0-9](?:\.\d+)?"[^>]*font-weight="600"'
    r'>[a-z]</text>')


def drop_letter(body):
    """Remove a panel's own letter."""
    return LETTER.sub("", body, count=1)


BOLD_CHAR = re.compile(r"<text\b([^>]*)>([a-z])</text>")


def take_letter(body, path):
    """-> (body without its letter, the letter). Exits unless there is exactly one."""
    hits = []
    for m in BOLD_CHAR.finditer(body):
        a = m.group(1)
        x = re.search(r"""\bx=["'](-?[\d.]+)["']""", a)
        y = re.search(r"""\by=["'](-?[\d.]+)["']""", a)
        if (re.search(r"""font-weight=["']600["']""", a) and x and y
                and float(x.group(1)) < 50 and float(y.group(1)) < 50):
            hits.append(m)
    if len(hits) != 1:
        sys.exit(f"{path}: expected one panel letter in the top-left corner, "
                 f"found {len(hits)}")
    m = hits[0]
    return body[:m.start()] + body[m.end():], m.group(2)


LETTER_PT = 13.0


def stamp_letter(letter, x, y, font):
    """The one letter a composite draws, in output units, at its panel's corner."""
    from lib_palette import pt
    return (f'<text x="{x:g}" y="{y:g}" font-size="{pt(LETTER_PT):.2f}" fill="#0b0b0b" '
            f"text-anchor='start' font-weight='600' "
            f"font-family='{font}'>{letter}</text>")


def key_row(names, right, y, font, size_pt=11.0):
    """-> SVG for a colour key, right-aligned to end at `right`."""
    from lib_palette import INK_SECONDARY, assign, pt
    size = pt(size_pt)
    pairs = [n.split("=", 1) if "=" in n else (n, n) for n in names]
    names = [lab for lab, _ in pairs]
    by_cat = assign([cat for _, cat in pairs])
    colour = {lab: by_cat[cat] for lab, cat in pairs}
    sw, gap, pad = 0.75 * size, 0.35 * size, 1.2 * size
    widths = [sw + gap + size * 0.56 * len(n) for n in names]
    total = sum(widths) + pad * (len(names) - 1)
    if total > right - 8:
        k = (right - 8) / total
        size, sw, gap, pad = size * k, sw * k, gap * k, pad * k
        widths = [w * k for w in widths]
    x = right - sum(widths) - pad * (len(names) - 1)
    parts = []
    for n, w in zip(names, widths):
        parts.append(f'<rect x="{x:.1f}" y="{y - sw:.1f}" width="{sw:.1f}" '
                     f'height="{sw:.1f}" rx="2" fill="{colour[n]}"/>')
        parts.append(f'<text x="{x + sw + gap:.1f}" y="{y:.1f}" '
                     f'font-size="{size:.2f}" fill="{INK_SECONDARY}" '
                     f"text-anchor='start' font-family='{font}'>{n}</text>")
        x += w + pad
    return parts


def parse(path):
    """-> (body, viewBox_w, viewBox_h, out_w, font_family). Strips the <svg>."""
    src = open(path).read()
    m = re.search(r"<svg[^>]*>", src)
    if not m:
        sys.exit(f"{path}: no <svg> element")
    head = m.group(0)
    vb = re.search(r'viewBox="([^"]+)"', head)
    if not vb:
        sys.exit(f"{path}: no viewBox")
    _x, _y, w, h = (float(v) for v in vb.group(1).split())
    ow = re.search(r'width="([\d.]+)"', head)
    ff = re.search(r"""font-family=(["'])(.*?)\1""", head, re.S)
    body = src[m.end():]
    body = body[:body.rfind("</svg>")]
    font = ff.group(2) if ff else None
    if font is None:
        inner = {s for _q, s in re.findall(r"""<g [^>]*font-family=(["'])(.*?)\1""",
                                           body, re.S)}
        if len(inner) == 1:
            font = inner.pop()
    return body, w, h, float(ow.group(1)) if ow else w, font


def compose_row(paths, out, gap=0.0, natural=False):
    """Place panels side by side at a common height (`natural`: own scale). -> (w, h)."""
    panels = [parse(p) for p in paths]
    if natural:
        scaled = [(body, ow / w, ow, ff) for body, w, h, ow, ff in panels]
        target = max(h * ow / w for _b, w, h, ow, _f in panels)
    else:
        target = max(p[2] for p in panels)
        scaled = [(body, target / h, w * target / h, ff)
                  for body, w, h, _ow, ff in panels]
    W = sum(s[2] for s in scaled) + gap * (len(scaled) - 1)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" '
             f'height="{target:.0f}" viewBox="0 0 {W:.0f} {target:.0f}">',
             f'<rect width="{W:.0f}" height="{target:.0f}" fill="#ffffff"/>']
    x = 0.0
    for (body, k, w, ff) in scaled:
        parts.append(f'<g transform="translate({x:.2f} 0) scale({k:.6f})"'
                     + (f" font-family='{ff}'" if ff else "") + ">")
        parts.append(body)
        parts.append("</g>")
        x += w + gap
    parts.append("</svg>")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(parts) + "\n")
    print(f"  wrote {out}  ({W:.0f}x{target:.0f})")
    for p, (_b, w, h, _o, _f), (_b2, k, sw_, _f2) in zip(paths, panels, scaled):
        print(f"    {os.path.basename(p):<28s} {w:>6.0f} -> {sw_:>6.0f}  "
              f"(x{k:.3f})")
    return W, target


def compose(paths, out, gap=0.0, width=None, letter=None, key=None, natural=False,
            reletter=False, pad_top=0.0):
    """Stack panels. -> (width, height)."""
    panels = [parse(p) for p in paths]
    own = [None] * len(panels)
    if reletter:
        taken = [take_letter(b, p) for (b, *_rest), p in zip(panels, paths)]
        own = [ch for _b, ch in taken]
        panels = [(b2, w, h, ow, ff)
                  for (b2, _ch), (_b, w, h, ow, ff) in zip(taken, panels)]
    if letter is not None:
        panels = [(drop_letter(b), w, h, ow, ff) for b, w, h, ow, ff in panels]
    if natural:
        target = max(p[3] for p in panels)
        scaled = [(body, ow / w, h * ow / w, ff) for body, w, h, ow, ff in panels]
    else:
        target = width or max(p[3] for p in panels)
        scaled = [(body, target / w, h * target / w, ff)
                  for body, w, h, _ow, ff in panels]
    H = pad_top + sum(s[2] for s in scaled) + gap * (len(scaled) - 1)

    fonts = {s[3] for s in scaled if s[3]}
    if len(fonts) > 1:
        print(f"  note: {len(fonts)} different font stacks across the panels")

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{target:.0f}" '
             f'height="{H:.0f}" viewBox="0 0 {target:.0f} {H:.0f}">',
             f'<rect width="{target:.0f}" height="{H:.0f}" fill="#ffffff"/>']
    y, tops = pad_top, []
    for (body, k, h, ff) in scaled:
        tops.append(y)
        # single-quoted: the font stack itself contains double quotes
        parts.append(f'<g transform="translate(0 {y:.2f}) scale({k:.6f})"'
                     + (f" font-family='{ff}'" if ff else "") + ">")
        parts.append(body)
        parts.append("</g>")
        y += h + gap
    # after the panels, which paint their own white background
    if letter:
        parts.append(stamp_letter(letter, 22, 32, next(iter(fonts), "")))
    first = next((s[3] for s in scaled if s[3]), "")
    for ch, top, s in zip(own, tops, scaled):
        if ch:
            parts.append(stamp_letter(ch, 22, top + 32, s[3] or first))
    if key:
        parts.extend(key_row(key, target - 22, 30, next(iter(fonts), "")))
    parts.append("</svg>")

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(parts) + "\n")
    print(f"  wrote {out}  ({target:.0f}x{H:.0f})")
    for p, (_b, _w, h, _o, _f), (_b2, k, sh, _f2) in zip(paths, panels, scaled):
        print(f"    {os.path.basename(p):<28s} {h:>6.0f} -> {sh:>6.0f}  "
              f"(x{k:.3f})")
    return target, H


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("svgs", nargs="+", help="panels, top to bottom")
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=float, default=None,
                    help="output width; default is the widest input")
    ap.add_argument("--letter", default=None,
                    help="stamp this letter and drop the panels' own")
    ap.add_argument("--row", action="store_true",
                    help="place panels side by side at a common height")
    ap.add_argument("--key", default=None, metavar="A,B,C",
                    help="colour key categories; Label=Category allowed")
    ap.add_argument("--pad-top", type=float, default=0.0,
                    help="blank strip above the panels, in output units")
    args = ap.parse_args(argv)
    if args.row:
        if args.width or args.letter or args.key or args.pad_top:
            ap.error("--row takes no --width, --letter, --key or --pad-top")
        compose_row(args.svgs, args.out)
    else:
        compose(args.svgs, args.out, width=args.width, letter=args.letter,
                key=args.key.split(",") if args.key else None, pad_top=args.pad_top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
