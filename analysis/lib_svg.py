"""The SVG canvas every figure draws on."""

from lib_palette import AXIS, INK, INK_MUTED, SURFACE, TEXT_BOOST, boost_type


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Canvas:
    def __init__(self, w, h, font, font_scale=1.0, out_w=None):
        self.w, self.h = w, h
        self.fs = font_scale
        ow = out_w or w
        oh = h * ow / w
        self.parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{ow:.0f}" '
            f'height="{oh:.0f}" viewBox="0 0 {w} {h}" font-family=\'{font}\'>',
            f'<rect width="{w}" height="{h}" fill="{SURFACE}"/>',
        ]

    def add(self, s):
        self.parts.append(s)

    def rect(self, x, y, w, h, fill="none", stroke="none", sw=1.0, rx=0,
             fo=None, so=None):
        a = (f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0.4):.1f}" '
             f'height="{max(h, 0.4):.1f}" fill="{fill}" stroke="{stroke}" '
             f'stroke-width="{sw}"')
        if rx:
            a += f' rx="{rx}"'
        if fo is not None:
            a += f' fill-opacity="{fo}"'
        if so is not None:
            a += f' stroke-opacity="{so}"'
        self.add(a + "/>")

    def line(self, x1, y1, x2, y2, stroke=INK_MUTED, sw=1.0, so=None):
        a = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
             f'stroke="{stroke}" stroke-width="{sw}"')
        if so is not None:
            a += f' stroke-opacity="{so}"'
        self.add(a + "/>")

    def text(self, x, y, s, size=10, fill=INK, anchor="start", weight=None,
             italic=False, spacing=None, rot=None, fit=False):
        """`fit` marks text that must stay inside its box, so it is not boosted."""
        a = (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size * self.fs:.2f}" '
             f'fill="{fill}" text-anchor="{anchor}"' + (' data-fit="1"' if fit else ''))
        if weight:
            a += f' font-weight="{weight}"'
        if italic:
            a += ' font-style="italic"'
        if spacing:
            a += f' letter-spacing="{spacing}"'
        if rot:
            # rotate about the anchor
            a += f' transform="rotate({rot:.1f} {x:.1f} {y:.1f})"'
        self.add(a + f">{esc(s)}</text>")

    def arrow(self, x0, y0, x1, y1, xmid=None, stroke=INK_MUTED, sw=1.4, head=5.0):
        """Orthogonal connector with rounded corners and a filled head at (x1,y1)."""
        tip = x1
        x1 -= head * 1.25
        if abs(y1 - y0) < 0.5:
            self.add(f'<path d="M {x0:.1f} {y0:.1f} L {x1:.1f} {y1:.1f}" '
                     f'fill="none" stroke="{stroke}" stroke-width="{sw}"/>')
        else:
            xm = xmid if xmid is not None else (x0 + x1) / 2
            r = min(7.0, abs(y1 - y0) / 2, max(abs(xm - x0), 1), max(abs(x1 - xm), 1))
            s = 1 if y1 > y0 else -1
            self.add(
                f'<path d="M {x0:.1f} {y0:.1f} L {xm - r:.1f} {y0:.1f} '
                f'Q {xm:.1f} {y0:.1f} {xm:.1f} {y0 + s * r:.1f} '
                f'L {xm:.1f} {y1 - s * r:.1f} '
                f'Q {xm:.1f} {y1:.1f} {xm + r:.1f} {y1:.1f} '
                f'L {x1:.1f} {y1:.1f}" fill="none" stroke="{stroke}" '
                f'stroke-width="{sw}" stroke-linecap="round"/>')
        self.add(f'<path d="M {tip:.1f} {y1:.1f} L {tip - head * 1.3:.1f} '
                 f'{y1 - head * 0.62:.1f} L {tip - head * 1.3:.1f} '
                 f'{y1 + head * 0.62:.1f} Z" fill="{stroke}"/>')

    def chip(self, x, y, w, h, title, sub=None, colour=None, bold=True):
        """A rounded box; tinted in `colour` if given, hairline neutral if not."""
        if colour:
            self.rect(x, y, w, h, fill=colour, fo=0.10, stroke=colour, sw=1.6, rx=7)
        else:
            self.rect(x, y, w, h, fill=SURFACE, stroke=AXIS, sw=1.1, rx=7)
        cx = x + w / 2
        if sub:
            # stack the two lines by their boosted printed heights
            t, u = 11 * self.fs * TEXT_BOOST, 8.8 * self.fs * TEXT_BOOST
            top = y + (h - (t + u * 1.05)) / 2
            self.text(cx, top + 0.8 * t, title, 11, INK, "middle",
                      "600" if bold else None)
            self.text(cx, top + t + 0.85 * u, sub, 8.8, INK_MUTED, "middle")
        else:
            self.text(cx, y + h / 2 + 4, title, 11, INK, "middle",
                      "600" if bold else None)

    def out(self):
        return boost_type("\n".join(self.parts + ["</svg>"]) + "\n")
