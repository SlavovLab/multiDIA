#!/usr/bin/env python3
"""One definition of the figure's colours and ink, shared by every panel."""

SLOTS = ["#eb6834", "#2a78d6", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

UNION = "#4a3aa7"
UNION_LABEL = "All"

DIVERGING_HIGH = ["#fc8e87", "#fa5251", "#dd1e2d", "#970f1b"]
DIVERGING_MID = "#f0efec"


def _relative_luminance(hex_fill):
    """WCAG relative luminance: gamma-expanded, not raw sRGB."""
    def lin(c):
        c /= 255.0
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(int(hex_fill[i:i + 2], 16)) for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def ink_on(hex_fill):
    """White or ink for a label on a filled tile, whichever has the higher WCAG contrast."""
    lum = _relative_luminance(hex_fill)
    return "#ffffff" if 1.05 / (lum + 0.05) > (lum + 0.05) / 0.05 else INK

INK = "#0b0b0b"
INK_SECONDARY = "#1f1e1c"
INK_MUTED = "#2c2b28"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
STRIP_FILL = "#f0efec"
SURFACE = "#ffffff"

FONT = 'system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif'

# Font-size multiplier applied at write-out; `data-fit` text is exempt.
TEXT_BOOST = 1.3


def boost_type(svg):
    """-> the SVG with every font-size scaled by TEXT_BOOST, except fitted text."""
    import re

    def one(m):
        tag = m.group(0)
        if 'data-fit="1"' in tag:
            return tag
        return re.sub(r'font-size="([\d.]+)"',
                      lambda f: f'font-size="{float(f.group(1)) * TEXT_BOOST:.2f}"',
                      tag)
    return re.sub(r"<text\b[^>]*>", one, svg)

# Common drawing width, so composed panels share one type scale.
PANEL_W = 720.0


DISPLAY = {"gluc": "Glu-C", "lysc": "Lys-C", "trypsin": "Trypsin",
           "aspn": "Asp-N", "lysn": "Lys-N", "chymotrypsin": "Chymotrypsin",
           "argc": "Arg-C", "all": "All proteases"}


def display(name):
    """The reader-facing spelling of a protease or the union label."""
    return DISPLAY.get(canon(name), str(name))


def canon(name):
    """Fold a protease name so 'Glu-C', 'GluC' and 'gluc' are one category."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def assign(categories):
    """Map category names to colours, stably and independently of the data."""
    singles = [c for c in categories if canon(c) != canon(UNION_LABEL)]
    slot = {k: SLOTS[i % len(SLOTS)]
            for i, k in enumerate(sorted({canon(c) for c in singles}))}
    out = {c: slot[canon(c)] for c in singles}
    for c in categories:
        if canon(c) == canon(UNION_LABEL):
            out[c] = UNION
    return out
