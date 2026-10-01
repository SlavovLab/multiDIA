#!/usr/bin/env python3
"""One definition of the figure's colours and ink, shared by every panel."""

SLOTS = ["#eb6834", "#2a78d6", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

UNION = "#4a3aa7"
UNION_LABEL = "All"

# Digest-count ramp; index 0 is one digest.
DEPTH = ["#b9b1e4", "#7a6cc9", "#4a3aa7"]

# Peptide-depth ramps (1, 2, 3, 4+), light to dark, per categorical hue.
RAMPS = {
    "#1baf7a": ["#5fc7a2", "#32b787", "#189e6e", "#137a55"],
    "#2a78d6": ["#90b9ea", "#4c8edd", "#266cc1", "#1d5496"],
    "#eb6834": ["#f29e7d", "#ed7a4c", "#d45e2f", "#a44924"],
    "#4a3aa7": ["#aea6d7", "#6b5eb7", "#433496", "#342975"],
}
EMPTY = "#f0efec"
EMPTY_DARK = "#000000"

# The same ramps run dark to light, for a black background.
RAMPS_ON_DARK = {
    "#1baf7a": ["#189e6e", "#36b88a", "#67caa6", "#98dbc3"],
    "#2a78d6": ["#266cc1", "#4388db", "#71a5e4", "#9fc2ed"],
    "#eb6834": ["#d45e2f", "#ed7a4c", "#f29a78", "#f6bba4"],
    "#4a3aa7": ["#5c4eb0", "#776bbd", "#9289ca", "#aea6d7"],
}


DIVERGING_LOW = ["#86b6ef", "#5598e7", "#2a78d6", "#184f95"]
DIVERGING_HIGH = ["#fc8e87", "#fa5251", "#dd1e2d", "#970f1b"]
DIVERGING_MID = "#f0efec"
MISSING = "#ffffff"             # not measured; drawn hatched


def diverging(x, breaks=(0.5, 1.0, 2.0, 3.0)):
    """Binned colour for a signed value; None -> None (draw as missing)."""
    if x is None:
        return None
    m = abs(x)
    if m < breaks[0]:
        return DIVERGING_MID
    arm = DIVERGING_HIGH if x > 0 else DIVERGING_LOW
    for i, b in enumerate(breaks[1:]):
        if m < b:
            return arm[i]
    return arm[-1]


def ramp_for(colour, on_dark=False):
    """The depth ramp belonging to a categorical colour."""
    table = RAMPS_ON_DARK if on_dark else RAMPS
    return table.get(colour, table[UNION])


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
