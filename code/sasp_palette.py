"""Validated chart palette (see the dataviz skill's reference instance).

Categorical slots are assigned in fixed order and never cycled; sequential
encoding is a single hue light->dark; signed bias uses the blue<->red diverging
pair with a NEUTRAL GRAY midpoint (never a hue at the midpoint).
"""
from matplotlib.colors import LinearSegmentedColormap

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

SEQ_STEPS = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("sasp_seq", SEQ_STEPS)

# diverging: blue pole <-> neutral gray <-> red pole
DIV = LinearSegmentedColormap.from_list(
    "sasp_div", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec",
                 "#f0a3a2", "#e34948", "#8f1f1e"])

STATUS = dict(good="#0ca30c", warning="#fab219",
              serious="#ec835a", critical="#d03b3b")


def apply_style(mpl):
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
        "axes.linewidth": 0.8, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.axisbelow": True,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.labelsize": 9, "axes.titlesize": 9.5,
        "font.size": 9, "legend.frameon": False, "legend.fontsize": 8,
        "lines.linewidth": 2.0, "lines.markersize": 5,
        "figure.dpi": 130, "savefig.dpi": 200,
        "text.color": INK,
    })
