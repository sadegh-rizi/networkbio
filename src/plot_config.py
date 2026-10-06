"""Shared plotting style and figure saving for the networkbio project.

Every figure goes through `save_figure`, which writes three files: a vector
PDF plus a presentation-resolution and a manuscript-resolution PNG. The DPI
constants are deliberate placeholders; tune them once a venue is known (see
`doc/agent-rules.md`, "Figures and tables").

Colors are defined once here so that two panels do not need two legends.
Add categories as the analysis defines them; do not hard-code a palette in
an individual figure.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Placeholders, not specifications. Adjust per venue.
PRESENTATION_DPI = 150
MANUSCRIPT_DPI = 600

# Okabe-Ito, colorblind-safe. Base for any new category mapping.
OKABE_ITO = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#000000",  # black
]

GROUP_COLORS = {
    "Ctrl": "#0072B2",
    "PMS": "#D55E00",
}

TREATMENT_COLORS = {
    "vehicle": "#56B4E9",
    "simvastatin": "#E69F00",
}

# Network node roles (inputs, intermediates, measured nodes).
NODE_ROLE_COLORS = {
    "input": "#CC79A7",
    "intermediate": "#999999",
    "measured_up": "#D55E00",
    "measured_down": "#0072B2",
}


def apply_plot_style() -> None:
    """Set project-wide Matplotlib defaults. Call once per notebook."""
    mpl.rcParams.update(
        {
            "figure.dpi": 100,
            "savefig.bbox": "tight",
            "savefig.transparent": False,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(color=OKABE_ITO),
            "legend.frameon": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,  # embed as TrueType so text stays editable
            "ps.fonttype": 42,
        }
    )


def validate_categories(observed: list[str], expected: dict[str, str], name: str) -> None:
    """Fail loudly if a column contains a category with no assigned color.

    Silently falling back to a default color is how two figures end up using
    different colors for the same group.
    """
    missing = sorted(set(observed) - set(expected))
    if missing:
        raise ValueError(
            f"{name}: no color defined for {missing}. "
            f"Add them to the mapping in src/plot_config.py rather than "
            f"overriding the palette in this figure."
        )


def save_figure(fig: plt.Figure, path: Path) -> list[Path]:
    """Write a figure as vector PDF plus presentation and manuscript PNGs.

    `path` is the stem, without extension. Produces `<stem>.pdf`,
    `<stem>_presentation.png`, and `<stem>_manuscript.png` in the same
    folder, per `doc/agent-rules.md`. Returns the paths written.
    """
    path = Path(path)
    if path.suffix:
        raise ValueError(f"pass a stem without an extension, got {path}")
    path.parent.mkdir(parents=True, exist_ok=True)

    written = [
        path.with_suffix(".pdf"),
        path.with_name(f"{path.name}_presentation.png"),
        path.with_name(f"{path.name}_manuscript.png"),
    ]
    fig.savefig(written[0])
    fig.savefig(written[1], dpi=PRESENTATION_DPI)
    fig.savefig(written[2], dpi=MANUSCRIPT_DPI)
    return written
