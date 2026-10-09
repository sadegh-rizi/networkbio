"""Draw a CORNETO/CARNIVAL solution on top of its signed prior-knowledge network.

Shared by the toy figure and, later, the stage-04 real-data figures, so the
visual encoding is defined once:

- Grey dashed edges: PKN edges the solver did not select.
- Dark edges: selected edges. Arrow head = activation (PKN sign +1), flat bar
  = inhibition (PKN sign -1). With several samples, line width and the label
  give the fraction of samples that use the edge.
- Node fill: inferred node state (+1 up, -1 down, white = 0 / not used).
  Colours come from `plot_config.NODE_ROLE_COLORS`.
- Node outline: input nodes (pink, thick), measured nodes (black, thick).
  Measured nodes also carry the observed sign underneath, so the fit can be
  read off the figure without colour alone.
- Node shape: circle = protein or gene, diamond = metabolite.

Edge-value convention (CORNETO CarnivalFlow): the value of a selected edge is
the sign of the signal it transmits, i.e. state(source) * PKN sign, and it
equals state(target). `node_states` checks this and raises if it does not
hold, so a misread solution table fails loudly instead of being drawn.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch

from plot_config import NODE_ROLE_COLORS

Edge = tuple[str, str]
METABOLITE_PREFIX = "Metab__"
UNUSED_FILL = "#FFFFFF"
UNSELECTED_EDGE = "#BBBBBB"
SELECTED_EDGE = "#333333"
NODE_SIZE = 950  # scatter marker area in pt^2
EDGE_SHRINK = 17  # points; roughly the marker radius


def pkn_signs(pkn: Iterable[tuple[str, int, str]]) -> dict[Edge, int]:
    """Map (source, target) -> sign from SIF-style (source, sign, target) tuples."""
    signs: dict[Edge, int] = {}
    for src, sign, tgt in pkn:
        if (src, tgt) in signs and signs[(src, tgt)] != sign:
            raise ValueError(f"edge {src}->{tgt} has both signs in the PKN")
        signs[(src, tgt)] = int(sign)
    return signs


def node_states(selected: pd.DataFrame, signs: Mapping[Edge, int]) -> dict[str, int]:
    """Infer node states for one sample from its selected edges.

    `selected` has columns source, target, edge_value. A target's state is the
    value of its incoming selected edge; a node with no incoming selected
    edge (a source of the solution) gets edge_value * PKN sign of an outgoing
    edge. Raises if edges disagree or an edge is not in the PKN.
    """
    states: dict[str, int] = {}

    def _set(node: str, value: int) -> None:
        if states.setdefault(node, value) != value:
            raise ValueError(f"inconsistent state for {node}: {states[node]} vs {value}")

    for row in selected.itertuples(index=False):
        edge = (row.source, row.target)
        if edge not in signs:
            raise ValueError(f"selected edge {edge} is not in the PKN")
        value = int(row.edge_value)
        _set(row.target, value)
        _set(row.source, value * signs[edge])
    return states


def layered_positions(nodes: Iterable[str], edges: Iterable[Edge]) -> dict[str, tuple[float, float]]:
    """Left-to-right layered layout for a directed graph.

    Layer = longest path from a node without incoming edges (BFS relaxation,
    capped at the node count so cycles cannot loop forever; in a cyclic graph
    the layer is then only approximate). Within a layer, nodes are ordered by
    the mean y of their predecessors, which removes most crossings in small
    graphs. Returns node -> (x, y).
    """
    nodes = list(dict.fromkeys(nodes))
    edges = list(edges)
    preds: dict[str, list[str]] = defaultdict(list)
    for s, t in edges:
        preds[t].append(s)
    layer = {n: 0 for n in nodes}
    for _ in range(len(nodes)):
        changed = False
        for s, t in edges:
            if layer[t] < layer[s] + 1 and layer[s] + 1 < len(nodes):
                layer[t] = layer[s] + 1
                changed = True
        if not changed:
            break

    pos: dict[str, tuple[float, float]] = {}
    for depth in sorted(set(layer.values())):
        members = [n for n in nodes if layer[n] == depth]

        def _key(n: str) -> tuple[float, int]:
            ys = [pos[p][1] for p in preds[n] if p in pos]
            return (-(sum(ys) / len(ys)) if ys else 0.0, nodes.index(n))

        members.sort(key=_key)
        k = len(members)
        for i, n in enumerate(members):
            pos[n] = (float(depth), (k - 1) / 2 - i)
    return pos


def edge_table(pkn: Iterable[tuple[str, int, str]], solutions: pd.DataFrame) -> pd.DataFrame:
    """One row per PKN edge with how often each sample's solution selected it.

    `solutions` has columns sample, source, target, edge_value. Output columns:
    source, target, pkn_sign, n_samples_selected, n_samples, frequency,
    samples (comma-separated). Suitable as a Cytoscape edge table.
    """
    signs = pkn_signs(pkn)
    samples = list(dict.fromkeys(solutions["sample"]))
    used = solutions.groupby(["source", "target"])["sample"].apply(lambda s: sorted(set(s)))
    rows = []
    for (s, t), sign in signs.items():
        hit = used.get((s, t), [])
        rows.append(
            {
                "source": s,
                "target": t,
                "pkn_sign": sign,
                "n_samples_selected": len(hit),
                "n_samples": len(samples),
                "frequency": len(hit) / len(samples) if samples else 0.0,
                "samples": ",".join(hit),
            }
        )
    return pd.DataFrame(rows)


def draw_signed_network(
    ax: plt.Axes,
    pkn: list[tuple[str, int, str]],
    solutions: pd.DataFrame,
    *,
    inputs: Iterable[str],
    measured: Mapping[str, int],
    labels: Mapping[str, str] | None = None,
    pos: Mapping[str, tuple[float, float]] | None = None,
    title: str = "",
) -> dict[str, tuple[float, float]]:
    """Draw the PKN with one or more samples' selected edges overlaid.

    With several samples, edges are weighted by selection frequency and node
    fill is the state shared by all samples that use the node (hatched if the
    samples disagree). Returns the positions used, so panels can share them.
    """
    signs = pkn_signs(pkn)
    nodes = list(dict.fromkeys([n for s, _, t in pkn for n in (s, t)]))
    pos = dict(pos) if pos is not None else layered_positions(nodes, signs)
    labels = dict(labels or {})
    inputs = set(inputs)

    per_sample = {
        name: node_states(df[["source", "target", "edge_value"]], signs)
        for name, df in solutions.groupby("sample", sort=False)
    }
    etab = edge_table(pkn, solutions).set_index(["source", "target"])
    n_samples = max(int(etab["n_samples"].max()), 1)

    for (s, t), sign in signs.items():
        freq = float(etab.loc[(s, t), "frequency"])
        selected = freq > 0
        arrow = FancyArrowPatch(
            pos[s],
            pos[t],
            arrowstyle="-|>" if sign > 0 else "-[, widthB=0.45, lengthB=0",
            mutation_scale=12 if selected else 9,
            shrinkA=EDGE_SHRINK,
            shrinkB=EDGE_SHRINK + (0 if sign > 0 else 4),  # keep the bar clear of the node ring
            color=SELECTED_EDGE if selected else UNSELECTED_EDGE,
            linewidth=(0.8 + 2.2 * freq) if selected else 0.8,
            linestyle="-" if selected else (0, (3, 2)),
            zorder=2 if selected else 1,
        )
        ax.add_patch(arrow)
        if selected and n_samples > 1:
            (x0, y0), (x1, y1) = pos[s], pos[t]
            ax.text(
                (x0 + x1) / 2,
                (y0 + y1) / 2 + 0.12,
                f"{int(etab.loc[(s, t), 'n_samples_selected'])}/{n_samples}",
                fontsize=6.5,
                ha="center",
                va="bottom",
                color=SELECTED_EDGE,
                zorder=4,
            )

    fill_for = {1: NODE_ROLE_COLORS["measured_up"], -1: NODE_ROLE_COLORS["measured_down"]}
    for n in nodes:
        values = {st[n] for st in per_sample.values() if n in st}
        mixed = len(values) > 1
        state = values.pop() if len(values) == 1 else 0
        fill = fill_for.get(state, UNUSED_FILL) if not mixed else UNUSED_FILL
        if n in inputs:
            edge_c, lw = NODE_ROLE_COLORS["input"], 2.6
        elif n in measured:
            edge_c, lw = "#000000", 2.0
        else:
            edge_c, lw = NODE_ROLE_COLORS["intermediate"], 0.9
        ax.scatter(
            *pos[n],
            s=NODE_SIZE,
            marker="D" if n.startswith(METABOLITE_PREFIX) else "o",
            facecolor=fill,
            edgecolor=edge_c,
            linewidth=lw,
            hatch="////" if mixed else None,
            zorder=3,
        )
        used = any(n in st for st in per_sample.values())
        ax.text(
            *pos[n],
            labels.get(n, n),
            ha="center",
            va="center",
            fontsize=7 if len(labels.get(n, n)) > 4 else 8,
            color="#FFFFFF" if fill != UNUSED_FILL else ("#000000" if used else "#888888"),
            zorder=4,
        )
        if n in measured:
            obs = "obs \u2191" if measured[n] > 0 else ("obs \u2193" if measured[n] < 0 else "obs 0")
            dy = 0.62 if n.startswith(METABOLITE_PREFIX) else 0.40  # the diamond is taller
            ax.text(pos[n][0], pos[n][1] - dy, obs, ha="center", va="top", fontsize=6.5, zorder=4)

    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.5, max(xs) + 0.5)
    ax.set_ylim(min(ys) - 0.9, max(ys) + 0.55)
    ax.set_title(title, loc="left")
    ax.set_axis_off()
    return pos


def network_legend_handles() -> list:
    """Legend entries for the encoding used by `draw_signed_network`."""
    marker = dict(marker="o", linestyle="none", markersize=9)
    return [
        Patch(facecolor=NODE_ROLE_COLORS["measured_up"], label="state +1 (up)"),
        Patch(facecolor=NODE_ROLE_COLORS["measured_down"], label="state \u22121 (down)"),
        Patch(facecolor=UNUSED_FILL, edgecolor="#999999", label="state 0 / not used"),
        Line2D([], [], markerfacecolor="white", markeredgecolor=NODE_ROLE_COLORS["input"],
               markeredgewidth=2.6, label="input node", **marker),
        Line2D([], [], markerfacecolor="white", markeredgecolor="#000000",
               markeredgewidth=2.0, label="measured node (obs = observed sign)", **marker),
        Line2D([], [], marker="D", linestyle="none", markersize=8, markerfacecolor="white",
               markeredgecolor="#999999", label="metabolite"),
        Line2D([], [], color=SELECTED_EDGE, linewidth=2, label="selected edge (\u2192 activation, \u22a3 inhibition)"),
        Line2D([], [], color=UNSELECTED_EDGE, linewidth=0.8, linestyle=(0, (3, 2)), label="PKN edge, not selected"),
    ]
