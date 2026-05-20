#!/usr/bin/env python3
"""
plot_cell_migration_composite.py
---------------------------------
Publication-quality composite figure: a single cell's migration trajectory
overlaid on the steady-state glucose concentration gradient.

For each selected time frame the cell is drawn as a filled circle coloured
by either its motility speed fraction ("velocity") or its growth rate.
The glucose field from the final simulation frame forms the background.

Because the domain is very wide (700 µm) relative to the cell diameter
(~2 µm), cells are drawn at an inflated display radius (--cell-radius, default
8 µm) so they are visible on the printed figure; this should be noted in the
figure caption.  The y-axis is kept at its true extent (±10 µm) and the axes
aspect ratio is set to "auto" so the full spatial context is preserved.

Usage
-----
    python plot_cell_migration_composite.py [options]

    --output-dir  DIR    PhysiCell output directory          [../output]
    --results-dir DIR    Where to write the figure           [../results/composite]
    --color-by    STR    "velocity" | "growth_rate"          [velocity]
    --n-frames    N      Number of time snapshots to overlay [12]
    --dpi         N      Output resolution                   [300]
    --cell-radius R      Display radius (µm, None=auto)      [auto]
    --x-scale     STR    "linear" | "symlog"                 [linear]
    --glc-frame   STR    "last" | "first"                    [last]
    --no-glc-bg          Disable glucose background (cells only)
    --glc-cmap STR       Colormap for glucose background      [Blues]
    --cell-cmap STR      Colormap for cell metric             [auto]
    --cell-edgecolor STR Cell outline color                   [dark green]
    --cell-linewidth R   Cell outline width                   [0.6]
    --glc-alpha R        Glucose background opacity           [0.85]

Recommended journal settings
----------------------------
    For a motility figure with an alive/proliferative visual language:
        --color-by velocity --glc-cmap Blues --cell-cmap summer

    For a growth-focused figure:
        --color-by growth_rate --glc-cmap Blues --cell-cmap Greens
"""

from __future__ import annotations

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")  # non-interactive backend – safe on HPC clusters
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.interpolate import griddata
from pctk import multicellds

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT   = _SCRIPT_DIR.parent
DEFAULT_OUTPUT_DIR  = str(_REPO_ROOT / "output")
DEFAULT_RESULTS_DIR = str(_REPO_ROOT / "results" / "composite")

# ---------------------------------------------------------------------------
# Microenvironment column indices
# (0-2 = x/y/z; 3 = voxel volume; 4 = oxygen; 5 = glucose; …)
# ---------------------------------------------------------------------------
MICRO_IDX_GLUCOSE = 5

# ---------------------------------------------------------------------------
# Colour-metric mapping
# ---------------------------------------------------------------------------
_METRIC_COL = {
    "velocity":    "migration_speed",
    "growth_rate": "growth_rate",
}
_METRIC_LABEL = {
    "velocity":    "Migration speed (µm/min)",
    "growth_rate": "Growth rate (h⁻¹)",
}
_METRIC_CMAP = {
    # Green/yellow palette avoids the “black = dead” impression while staying
    # visible on a blue glucose background.
    "velocity":    "summer",
    "growth_rate": "Greens",
}


# ===========================================================================
# Helper utilities (reused from dynamics_analysis.py patterns)
# ===========================================================================

def _get_cell_columns(output_folder: str) -> list[str]:
    """Parse cell-matrix column labels from initial.xml."""
    xml_fname = os.path.join(output_folder, "initial.xml")
    tree = ET.parse(xml_fname)
    root = tree.getroot()

    def _find_simplified(node):
        for child in node:
            if child.tag == "simplified_data" and child.attrib.get("source") == "PhysiCell":
                return child
        for child in node:
            r = _find_simplified(child)
            if r is not None:
                return r

    node   = _find_simplified(root.find("cellular_information"))
    labels = node.find("labels")

    cols: list[str] = []
    for child in labels:
        name = child.text.strip()
        size = int(child.attrib["size"])
        if size == 1:
            cols.append(name)
        elif size == 3:
            cols.extend([f"x_{name}", f"y_{name}", f"z_{name}"])
        elif size == 2:
            cols.extend([f"x_{name}", f"y_{name}"])
        else:
            cols.extend([f"{name}_{i}" for i in range(size)])
    return cols


def _ensure_rows(a: np.ndarray, ncols: int) -> np.ndarray:
    """Ensure the matrix is (n_cells, n_features)."""
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError(f"Expected 2-D matrix, got {a.shape}")
    if a.shape[1] == ncols:
        return a
    if a.shape[0] == ncols:
        return a.T
    raise ValueError(f"Cannot orient matrix: shape={a.shape}, expected ncols={ncols}")


def _get_domain_bounds(reader: multicellds.MultiCellDS) -> tuple[tuple, tuple]:
    """Read domain extents from initial.xml bounding_box."""
    root    = reader._tree.getroot()
    bb_text = root.find(".//microenvironment/domain/mesh/bounding_box").text
    xmin, ymin, _zmin, xmax, ymax, _zmax = map(float, bb_text.split())
    return (xmin, xmax), (ymin, ymax)


# ===========================================================================
# Publication style
# ===========================================================================

def _trimmed_cmap(name: str, vmin: float = 0.15, vmax: float = 0.95):
    """Return a colormap with extreme dark/white ends removed.

    This is useful for cell overlays: it avoids very dark cells that read as
    dead/necrotic and avoids nearly-white cells that disappear on pale regions.
    """
    base = mpl.colormaps[name]
    vmin = max(0.0, min(1.0, float(vmin)))
    vmax = max(0.0, min(1.0, float(vmax)))
    if vmax <= vmin:
        vmin, vmax = 0.0, 1.0
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"{name}_trimmed_{vmin:.2f}_{vmax:.2f}",
        base(np.linspace(vmin, vmax, 256)),
    )


def _apply_publication_style() -> None:
    """Apply Nature-style matplotlib rcParams (7 pt, sans-serif)."""
    mpl.rcParams.update({
        # Font
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":         7,
        "axes.titlesize":    8,
        "axes.labelsize":    7,
        "xtick.labelsize":   6,
        "ytick.labelsize":   6,
        "legend.fontsize":   6,
        # Lines
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size":  3,
        "ytick.major.size":  3,
        "lines.linewidth":   1.0,
        # Layout
        "figure.dpi":        300,
        "savefig.dpi":       300,
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.02,
        # Colorbars
        "image.interpolation": "bilinear",
    })


# ===========================================================================
# Data loading
# ===========================================================================

def _load_glucose_background(
    reader: multicellds.MultiCellDS,
    frame: str,
    domain_x: tuple[float, float],
    domain_y: tuple[float, float],
    grid_nx: int = 700,
    grid_ny: int = 100,
) -> tuple[np.ndarray, float, float]:
    """
    Interpolate the glucose field onto a regular grid.

    Returns
    -------
    Gi        : (grid_ny, grid_nx) float array
    glc_vmin  : float
    glc_vmax  : float
    """
    all_frames = list(reader.microenvironment_as_matrix_iterator())
    if not all_frames:
        raise RuntimeError("No microenvironment data found.")

    if frame == "first":
        _t, m = all_frames[0]
    else:  # "last"
        _t, m = all_frames[-1]

    vox_x = m[0, :]
    vox_y = m[1, :]
    glc   = m[MICRO_IDX_GLUCOSE, :]

    xi     = np.linspace(domain_x[0], domain_x[1], grid_nx)
    yi     = np.linspace(domain_y[0], domain_y[1], grid_ny)
    Xi, Yi = np.meshgrid(xi, yi)
    Gi     = griddata((vox_x, vox_y), glc, (Xi, Yi), method="linear")

    # Use the colour scale from the displayed frame only so the gradient
    # contrast is maximal (using all frames would average out most variation).
    glc_vmin = float(np.nanmin(glc))
    glc_vmax = float(np.nanmax(glc))
    # Protect against degenerate (uniform) fields
    if glc_vmax <= glc_vmin:
        glc_vmax = glc_vmin + 1e-9

    return Gi, glc_vmin, glc_vmax


def _load_cell_snapshots(
    reader: multicellds.MultiCellDS,
    n_frames: int,
    color_col: str,
    frame_spacing: str = "linear",
) -> list[dict]:
    """
    Load per-cell data for N evenly-spaced time frames.

    Returns a list of dicts with keys:
        time_h, x, y, total_volume, <color_col>
    (each value is a 1-D numpy array over all cells in that snapshot).
    """
    cols     = _get_cell_columns(reader._output_folder)
    name2idx = {n: i for i, n in enumerate(cols)}

    needed = ["x_position", "y_position", "total_volume", color_col]
    missing = [c for c in needed if c not in name2idx]
    if missing:
        raise ValueError(
            f"Required columns not found in cell matrix: {missing}\n"
            f"Available columns include: {cols[:30]}"
        )

    all_items = list(reader.cells_as_matrix_iterator())
    if not all_items:
        raise RuntimeError("No cell data found in output directory.")

    n_total = len(all_items)
    if n_total <= 1:
        indices = np.array([0], dtype=int)
    elif frame_spacing == "log":
        # Denser sampling in the early (fast-changing) period.
        raw     = np.geomspace(1, n_total, min(n_frames, n_total))
        indices = np.unique(np.clip(np.round(raw).astype(int) - 1, 0, n_total - 1))
    else:
        # Linear: evenly spaced across the full simulation time.
        indices = np.unique(np.linspace(0, n_total - 1, min(n_frames, n_total)).astype(int))

    snapshots = []
    for idx in indices:
        t_min, a = all_items[idx]
        a        = _ensure_rows(a, len(cols))
        time_h   = t_min / 60.0
        snapshots.append({
            "time_h":       time_h,
            "x":            a[:, name2idx["x_position"]],
            "y":            a[:, name2idx["y_position"]],
            "total_volume": a[:, name2idx["total_volume"]],
            color_col:      a[:, name2idx[color_col]],
        })

    print(f"  Loaded {len(snapshots)} cell snapshots "
          f"(t = {snapshots[0]['time_h']:.2f} – {snapshots[-1]['time_h']:.2f} h)")
    return snapshots


# ===========================================================================
# Main figure
# ===========================================================================

def plot_composite(
    output_dir:        str,
    out_stem:          str,
    color_by:          str  = "velocity",
    n_frames:          int  = 12,
    glc_frame:         str  = "last",
    dpi:               int  = 300,
    cell_radius:       float | None = None,
    x_scale:           str  = "linear",
    show_glc_bg:       bool  = True,
    figsize:           tuple = (7.09, 2.4),  # inches: Nature double-column
    frame_spacing:     str  = "linear",
    glc_cmap:          str  = "Blues",
    cell_cmap:         str | None = None,
    cell_edgecolor:    str  = "#1b4332",
    cell_linewidth:    float = 0.6,
    cell_alpha:        float = 0.97,
    glc_alpha:         float = 0.85,
    cell_cmap_min:     float = 0.15,
    cell_cmap_max:     float = 0.95,
) -> None:
    """
    Build and save the composite figure.

    Parameters
    ----------
    output_dir    : PhysiCell output folder
    out_stem      : output path without extension (both .svg and .pdf are saved)
    color_by      : "velocity" | "growth_rate"
    n_frames      : number of time snapshots to overlay
    glc_frame     : "last" | "first"  – which frame supplies the glucose background
    dpi           : output resolution
    cell_radius   : display radius of each cell circle (µm).
                    None = physical radius per cell from total_volume: r = (3V/4π)^(1/3).
    x_scale       : "linear" | "symlog"
    show_glc_bg   : whether to show the glucose concentration background
    figsize       : figure size in inches (width, height)
    glc_cmap      : background colormap; recommended: "Blues" or "Greys"
    cell_cmap     : cell colormap; default uses metric-specific mapping
    cell_edgecolor: cell outline color for contrast against the background
    cell_linewidth: cell outline width
    cell_alpha    : cell fill opacity
    glc_alpha     : glucose background opacity
    cell_cmap_min : lower fraction of cell colormap to use
    cell_cmap_max : upper fraction of cell colormap to use
    """
    if color_by not in _METRIC_COL:
        raise ValueError(f"--color-by must be one of {list(_METRIC_COL)}")

    color_col   = _METRIC_COL[color_by]
    color_label = _METRIC_LABEL[color_by]
    color_cmap  = cell_cmap if cell_cmap is not None else _METRIC_CMAP[color_by]

    # -----------------------------------------------------------------------
    # Load data
    # -----------------------------------------------------------------------
    print(f"[composite] Loading output from: {output_dir}")
    reader   = multicellds.MultiCellDS(output_dir)
    domain_x, domain_y = _get_domain_bounds(reader)
    print(f"  Domain: x {domain_x}, y {domain_y}")

    snapshots = _load_cell_snapshots(reader, n_frames, color_col, frame_spacing)

    # Determine colour scale across all frames so the same scale is used
    # for every circle in the composite plot.
    all_metric = np.concatenate([s[color_col] for s in snapshots])
    # Linear norm: migration speed fraction is in [0, 1]
    c_vmin = 0.0
    c_vmax = 1.0

    norm_cells = Normalize(vmin=c_vmin, vmax=c_vmax)
    cmap_cells = _trimmed_cmap(color_cmap, cell_cmap_min, cell_cmap_max)

    Gi = glc_vmin = glc_vmax = None
    if show_glc_bg:
        print(f"  Loading glucose background from {glc_frame} frame …")
        Gi, glc_vmin, glc_vmax = _load_glucose_background(
            reader, glc_frame, domain_x, domain_y,
            grid_nx=700, grid_ny=200,
        )
        print(f"  Glucose colour range: [{glc_vmin:.4f}, {glc_vmax:.4f}] mM")

    # -----------------------------------------------------------------------
    # Apply publication style
    # -----------------------------------------------------------------------
    _apply_publication_style()

    # -----------------------------------------------------------------------
    # Figure layout
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    # Two narrow colorbars on the right; constrained_layout handles spacing.
    if show_glc_bg:
        gs = fig.add_gridspec(
            1, 3,
            width_ratios=[1, 0.04, 0.04],
        )
        ax       = fig.add_subplot(gs[0, 0])
        cax_glc  = fig.add_subplot(gs[0, 1])
        cax_cell = fig.add_subplot(gs[0, 2])
    else:
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.04])
        ax       = fig.add_subplot(gs[0, 0])
        cax_cell = fig.add_subplot(gs[0, 1])
        cax_glc  = None

    # -----------------------------------------------------------------------
    # Glucose background
    # -----------------------------------------------------------------------
    if show_glc_bg and Gi is not None:
        norm_glc  = Normalize(vmin=glc_vmin, vmax=glc_vmax)
        cmap_glc  = mpl.colormaps[glc_cmap]
        im = ax.imshow(
            Gi,
            origin="lower",
            extent=[domain_x[0], domain_x[1], domain_y[0], domain_y[1]],
            vmin=glc_vmin,
            vmax=glc_vmax,
            cmap=cmap_glc,
            interpolation="bilinear",
            aspect="auto",
            alpha=glc_alpha,
            zorder=1,
        )
        cb_glc = fig.colorbar(ScalarMappable(norm=norm_glc, cmap=cmap_glc),
                               cax=cax_glc)
        cb_glc.set_label("Glucose (mM)", labelpad=3)
        cb_glc.ax.tick_params(labelsize=5, length=2)
        cb_glc.outline.set_linewidth(0.4)


    # -----------------------------------------------------------------------
    # Cell circles (one circle per cell per selected frame)
    # -----------------------------------------------------------------------
    # Colour time annotations alternately above/below to avoid crowding.
    n_snap   = len(snapshots)
    cmap_t   = mpl.colormaps["Greys"].resampled(n_snap + 2)  # light grey → dark grey for time labels
    edgecolors = ["#222222"] * n_snap
    # Minimum x-separation (data units) between consecutive top labels.
    _label_min_sep = (domain_x[1] - domain_x[0]) / max(n_snap, 1) * 0.5
    _last_label_x  = -np.inf

    for k, snap in enumerate(snapshots):
        xs       = snap["x"]
        ys       = snap["y"]
        vols     = snap["total_volume"]
        metric   = snap[color_col]
        time_h   = snap["time_h"]

        for i in range(len(xs)):
            face_color = cmap_cells(norm_cells(metric[i]))
            _r = cell_radius if cell_radius is not None \
                else (3.0 * vols[i] / (4.0 * np.pi))**(1.0 / 3.0)
            # Draw a thin white halo under the cell and then a dark outline.
            # This keeps both low- and high-speed cells visible on any background.
            halo = Circle(
                (xs[i], ys[i]),
                radius=_r * 1.08,
                facecolor="none",
                edgecolor="white",
                linewidth=max(cell_linewidth * 1.6, 0.8),
                alpha=0.95,
                zorder=3 + k,
            )
            ax.add_patch(halo)

            circ = Circle(
                (xs[i], ys[i]),
                radius=_r,
                facecolor=face_color,
                edgecolor=cell_edgecolor,
                linewidth=cell_linewidth,
                alpha=cell_alpha,
                zorder=4 + k,
            )
            ax.add_patch(circ)

        # Time label above the top spine, at the cell centroid x position.
        # get_xaxis_transform() blends data-x with axes-fraction-y.
        cx_mean = xs.mean()
        if abs(cx_mean - _last_label_x) >= _label_min_sep or k == 0 or k == n_snap - 1:
            trans = ax.get_xaxis_transform()
            ax.plot([cx_mean, cx_mean], [1.0, 1.025], transform=trans,
                    color="#888888", lw=0.5, clip_on=False, zorder=11)
            ax.text(cx_mean, 1.03, f"{time_h:.1f} h",
                    transform=trans,
                    ha="center", va="bottom", fontsize=5,
                    color="#333333", clip_on=False, zorder=12)
            _last_label_x = cx_mean

    # -----------------------------------------------------------------------
    # Cell colorbar
    # -----------------------------------------------------------------------
    sm_cells = ScalarMappable(norm=norm_cells, cmap=cmap_cells)
    sm_cells.set_array([])
    cb_cells = fig.colorbar(sm_cells, cax=cax_cell)
    cb_cells.set_label(color_label, labelpad=3)
    cb_cells.ax.tick_params(labelsize=5, length=2)
    cb_cells.outline.set_linewidth(0.4)

    # -----------------------------------------------------------------------
    # Axes formatting
    # -----------------------------------------------------------------------
    # Use the full domain extent so the complete glucose gradient is visible.
    ax.set_xlim(domain_x[0], domain_x[1])
    ax.set_ylim(domain_y[0] - 2, domain_y[1] + 2)
    # Equal aspect so cells (drawn in data units) appear as circles, not ovals.
    ax.set_aspect("equal", adjustable="box")

    if x_scale == "symlog":
        # linthresh: linear region around 0 of width ~20 µm
        ax.set_xscale("symlog", linthresh=20, linscale=0.5)
        ax.xaxis.set_minor_locator(ticker.NullLocator())

    ax.set_xlabel("Position along gradient axis (µm)")
    ax.set_ylabel("y position (µm)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["left"].set_linewidth(0.6)
    ax.tick_params(axis="both", which="major", labelsize=6, length=3, width=0.6)

    # Legend


    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    os.makedirs(os.path.dirname(out_stem) if os.path.dirname(out_stem) else ".", exist_ok=True)

    for ext in ("svg", "png"):
        fpath = f"{out_stem}.{ext}"
        fig.savefig(fpath, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
        print(f"  Saved: {fpath}")

    plt.close(fig)


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Composite publication figure: cell trajectory + glucose gradient.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--output-dir",  default=DEFAULT_OUTPUT_DIR,
                   help="PhysiCell output directory")
    p.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                   help="Directory for output figures")
    p.add_argument("--color-by",    default="velocity",
                   choices=["velocity", "growth_rate"],
                   help="Cell colour metric")
    p.add_argument("--n-frames",    type=int, default=12,
                   help="Number of time snapshots to overlay")
    p.add_argument("--dpi",         type=int, default=300,
                   help="Output resolution in DPI")
    p.add_argument("--cell-radius", type=float, default=None,
                   help="Display radius of each cell (µm). Default: physical radius from total_volume.")
    p.add_argument("--x-scale",     default="linear",
                   choices=["linear", "symlog"],
                   help="X-axis scaling")
    p.add_argument("--glc-frame",   default="last",
                   choices=["last", "first"],
                   help="Which frame to use for the glucose background")
    p.add_argument("--no-glc-bg",   action="store_true",
                   help="Disable the glucose concentration background")
    p.add_argument("--figsize",     nargs=2, type=float, default=[7.09, 2.4],
                   metavar=("W", "H"),
                   help="Figure size in inches (width height)")
    p.add_argument("--frame-spacing", default="linear",
                   choices=["linear", "log"],
                   help="Frame sampling: 'linear'=evenly spaced, 'log'=denser early frames")
    p.add_argument("--glc-cmap", default="Blues",
                   help="Glucose background colormap. Good options: Blues, Greys, cividis")
    p.add_argument("--cell-cmap", default=None,
                   help="Cell metric colormap. Default: summer for velocity, Greens for growth_rate")
    p.add_argument("--cell-edgecolor", default="#1b4332",
                   help="Cell outline color. Default is dark green rather than black.")
    p.add_argument("--cell-linewidth", type=float, default=0.6,
                   help="Cell outline width")
    p.add_argument("--cell-alpha", type=float, default=0.97,
                   help="Cell fill opacity")
    p.add_argument("--cell-cmap-min", type=float, default=0.15,
                   help="Lower fraction of the cell colormap to use; avoids near-black/near-white extremes")
    p.add_argument("--cell-cmap-max", type=float, default=0.95,
                   help="Upper fraction of the cell colormap to use; avoids near-black/near-white extremes")
    p.add_argument("--glc-alpha", type=float, default=0.85,
                   help="Glucose background opacity")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)

    out_stem = os.path.join(
        args.results_dir,
        f"cell_migration_composite_{args.color_by}",
    )

    plot_composite(
        output_dir    = args.output_dir,
        out_stem      = out_stem,
        color_by      = args.color_by,
        n_frames      = args.n_frames,
        glc_frame     = args.glc_frame,
        dpi           = args.dpi,
        cell_radius   = args.cell_radius,
        x_scale       = args.x_scale,
        show_glc_bg   = not args.no_glc_bg,
        figsize       = tuple(args.figsize),
        frame_spacing = args.frame_spacing,
        glc_cmap      = args.glc_cmap,
        cell_cmap     = args.cell_cmap,
        cell_edgecolor= args.cell_edgecolor,
        cell_linewidth= args.cell_linewidth,
        cell_alpha    = args.cell_alpha,
        glc_alpha     = args.glc_alpha,
        cell_cmap_min = args.cell_cmap_min,
        cell_cmap_max = args.cell_cmap_max,
    )


if __name__ == "__main__":
    main()
