# ABOUTME: Human-comprehension CLI. Runs a Problem through the solver/harness
# once and emits a visual dashboard (mesh, source, kappa, exact T, computed
# T_h, error field) plus a convergence plot. Outputs land under
# artifacts/inspect/<problem_name>/. Not exercised by pytest - this is the
# look-at-your-work tool, separate from the failure-only artifacts under
# tests/_artifacts/.
#
# Usage:
#   uv run python -m scripts.inspect problem_01
#   uv run python -m scripts.inspect problem_01 --mesh-size 0.05 --no-convergence

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # static PNG only.
import matplotlib.pyplot as plt
import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.harness.study import run_refinement_study  # noqa: E402
from src.problems.problem_01_manufactured import Problem01Manufactured  # noqa: E402
from src.problems.problem_02_slab import Problem02Slab  # noqa: E402
from src.solver.solve_scalar import solve_scalar  # noqa: E402

logger = logging.getLogger("inspect")


# ============================================================================
# CONFIG
# ============================================================================

ARTIFACT_ROOT = _PROJECT_ROOT / "artifacts" / "inspect"

# Registry. New problems land here when their submission accepts.
PROBLEMS: dict[str, type] = {
    "problem_01": Problem01Manufactured,
    "problem_02": Problem02Slab,
}


# ============================================================================
# PANELS
# ============================================================================


# Labeling convention: Part 1 verification problems are non-dimensional
# (physics.md § Symbol glossary, Q row). x, y, T, kappa, Q, and the error norms
# all carry the dimensionless marker [-]. The y-axis label "error norm [-]"
# on the convergence plot stays generic because L^2 (K m) and H^1-seminorm (K)
# have different physical dimensions in general - this is documented next to
# _panel_convergence.
_DIM = r"$[-]$"


def _label_axes(ax) -> None:
    """Standard x/y axis labels for non-dimensional 2D verification panels."""
    ax.set_xlabel(rf"$x$  {_DIM}")
    ax.set_ylabel(rf"$y$  {_DIM}")


def _panel_mesh(ax, mesh, pin_dof: int | None) -> None:
    """Triangulation + pinned-DOF marker + boundary-edge color by tag."""
    ax.triplot(mesh.p[0], mesh.p[1], mesh.t.T, linewidth=0.3, color="0.6")
    # Boundary edges, colored per physical tag. Names come from gmsh physical
    # group names (e.g. "bottom", "left") so the legend reads the way you
    # built the geometry.
    palette = plt.get_cmap("tab10")
    if getattr(mesh, "boundaries", None):
        for i, (name, facets) in enumerate(mesh.boundaries.items()):
            for facet in facets:
                a, b = mesh.facets[:, facet]
                ax.plot(
                    [mesh.p[0, a], mesh.p[0, b]],
                    [mesh.p[1, a], mesh.p[1, b]],
                    color=palette(i % 10),
                    linewidth=1.6,
                    label=name if facet == facets[0] else None,
                )
    if pin_dof is not None:
        ax.plot(
            mesh.p[0, pin_dof], mesh.p[1, pin_dof],
            marker="*", markersize=14, color="red", markeredgecolor="black",
            label=f"pin (DOF {pin_dof})",
        )
    ax.set_aspect("equal")
    ax.set_title(
        f"mesh — {mesh.t.shape[1]} elements, {mesh.p.shape[1]} nodes\n"
        f"(boundary tags colored; star = nullspace pin)"
    )
    _label_axes(ax)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.85)


def _panel_field(
    ax, mesh, values, *, title: str, cbar_label: str,
    cmap: str = "viridis", symmetric: bool = False,
) -> None:
    """tricontourf of a nodal field with labeled axes and colorbar."""
    if symmetric:
        vmax = float(np.max(np.abs(values))) or 1.0
        levels = np.linspace(-vmax, vmax, 21)
        kw = {"cmap": cmap, "levels": levels}
    else:
        kw = {"cmap": cmap, "levels": 21}
    tcf = ax.tricontourf(mesh.p[0], mesh.p[1], mesh.t.T, values, **kw)
    cbar = plt.colorbar(tcf, ax=ax, shrink=0.85)
    cbar.set_label(cbar_label)
    ax.set_aspect("equal")
    ax.set_title(title)
    _label_axes(ax)


def _panel_kappa(ax, mesh, problem) -> None:
    """Per-element kappa as P0 field via tripcolor.

    Uniformity-aware colorbar: when min(kappa) == max(kappa), matplotlib's
    default auto-range would expand to ~+/-10% of the value and the constant
    field would *look* variable. We detect uniformity and either pin a tight
    range around the value (and annotate "uniform") or use the data range
    when kappa actually varies (Problem 2+).
    """
    subdomains = getattr(mesh, "subdomains", None) or {}
    n_elem = mesh.t.shape[1]
    kappa = np.full(n_elem, np.nan)
    # Walk real surface tags; skip gmsh metadata leaks. Dispatch by name to
    # mirror src/solver/solve_scalar.py:_build_kappa_p0_values - Problem.kappa
    # takes the subdomain name string, not an integer tag.
    for name, idx in subdomains.items():
        if name.startswith("gmsh:"):
            continue
        kappa[idx] = float(problem.kappa(name))
    if np.isnan(kappa).any():
        # No real subdomain tags exported (single-region geometries with no
        # physical groups). Fall back to a single Problem.kappa("") call,
        # mirroring the solver's no-tags branch.
        kappa = np.nan_to_num(kappa, nan=float(problem.kappa("")))

    kmin, kmax = float(kappa.min()), float(kappa.max())
    is_uniform = np.isclose(kmin, kmax)
    if is_uniform:
        # Tight symmetric pad around the value so the colorbar reads as a
        # single tick, not a fake gradient.
        pad = max(abs(kmin) * 1e-3, 1e-12)
        vmin, vmax = kmin - pad, kmin + pad
    else:
        vmin, vmax = kmin, kmax

    tpc = ax.tripcolor(
        mesh.p[0], mesh.p[1], mesh.t.T, facecolors=kappa, edgecolors="none",
        cmap="cividis", vmin=vmin, vmax=vmax,
    )
    cbar = plt.colorbar(tpc, ax=ax, shrink=0.85)
    cbar.set_label(rf"$\kappa$  {_DIM}")
    if is_uniform:
        # Replace the (now meaningless) tight tick range with a single label.
        cbar.set_ticks([kmin])
        cbar.set_ticklabels([f"{kmin:.4g}"])
    ax.set_aspect("equal")
    span = (
        f"uniform $\\kappa = {kmin:.3g}$"
        if is_uniform
        else f"range $\\kappa \\in [{kmin:.3g}, {kmax:.3g}]$"
    )
    ax.set_title(f"conductivity $\\kappa(x, y)$  (P0 per element)\n{span}")
    _label_axes(ax)


def _panel_convergence(ax, study, problem_name: str) -> None:
    # Units: Part 1 verification problems are non-dimensional by stated
    # convention (physics.md § Symbol glossary, Q row: "in practice
    # dimensionless in the warm-up tests"). The unit square has no physical
    # length scale, so `h`, T, and the error norms here are all dimensionless.
    # L^2 and H^1-seminorm have *different* dimensions in the physical case
    # (K m vs K), so the y-axis stays generic and the legend names each norm.
    h = np.array([lvl.mesh_size for lvl in study.levels])
    l2 = np.array([lvl.l2_error for lvl in study.levels])
    h1 = np.array([lvl.h1_error for lvl in study.levels])
    ax.loglog(h, l2, "o-", label=rf"$\|T_h - T\|_{{L^2}}$ (fit {study.l2_rate:.2f})")
    ax.loglog(h, h1, "s-", label=rf"$|T_h - T|_{{H^1}}$ (fit {study.h1_rate:.2f})")
    # Reference slopes anchored at the coarsest point.
    ax.loglog(h, l2[0] * (h / h[0]) ** 2, "--", color="0.5", label="slope 2")
    ax.loglog(h, h1[0] * (h / h[0]) ** 1, ":", color="0.5", label="slope 1")
    ax.set_xlabel(r"$h$  (characteristic mesh length, $[-]$)")
    ax.set_ylabel(r"error norm  $[-]$")
    ax.set_title(f"convergence  |  {problem_name}  (non-dimensional verification)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)


# ============================================================================
# DASHBOARD
# ============================================================================


def _interface_x_or_none(problem) -> float | None:
    """Return the vertical interface x-coordinate, if the problem has one.

    Problem 2 has a derivative kink at x=1/2; the error panel annotates it
    so the visible structure reads as theory-confirming (inspector.md Rule
    6). When a third problem grows its own interface, replace this isinstance
    branch with a Protocol-level accessor; one if-statement is not worth a
    hook system yet.
    """
    if isinstance(problem, Problem02Slab):
        from src.geometry.rectangle_split import _X_INTERFACE  # local import

        return float(_X_INTERFACE)
    return None


def _two_panel_layout(mesh) -> tuple[tuple[int, int], tuple[float, float]]:
    """Pick (nrows, ncols) and figsize so a two-panel figure stays legible.

    Wide-and-flat geometries (Problem 2's 1x0.1 slab) are unreadable when
    panels sit side-by-side: ax.set_aspect("equal") crushes them next to
    their colorbars. Stacking vertically gives each panel the full figure
    width. For roughly-square geometries (Problem 1's unit square),
    side-by-side reads more naturally.
    """
    xspan = float(np.ptp(mesh.p[0]))
    yspan = float(np.ptp(mesh.p[1]))
    if xspan / max(yspan, 1e-12) > 3.0:
        return (2, 1), (11.0, 5.5)
    return (1, 2), (12.0, 5.0)


def _save(fig, out: Path) -> None:
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out)


def _emit_dashboard(out_dir: Path, problem, mesh_size: float) -> None:
    """Run one solve at ``mesh_size`` and render three two-panel PNGs.

    Split into setup (mesh + kappa), problem statement (source + exact),
    and result (computed + pointwise error). The 2x3 single-figure layout
    is unreadable for wide-flat geometries; the split also reads better
    when a panel is shown standalone in a writeup.
    """
    sr = solve_scalar(problem, mesh_size)
    mesh = sr.basis.mesh
    x, y = mesh.p[0], mesh.p[1]
    source_vals = problem.source(x, y)
    exact_vals = problem.exact_solution(x, y)
    err_vals = sr.solution - exact_vals
    err_max = float(np.max(np.abs(err_vals)))

    shape, figsize = _two_panel_layout(mesh)
    base_suptitle = f"{problem.name}  |  h = {mesh_size:.4g}  |  non-dimensional verification"

    # --- 01_setup: mesh + kappa --------------------------------------------
    fig, axes = plt.subplots(*shape, figsize=figsize)
    _panel_mesh(axes[0], mesh, sr.pin_dof)
    _panel_kappa(axes[1], mesh, problem)
    fig.suptitle(f"setup  ::  {base_suptitle}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, out_dir / "01_setup.png")

    # --- 02_solution: exact vs computed ------------------------------------
    # Same colormap and (implicit) same value range - the comparison is
    # apples-to-apples by construction since T_h is converging to T.
    fig, axes = plt.subplots(*shape, figsize=figsize)
    _panel_field(
        axes[0], mesh, exact_vals,
        title=r"exact $T(x, y)$  (manufactured solution)",
        cbar_label=rf"$T$  {_DIM}",
        cmap="viridis",
    )
    _panel_field(
        axes[1], mesh, sr.solution,
        title=r"computed $T_h(x, y)$  (P1 FEM)",
        cbar_label=rf"$T_h$  {_DIM}",
        cmap="viridis",
    )
    fig.suptitle(f"solution  ::  {base_suptitle}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, out_dir / "02_solution.png")

    # --- 03_diagnostic: source + pointwise error ---------------------------
    # Both are scalar fields with similar visual weight; pairing them puts
    # "what was the PDE driving" next to "how far off is the discretization".
    fig, axes = plt.subplots(*shape, figsize=figsize)
    _panel_field(
        axes[0], mesh, source_vals,
        title=r"source $Q(x, y)$  (PDE right-hand side)",
        cbar_label=rf"$Q$  {_DIM}",
        cmap="magma",
    )
    _panel_field(
        axes[1], mesh, err_vals,
        title=r"pointwise error $T_h - T$",
        cbar_label=rf"$T_h - T$  {_DIM}",
        cmap="RdBu_r", symmetric=True,
    )
    # Annotate the error panel: star at the pinned DOF (error = 0 by
    # construction) and 'x' at the |error| max. The visible asymmetry of the
    # field is the pin doing its job - without these markers the reader has
    # to recompute it.
    err_panel = axes[1]
    if sr.pin_dof is not None:
        err_panel.plot(
            mesh.p[0, sr.pin_dof], mesh.p[1, sr.pin_dof],
            marker="*", markersize=14, color="white", markeredgecolor="black",
            linestyle="none", label=f"pin (DOF {sr.pin_dof}, err=0)",
        )
    imax = int(np.argmax(np.abs(err_vals)))
    err_panel.plot(
        mesh.p[0, imax], mesh.p[1, imax],
        marker="x", markersize=10, color="black", markeredgewidth=2,
        linestyle="none", label=f"max |err| = {abs(err_vals[imax]):.2e}",
    )
    # Problem-specific Rule-6 annotation: Problem 2's exact solution has a
    # derivative kink at x=1/2 where the parabolic left region meets the
    # constant right region. The P1 error field shows structure aligned with
    # the interface; drawing the line tells the reader "this is the kink,
    # not a bug" before they propose investigating it.
    interface_x = _interface_x_or_none(problem)
    if interface_x is not None:
        y_lo, y_hi = float(mesh.p[1].min()), float(mesh.p[1].max())
        err_panel.plot(
            [interface_x, interface_x], [y_lo, y_hi],
            color="black", linestyle="--", linewidth=1.0,
            label=rf"interface $x={interface_x:g}$ (kink)",
        )
    err_panel.legend(loc="upper right", fontsize=7, framealpha=0.85)
    fig.suptitle(
        f"diagnostic  ::  {base_suptitle}  |  $\\|T_h - T\\|_\\infty$ = {err_max:.3e}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, out_dir / "03_diagnostic.png")


def _emit_convergence(out_dir: Path, problem) -> None:
    """Run the full refinement study and render the convergence plot."""
    study = run_refinement_study(problem)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    _panel_convergence(ax, study, problem.name)
    fig.tight_layout()
    out = out_dir / "convergence.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out)


# ============================================================================
# CLI
# ============================================================================


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)-7s %(name)s :: %(message)s"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "problem", choices=sorted(PROBLEMS), help="Problem name to inspect."
    )
    parser.add_argument(
        "--mesh-size", type=float, default=None,
        help="Single mesh size for the dashboard. Default: finest from Problem.",
    )
    parser.add_argument(
        "--no-convergence", action="store_true",
        help="Skip the refinement study and convergence plot.",
    )
    args = parser.parse_args()

    problem = PROBLEMS[args.problem]()
    mesh_size = args.mesh_size if args.mesh_size is not None else problem.mesh_sizes()[-1]
    out_dir = ARTIFACT_ROOT / args.problem
    out_dir.mkdir(parents=True, exist_ok=True)

    _emit_dashboard(out_dir, problem, mesh_size)
    if not args.no_convergence:
        _emit_convergence(out_dir, problem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
