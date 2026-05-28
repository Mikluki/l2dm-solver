# ABOUTME: Layer 2 solver invariant checker. Runs solve_scalar(), re-assembles
# K and b using the *same* forms the solver uses (an honest shadow check —
# artifacts/inspect/CONVENTIONS.md § Technique 2), and writes
# artifacts/inspect/<problem>/internals.md with three sections:
# subdomain integrity, source verification, linear-system invariants.
#
# It intentionally imports the underscore-prefixed forms and helpers from
# src.solver.solve_scalar. They are private to the solver as a module, but
# diagnose is a sibling diagnostic tool that must use the literal same forms
# the solve does — otherwise the "honest shadow" property is lost and a
# load-assembly bug could hide. The privacy is a code-organisation hint, not
# a contract.
#
# Usage:
#   uv run python -m scripts.diagnose problem_01
#   uv run python -m scripts.diagnose problem_01 --mesh-size 0.1

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
from skfem import ElementTriP0, Functional

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.problems.problem_01_manufactured import Problem01Manufactured  # noqa: E402
from src.solver.solve_scalar import (  # noqa: E402
    _RESERVED_SUBDOMAIN_NAMES,
    _build_kappa_p0_values,
    _make_load,
    _resolve_dirichlet_dofs,
    _stiffness_kappa,
    solve_scalar,
)

logger = logging.getLogger("diagnose")


# ============================================================================
# CONFIG
# ============================================================================

ARTIFACT_ROOT = _PROJECT_ROOT / "artifacts" / "inspect"

# Problem registry. Mirrors scripts/inspect.py; kept separate so the two
# tools can diverge if needed (e.g. inspect drops a problem before diagnose
# learns about it, or vice versa).
PROBLEMS: dict[str, type] = {
    "problem_01": Problem01Manufactured,
}

# Thresholds per artifacts/inspect/CONVENTIONS.md § Layer 2.
_SYMMETRY_TOL = 1e-14
_RESIDUAL_TOL = 1e-10  # spsolve floor; loosened from the doc's ~1e-12 to
# tolerate larger right-hand-side magnitudes (Problem 1's b has entries ~10).
_SOURCE_TOL_FACTOR = 100.0  # |Δ| ≤ factor · h⁴ for smooth Q under P1 quadrature


# ============================================================================
# HELPERS
# ============================================================================


def _element_areas(mesh) -> np.ndarray:
    """Triangle areas via |det([p1-p0, p2-p0])| / 2."""
    p, t = mesh.p, mesh.t
    v0, v1, v2 = p[:, t[0]], p[:, t[1]], p[:, t[2]]
    e1 = v1 - v0
    e2 = v2 - v0
    return 0.5 * np.abs(e1[0] * e2[1] - e1[1] * e2[0])


def _tick(passed: bool | None) -> str:
    if passed is None:
        return "(skipped)"
    return "✓" if passed else "✗"


# ============================================================================
# TECHNIQUES
# ============================================================================


def _subdomain_integrity(mesh, problem) -> tuple[bool, list[str]]:
    """Technique 1: subdomain integrity report.

    Walks ``mesh.subdomains``, names starting with ``gmsh:`` or in the solver's
    reserved set are skipped. For each real subdomain, reports n_elem, area,
    fraction, and the κ value the Problem assigns by name. Verdict ✓ iff
    real-subdomain areas cover the domain to roundoff and every real
    subdomain has a Problem-resolvable κ.
    """
    areas = _element_areas(mesh)
    domain_area = float(areas.sum())
    rows: list[str] = []
    real_area = 0.0
    real_count = 0
    kappa_lookup_failed = False

    if not mesh.subdomains:
        # Single untagged region (Problem 1 with no surface tag exported).
        try:
            kval = float(problem.kappa(""))
            kstr = f"{kval:.4g}"
        except Exception as exc:  # noqa: BLE001
            kstr = f"ERR: {exc}"
            kappa_lookup_failed = True
        rows.append(
            f"| (untagged) | {mesh.t.shape[1]} | {domain_area:.4f} | 100.0% | "
            f"{kstr} | no |"
        )
        real_area = domain_area
        real_count = 1
    else:
        for name in mesh.subdomains:
            idx = np.asarray(mesh.subdomains[name])
            is_meta = name in _RESERVED_SUBDOMAIN_NAMES or name.startswith("gmsh:")
            if is_meta:
                rows.append(f"| `{name}` | — | — | — | — | yes (skipped) |")
                continue
            n_elem = int(idx.size)
            area = float(areas[idx].sum())
            fraction = area / domain_area if domain_area else float("nan")
            try:
                kval = float(problem.kappa(name))
                kstr = f"{kval:.4g}"
            except Exception as exc:  # noqa: BLE001
                kstr = f"ERR: {exc}"
                kappa_lookup_failed = True
            rows.append(
                f"| `{name}` | {n_elem} | {area:.4f} | {100*fraction:.1f}% | "
                f"{kstr} | no |"
            )
            real_area += area
            real_count += 1

    delta = abs(real_area - domain_area)
    passed = (delta < 1e-10) and (real_count >= 1) and not kappa_lookup_failed
    rows.append("")
    rows.append(
        f"**Coverage:** {real_area:.6f} / {domain_area:.6f}  "
        f"(Δ = {delta:.2e})  {_tick(passed)}"
    )
    if real_count < 1:
        rows.append("✗ No real subdomain found.")
    return passed, rows


def _source_verification(
    basis, problem, mesh_size: float
) -> tuple[bool | None, list[str]]:
    """Technique 2: ∫Q dA two ways.

    Returns (None, [skip-reason]) when the Problem doesn't declare
    ``source_integral()``. See CONVENTIONS § Technique 2 gotcha — a missing
    declaration is the *honest* outcome when the integral has no obviously-
    correct one-line form.
    """
    if not hasattr(problem, "source_integral"):
        return None, [
            "Skipped: problem does not declare `source_integral()`.",
            "(See CONVENTIONS § Technique 2 — declare only if a one-line",
            "closed form exists.)",
        ]

    analytic = float(problem.source_integral())

    @Functional
    def _q_int(w):
        return problem.source(w.x[0], w.x[1])

    fe_val = float(_q_int.assemble(basis))
    diff = abs(fe_val - analytic)
    floor = _SOURCE_TOL_FACTOR * mesh_size ** 4
    passed = diff <= floor

    return passed, [
        f"- Analytic ∫Q dA:        `{analytic:+.6e}`",
        f"- FE quadrature ∫Q dA:   `{fe_val:+.6e}`",
        f"- |Δ|:                   `{diff:.2e}`",
        f"- Threshold (100·h⁴):    `{floor:.2e}`",
        f"- Verdict: {_tick(passed)}",
    ]


def _linear_system_invariants(
    basis, K, b, u_h, pin_dof, dirichlet_dofs
) -> tuple[bool, list[str]]:
    """Technique 3: matrix symmetry + linear residual on the *unsubstituted* rows.

    The residual is checked only on rows whose equation was actually enforced
    by spsolve. The pin row and any Dirichlet rows were replaced by the
    substitution ``u_i = value`` during ``condense``; their original PDE
    equation ``K[i, :] @ u = b[i]`` is not solved and its residual lives at
    problem scale, not solver scale. Including those rows would mean the
    threshold is testing nothing useful.
    """
    sym_num = float(np.sqrt(float((K - K.T).multiply(K - K.T).sum())))
    sym_den = float(np.sqrt(float(K.multiply(K).sum())))
    sym_rel = sym_num / sym_den if sym_den else float("nan")

    enforced = np.ones(basis.N, dtype=bool)
    if pin_dof is not None:
        enforced[pin_dof] = False
    if dirichlet_dofs.size:
        enforced[dirichlet_dofs] = False
    full_residual = K @ u_h - b
    residual = float(np.max(np.abs(full_residual[enforced])))
    n_enforced = int(enforced.sum())
    n_substituted = basis.N - n_enforced

    sym_pass = sym_rel < _SYMMETRY_TOL
    res_pass = residual < _RESIDUAL_TOL
    passed = sym_pass and res_pass

    rows = [
        f"- `‖K − Kᵀ‖_F / ‖K‖_F`:  `{sym_rel:.2e}`  "
        f"(threshold `{_SYMMETRY_TOL:.0e}`)  {_tick(sym_pass)}",
        f"- `‖A·u − b‖_∞` (enforced rows): `{residual:.2e}`  "
        f"(threshold `{_RESIDUAL_TOL:.0e}`)  {_tick(res_pass)}",
        f"- DOFs (total / enforced / substituted): {basis.N} / "
        f"{n_enforced} / {n_substituted}",
    ]
    if pin_dof is not None:
        x_pin = float(basis.mesh.p[0, pin_dof])
        y_pin = float(basis.mesh.p[1, pin_dof])
        rows.append(f"- Pinned DOF: {pin_dof} at ({x_pin:.3f}, {y_pin:.3f})")
    if dirichlet_dofs.size:
        rows.append(f"- Dirichlet DOFs: {dirichlet_dofs.size}")
    return passed, rows


# ============================================================================
# REPORT
# ============================================================================


def diagnose(problem_name: str, mesh_size: float | None = None) -> Path:
    problem = PROBLEMS[problem_name]()
    if mesh_size is None:
        mesh_size = float(problem.mesh_sizes()[-1])
    out_dir = ARTIFACT_ROOT / problem_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Solve --------------------------------------------------------------
    sr = solve_scalar(problem, mesh_size)
    mesh = sr.basis.mesh

    # --- Re-assemble K and b with the same forms solve_scalar uses ----------
    basis_p0 = sr.basis.with_element(ElementTriP0())
    kappa_values = _build_kappa_p0_values(mesh, problem)
    kappa_field = basis_p0.interpolate(kappa_values)
    K = _stiffness_kappa.assemble(sr.basis, kappa=kappa_field)
    b = _make_load(problem.source).assemble(sr.basis)

    # --- Run the three Layer 2 techniques -----------------------------------
    sub_pass, sub_rows = _subdomain_integrity(mesh, problem)
    src_pass, src_rows = _source_verification(sr.basis, problem, mesh_size)
    dirichlet_dofs, _ = _resolve_dirichlet_dofs(sr.basis, problem.boundary_conditions())
    sys_pass, sys_rows = _linear_system_invariants(
        sr.basis, K, b, sr.solution, sr.pin_dof, dirichlet_dofs
    )

    # Source verification is skip-eligible; it doesn't fail the overall verdict
    # unless explicitly ✗.
    overall = sub_pass and (src_pass is not False) and sys_pass

    # --- Render internals.md ------------------------------------------------
    md = "\n".join(
        [
            f"# Internals report — `{problem.name}`",
            "",
            f"- Mesh size: h = {mesh_size:.4g}",
            f"- Mesh: {mesh.t.shape[1]} elements, {mesh.p.shape[1]} nodes",
            f"- **Overall: {_tick(overall)}**",
            "",
            "Generated by `scripts/diagnose.py`. See `CONVENTIONS.md` § Layer 2.",
            "",
            f"## Technique 1 — Subdomain integrity  {_tick(sub_pass)}",
            "",
            "| tag/name | n_elem | area | fraction | κ | metadata? |",
            "|----------|--------|------|----------|---|-----------|",
            *sub_rows,
            "",
            f"## Technique 2 — Source verification  {_tick(src_pass)}",
            "",
            *src_rows,
            "",
            f"## Technique 3 — Linear-system invariants  {_tick(sys_pass)}",
            "",
            *sys_rows,
            "",
        ]
    )
    out_path = out_dir / "internals.md"
    out_path.write_text(md, encoding="utf-8")
    logger.info("wrote %s  (overall %s)", out_path, _tick(overall))
    return out_path


# ============================================================================
# CLI
# ============================================================================


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)-7s %(name)s :: %(message)s"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("problem", choices=sorted(PROBLEMS))
    parser.add_argument(
        "--mesh-size",
        type=float,
        default=None,
        help="Mesh size to inspect. Default: finest from problem.mesh_sizes().",
    )
    args = parser.parse_args()
    diagnose(args.problem, args.mesh_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
