# ABOUTME: gmsh OCC-kernel builder for Problem 4: two inner disks of radius
# R_inner centred at (+/- d_sep/2, 0) embedded in a larger disk of radius
# R_outer. Three addDisk calls + occ.fragment([(2, outer)], [(2, innerA),
# (2, innerB)]) give conforming shared edges at both r=R_inner circles
# (ADR-0003). Surface and curve tags are recovered from the fragment map and
# from per-surface getBoundary lookups — never hard-coded — because fragment
# renumbers entities. synchronize() is called after fragment and before any
# getBoundary or addPhysicalGroup call (OCC sequencing rule). Two distinct
# inner subdomains "inner_disk_A" and "inner_disk_B" are emitted rather than
# collapsed (submission 0008 § Decisions 6), so an asymmetric-kappa variant
# slots in without re-tagging and so the mirror-symmetry forced-failure check
# stays sharp. Cache key: {R_inner, d_sep, R_outer, mesh_size} (submission
# 0008 § Decisions 7); kappa is not geometry and must not enter the key.

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import gmsh

from src.geometry.cache import cached_mesh

logger = logging.getLogger(__name__)


# ============================================================================
# CONSTANTS
# ============================================================================

GEOMETRY_NAME = "two_disks_in_disk"

# Physical-group names: strings are the keys scikit-fem sees in
# mesh.subdomains / mesh.boundaries.
INNER_DISK_A_NAME = "inner_disk_A"
INNER_DISK_B_NAME = "inner_disk_B"
OUTER_ANNULUS_NAME = "outer_annulus"
INNER_BOUNDARY_A_NAME = "inner_boundary_A"  # tagged for inspection; no BC attached
INNER_BOUNDARY_B_NAME = "inner_boundary_B"  # tagged for inspection; no BC attached
OUTER_BOUNDARY_NAME = "outer_boundary"

# Numeric tags (gmsh requires them; scikit-fem keys on the *name*, not the int).
_INNER_DISK_A_TAG = 1
_INNER_DISK_B_TAG = 2
_OUTER_ANNULUS_TAG = 3
_INNER_BOUNDARY_A_TAG = 10
_INNER_BOUNDARY_B_TAG = 11
_OUTER_BOUNDARY_TAG = 20


# ============================================================================
# DATA
# ============================================================================


@dataclass(frozen=True)
class TwoDisksInDiskSpec:
    """Parameter handle the solver uses to materialise the cached `.msh`.

    Only R_inner, d_sep, R_outer, and mesh_size enter the cache key (submission
    0008 § Decisions 7). Disk centres at (+/- d_sep/2, 0) are determined by
    d_sep so they need not be stored separately. kappa values are not geometry
    and must not pollute the hash.
    """

    R_inner: float
    d_sep: float
    R_outer: float
    mesh_size: float
    geometry_name: str = GEOMETRY_NAME

    def params(self) -> dict[str, object]:
        return {
            "R_inner": float(self.R_inner),
            "d_sep": float(self.d_sep),
            "R_outer": float(self.R_outer),
            "mesh_size": float(self.mesh_size),
        }


# ============================================================================
# FUNCTIONS
# ============================================================================


def _build_msh(
    R_inner: float,
    d_sep: float,
    R_outer: float,
    mesh_size: float,
    out_path: Path,
) -> None:
    """Run gmsh (OCC kernel) to produce a `.msh` with both r=R_inner circles
    as conforming shared edges (ADR-0003).

    addDisk(outer) + addDisk(innerA) + addDisk(innerB) are joined by a single
    occ.fragment call with the outer disk as the object and both inner disks
    as tools; post-boolean surface and curve tags are recovered from the
    outDimTagsMap so they remain correct even if fragment renumbers them.

    OCC sequencing: synchronize() MUST follow fragment and precede any
    getBoundary, getBoundingBox, or addPhysicalGroup call.
    """
    a = 0.5 * d_sep
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Three disks in the OCC kernel. Tags may change after fragment.
        outer_tag_init = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, R_outer, R_outer)
        inner_a_tag_init = gmsh.model.occ.addDisk(-a, 0.0, 0.0, R_inner, R_inner)
        inner_b_tag_init = gmsh.model.occ.addDisk(+a, 0.0, 0.0, R_inner, R_inner)

        # fragment(object, tools) creates conforming shared boundaries at both
        # inner circles. With one outer and two inner inputs, fragment returns
        # three fragments for the outer (annulus + the two inner disks) and
        # exactly one for each inner tool. Tag recovery uses out_map (not
        # centroid classification, ADR-0003 anti-pattern).
        _out_tags, out_map = gmsh.model.occ.fragment(
            [(2, outer_tag_init)],
            [(2, inner_a_tag_init), (2, inner_b_tag_init)],
        )

        # synchronize() BEFORE any getBoundary or addPhysicalGroup (OCC rule).
        gmsh.model.occ.synchronize()

        if len(out_map[1]) != 1:
            raise RuntimeError(
                f"Expected inner disk A to produce 1 fragment; got "
                f"{len(out_map[1])}. Fragment result: out_map[1]={out_map[1]}"
            )
        if len(out_map[2]) != 1:
            raise RuntimeError(
                f"Expected inner disk B to produce 1 fragment; got "
                f"{len(out_map[2])}. Fragment result: out_map[2]={out_map[2]}"
            )
        if len(out_map[0]) != 3:
            raise RuntimeError(
                f"Expected outer disk to produce 3 fragments (annulus + 2 "
                f"inner disks); got {len(out_map[0])}. Fragment result: "
                f"out_map[0]={out_map[0]}"
            )

        inner_a_surf = out_map[1][0][1]
        inner_b_surf = out_map[2][0][1]
        annulus_surf = next(
            t for _, t in out_map[0] if t not in (inner_a_surf, inner_b_surf)
        )

        # Identify boundary curves. Each inner disk has exactly one boundary
        # curve (its own circle). The annulus has three: the outer circle and
        # both shared inner circles. Use set membership to pick out each name.
        inner_a_curves = [
            abs(t)
            for _, t in gmsh.model.getBoundary([(2, inner_a_surf)], oriented=False)
        ]
        inner_b_curves = [
            abs(t)
            for _, t in gmsh.model.getBoundary([(2, inner_b_surf)], oriented=False)
        ]
        annulus_curves = [
            abs(t)
            for _, t in gmsh.model.getBoundary([(2, annulus_surf)], oriented=False)
        ]

        if len(inner_a_curves) != 1:
            raise RuntimeError(
                f"Inner disk A should have exactly 1 boundary curve; got "
                f"{inner_a_curves}"
            )
        if len(inner_b_curves) != 1:
            raise RuntimeError(
                f"Inner disk B should have exactly 1 boundary curve; got "
                f"{inner_b_curves}"
            )

        inner_a_circle = inner_a_curves[0]
        inner_b_circle = inner_b_curves[0]
        outer_circle = next(
            t for t in annulus_curves if t not in (inner_a_circle, inner_b_circle)
        )

        # Uniform target mesh size on all OCC points. gmsh refines near
        # curved boundaries by default; mesh_size is the target, not a cap.
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

        # Physical groups — string names are the keys scikit-fem sees.
        gmsh.model.addPhysicalGroup(
            2, [inner_a_surf], tag=_INNER_DISK_A_TAG, name=INNER_DISK_A_NAME
        )
        gmsh.model.addPhysicalGroup(
            2, [inner_b_surf], tag=_INNER_DISK_B_TAG, name=INNER_DISK_B_NAME
        )
        gmsh.model.addPhysicalGroup(
            2, [annulus_surf], tag=_OUTER_ANNULUS_TAG, name=OUTER_ANNULUS_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [inner_a_circle], tag=_INNER_BOUNDARY_A_TAG, name=INNER_BOUNDARY_A_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [inner_b_circle], tag=_INNER_BOUNDARY_B_TAG, name=INNER_BOUNDARY_B_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [outer_circle], tag=_OUTER_BOUNDARY_TAG, name=OUTER_BOUNDARY_NAME
        )

        gmsh.model.mesh.generate(2)
        gmsh.write(str(out_path))
        logger.debug(
            "two_disks_in_disk mesh written: R_inner=%.3f d_sep=%.3f "
            "R_outer=%.3f h=%.4f -> %s",
            R_inner,
            d_sep,
            R_outer,
            mesh_size,
            out_path,
        )
    finally:
        gmsh.finalize()


def build_two_disks_in_disk(
    R_inner: float, d_sep: float, R_outer: float
) -> Callable[[float], TwoDisksInDiskSpec]:
    """Return the geometry callable stored in Problem04TwoDisks.geometry().

    The returned callable maps mesh_size -> TwoDisksInDiskSpec; it is what
    the geometry-cache layer calls to produce the spec passed to materialise().
    """

    def _builder(mesh_size: float) -> TwoDisksInDiskSpec:
        return TwoDisksInDiskSpec(
            R_inner=R_inner,
            d_sep=d_sep,
            R_outer=R_outer,
            mesh_size=mesh_size,
        )

    return _builder


def materialise(spec: TwoDisksInDiskSpec, cache_dir: Path) -> Path:
    """Return the path to a (possibly cached) `.msh` for ``spec``."""
    return cached_mesh(
        cache_dir=cache_dir,
        geometry_name=spec.geometry_name,
        params=spec.params(),
        build=lambda out: _build_msh(
            spec.R_inner, spec.d_sep, spec.R_outer, spec.mesh_size, out
        ),
    )
