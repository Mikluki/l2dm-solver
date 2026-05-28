# ABOUTME: gmsh OCC-kernel builder for Problem 3: inner disk of radius R_inner
# embedded in a larger disk of radius R_outer. Uses addDisk × 2 + occ.fragment
# to produce a conforming shared edge at r=R_inner (ADR-0003). Post-fragment
# surface tags are recovered from the outDimTagsMap returned by fragment — never
# hard-coded — because fragment renumbers entities. Synchronize() is called
# after fragment and before any getBoundary or addPhysicalGroup call (OCC
# sequencing rule). No transfinite meshing: unstructured Delaunay (brief §
# Decisions 3). Cache key: {R_inner, R_outer, mesh_size}; kappa is not geometry
# and must not enter the key (brief § Decisions 5). Builder is cache-backed
# (ADR-0007).

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

GEOMETRY_NAME = "disk_in_disk"

# Physical-group names: strings are the keys scikit-fem sees in
# mesh.subdomains / mesh.boundaries.
INNER_DISK_NAME = "inner_disk"
OUTER_ANNULUS_NAME = "outer_annulus"
INNER_BOUNDARY_NAME = "inner_boundary"  # tagged for inspection; no BC attached
OUTER_BOUNDARY_NAME = "outer_boundary"

# Numeric tags (gmsh requires them; scikit-fem keys on the *name*, not the int).
_INNER_DISK_TAG = 1
_OUTER_ANNULUS_TAG = 2
_INNER_BOUNDARY_TAG = 10
_OUTER_BOUNDARY_TAG = 20


# ============================================================================
# DATA
# ============================================================================


@dataclass(frozen=True)
class DiskInDiskSpec:
    """Parameter handle the solver uses to materialise the cached `.msh`.

    Only R_inner, R_outer, and mesh_size enter the cache key; kappa values are
    not geometry and must not pollute the hash (brief § Decisions 5 requires
    cache hits across the kappa_1 sweep).
    """

    R_inner: float
    R_outer: float
    mesh_size: float
    geometry_name: str = GEOMETRY_NAME

    def params(self) -> dict[str, object]:
        return {
            "R_inner": float(self.R_inner),
            "R_outer": float(self.R_outer),
            "mesh_size": float(self.mesh_size),
        }


# ============================================================================
# FUNCTIONS
# ============================================================================


def _build_msh(
    R_inner: float, R_outer: float, mesh_size: float, out_path: Path
) -> None:
    """Run gmsh (OCC kernel) to produce a `.msh` with the r=R_inner circle as
    a conforming shared edge (ADR-0003).

    addDisk(outer) + addDisk(inner) are joined by occ.fragment; post-boolean
    surface and curve tags are recovered from the outDimTagsMap so they remain
    correct even if fragment renumbers them. No transfinite meshing — brief §
    Decisions 3 forbids it for disk geometries.

    OCC sequencing: synchronize() MUST follow fragment and precede any
    getBoundary, getBoundingBox, or addPhysicalGroup call.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Create both disks in the OCC kernel. addDisk(xc, yc, zc, rx, ry)
        # with rx=ry creates a circle. Tags may change after fragment.
        outer_tag_init = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, R_outer, R_outer)
        inner_tag_init = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, R_inner, R_inner)

        # fragment(object, tool) creates a conforming shared boundary at
        # r=R_inner (ADR-0003). outDimTagsMap[0] = new surfaces that replace
        # outer_tag_init (split into annulus + shared inner disk);
        # outDimTagsMap[1] = new surfaces that replace inner_tag_init (just
        # the shared inner disk, now the canonical copy).
        _out_tags, out_map = gmsh.model.occ.fragment(
            [(2, outer_tag_init)], [(2, inner_tag_init)]
        )

        # synchronize() BEFORE any getBoundary or addPhysicalGroup (OCC rule).
        gmsh.model.occ.synchronize()

        # Identify surfaces: inner disk is the single fragment of the inner
        # input; annulus is the outer-input fragment that is not the inner disk.
        if len(out_map[1]) != 1:
            raise RuntimeError(
                f"Expected inner disk to produce 1 fragment; got {len(out_map[1])}. "
                f"Fragment result: out_map[1]={out_map[1]}"
            )
        if len(out_map[0]) != 2:
            raise RuntimeError(
                f"Expected outer disk to produce 2 fragments (annulus + inner); "
                f"got {len(out_map[0])}. Fragment result: out_map[0]={out_map[0]}"
            )

        inner_surf = out_map[1][0][1]
        annulus_surf = next(t for _, t in out_map[0] if t != inner_surf)

        # Identify boundary curves. The inner disk has exactly one boundary
        # curve (the inner circle). The annulus has two: the shared inner
        # circle and the outer circle. Use set difference to find each.
        inner_curves = [
            abs(t) for _, t in gmsh.model.getBoundary([(2, inner_surf)], oriented=False)
        ]
        annulus_curves = [
            abs(t) for _, t in gmsh.model.getBoundary([(2, annulus_surf)], oriented=False)
        ]

        if len(inner_curves) != 1:
            raise RuntimeError(
                f"Inner disk should have exactly 1 boundary curve; got {inner_curves}"
            )

        inner_circle = inner_curves[0]
        outer_circle = next(t for t in annulus_curves if t != inner_circle)

        # Set mesh size on all entities: uniform target. gmsh will refine
        # near curved boundaries by default; mesh_size is the target, not a cap.
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

        # Physical groups — string names are the keys scikit-fem sees.
        gmsh.model.addPhysicalGroup(
            2, [inner_surf], tag=_INNER_DISK_TAG, name=INNER_DISK_NAME
        )
        gmsh.model.addPhysicalGroup(
            2, [annulus_surf], tag=_OUTER_ANNULUS_TAG, name=OUTER_ANNULUS_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [inner_circle], tag=_INNER_BOUNDARY_TAG, name=INNER_BOUNDARY_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [outer_circle], tag=_OUTER_BOUNDARY_TAG, name=OUTER_BOUNDARY_NAME
        )

        gmsh.model.mesh.generate(2)
        gmsh.write(str(out_path))
        logger.debug(
            "disk_in_disk mesh written: R_inner=%.3f R_outer=%.3f h=%.4f -> %s",
            R_inner,
            R_outer,
            mesh_size,
            out_path,
        )
    finally:
        gmsh.finalize()


def build_disk_in_disk(R_inner: float, R_outer: float) -> Callable[[float], DiskInDiskSpec]:
    """Return the geometry callable stored in Problem03Disk.geometry().

    The returned callable maps mesh_size -> DiskInDiskSpec; it is what the
    geometry-cache layer calls to produce the spec passed to materialise().
    """

    def _builder(mesh_size: float) -> DiskInDiskSpec:
        return DiskInDiskSpec(R_inner=R_inner, R_outer=R_outer, mesh_size=mesh_size)

    return _builder


def materialise(spec: DiskInDiskSpec, cache_dir: Path) -> Path:
    """Return the path to a (possibly cached) `.msh` for ``spec``."""
    return cached_mesh(
        cache_dir=cache_dir,
        geometry_name=spec.geometry_name,
        params=spec.params(),
        build=lambda out: _build_msh(spec.R_inner, spec.R_outer, spec.mesh_size, out),
    )
