# ABOUTME: gmsh builder for the unit-square geometry used by Problem 1. The
# builder is cache-backed (ADR-0007) and returns a parameter handle the solver
# can consume to materialise a scikit-fem MeshTri (subdomain tag 1, boundaries
# tagged 1..4 for x=0, x=1, y=0, y=1).

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import gmsh

from src.geometry.cache import cached_mesh

logger = logging.getLogger(__name__)


# ============================================================================
# CONSTANTS
# ============================================================================

GEOMETRY_NAME = "unit_square"

# Physical group tags. Surface tag = 1 (single subdomain). Boundary tags
# numbered counterclockwise from the bottom edge.
SUBDOMAIN_TAG = 1
BOUNDARY_TAGS = {
    "bottom": 1,  # y = 0
    "right": 2,  # x = 1
    "top": 3,  # y = 1
    "left": 4,  # x = 0
}


# ============================================================================
# DATA
# ============================================================================


@dataclass(frozen=True)
class UnitSquareSpec:
    """Parameter handle a solver consumes to locate the cached `.msh`.

    Held as data, not as a gmsh object, so the solver can decide when to
    materialise the mesh. ``mesh_size`` is the gmsh characteristic length used
    at every corner.
    """

    mesh_size: float
    geometry_name: str = GEOMETRY_NAME

    def params(self) -> dict[str, object]:
        return {"mesh_size": float(self.mesh_size)}


# ============================================================================
# FUNCTIONS
# ============================================================================


def _build_msh(mesh_size: float, out_path: Path) -> None:
    """Run gmsh to produce a `.msh` file for the unit square at ``out_path``."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Four corners with the requested characteristic length.
        p0 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, mesh_size)
        p1 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, mesh_size)
        p2 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, mesh_size)
        p3 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, mesh_size)

        # Four edges; orientation matches BOUNDARY_TAGS keys.
        l_bottom = gmsh.model.geo.addLine(p0, p1)
        l_right = gmsh.model.geo.addLine(p1, p2)
        l_top = gmsh.model.geo.addLine(p2, p3)
        l_left = gmsh.model.geo.addLine(p3, p0)

        loop = gmsh.model.geo.addCurveLoop([l_bottom, l_right, l_top, l_left])
        surf = gmsh.model.geo.addPlaneSurface([loop])

        gmsh.model.geo.synchronize()

        # Tag the surface and each boundary edge as physical groups so meshio
        # propagates them through to scikit-fem.
        gmsh.model.addPhysicalGroup(2, [surf], tag=SUBDOMAIN_TAG, name="bulk")
        gmsh.model.addPhysicalGroup(
            1, [l_bottom], tag=BOUNDARY_TAGS["bottom"], name="bottom"
        )
        gmsh.model.addPhysicalGroup(
            1, [l_right], tag=BOUNDARY_TAGS["right"], name="right"
        )
        gmsh.model.addPhysicalGroup(1, [l_top], tag=BOUNDARY_TAGS["top"], name="top")
        gmsh.model.addPhysicalGroup(1, [l_left], tag=BOUNDARY_TAGS["left"], name="left")

        gmsh.model.mesh.generate(2)
        # gmsh's default v4 ASCII; meshio reads it via from_meshio downstream.
        gmsh.write(str(out_path))
    finally:
        gmsh.finalize()


def build_unit_square(mesh_size: float) -> UnitSquareSpec:
    """Public builder callable returned by ``Problem.geometry()``.

    Returns a lightweight parameter handle; the solver materialises the
    cached `.msh` file via :func:`materialise`.
    """
    return UnitSquareSpec(mesh_size=mesh_size)


def materialise(spec: UnitSquareSpec, cache_dir: Path) -> Path:
    """Return the path to a (possibly cached) `.msh` for ``spec``."""
    return cached_mesh(
        cache_dir=cache_dir,
        geometry_name=spec.geometry_name,
        params=spec.params(),
        build=lambda out: _build_msh(spec.mesh_size, out),
    )
