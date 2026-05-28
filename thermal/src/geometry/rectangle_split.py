# ABOUTME: gmsh builder for the Problem-2 rectangle split at x=1/2. The
# interface line is a feature edge (ADR-0003) so every triangle lies wholly
# inside one subdomain; subdomains and boundaries are tagged as named gmsh
# physical groups so meshio propagates them to scikit-fem as keys in
# mesh.subdomains and mesh.boundaries. Builder is cache-backed (ADR-0007).

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import gmsh
import numpy as np

from src.geometry.cache import cached_mesh

logger = logging.getLogger(__name__)


# ============================================================================
# CONSTANTS
# ============================================================================

GEOMETRY_NAME = "rectangle_split"

# Domain dimensions per verification.md § Problem 2: [0,1] x [0, SLAB_HEIGHT].
SLAB_HEIGHT = 0.1
_X_INTERFACE = 0.5

# Physical-group names (strings - meshio/skfem keys subdomains by name).
LEFT_SUBDOMAIN_NAME = "left_bulk"
RIGHT_SUBDOMAIN_NAME = "right_bulk"
LEFT_EDGE_NAME = "left"
RIGHT_EDGE_NAME = "right"
BOTTOM_EDGE_NAME = "bottom"
TOP_EDGE_NAME = "top"

# Numeric tags (gmsh requires them; meshio carries them in field_data but
# scikit-fem keys on the *name*, not the integer).
_LEFT_SUBDOMAIN_TAG = 1
_RIGHT_SUBDOMAIN_TAG = 2
_LEFT_EDGE_TAG = 10
_RIGHT_EDGE_TAG = 20
_BOTTOM_EDGE_TAG = 30
_TOP_EDGE_TAG = 40


# ============================================================================
# DATA
# ============================================================================


@dataclass(frozen=True)
class RectangleSplitSpec:
    """Parameter handle the solver uses to materialise the cached `.msh`.

    Only ``mesh_size`` and the slab height enter the cache key; kappa is not
    geometry and must not pollute the hash (acceptance #6 in submission 0003
    requires cache hits across the kappa_2 sweep).
    """

    mesh_size: float
    height: float = SLAB_HEIGHT
    geometry_name: str = GEOMETRY_NAME

    def params(self) -> dict[str, object]:
        return {"mesh_size": float(self.mesh_size), "height": float(self.height)}


# ============================================================================
# FUNCTIONS
# ============================================================================


def _build_msh(mesh_size: float, height: float, out_path: Path) -> None:
    """Run gmsh to produce a `.msh` with the x=1/2 interface as a feature edge.

    Six corner points (two domain bottoms, two interface, two domain tops);
    the interface edge p_b_mid -> p_t_mid bounds *both* surfaces, which is
    what forces every element to respect the discontinuity (ADR-0003).

    The two surfaces are meshed *transfinitely* (structured triangulation):
    each side of the slab is a uniform grid of right triangles. This costs
    nothing on a rectangle and removes the unstructured-mesh-quality noise
    that otherwise makes the asymptotic L^2 convergence rate fluctuate
    outside the planner's "within 0.2 of 2.0" acceptance window for an
    exact solution that is piecewise polynomial of degree <= 2. Structured
    meshing on a rectangular domain is not a physics choice and does not
    deviate from any ADR; the alignment requirement (ADR-0003) is the only
    binding mesh constraint and it is satisfied trivially.
    """
    # Number of node points along each transfinite curve. gmsh wants the
    # *node* count (segments + 1); round up so the effective edge length
    # stays at or below ``mesh_size``.
    half_width = _X_INTERFACE
    n_along_x_half = max(2, int(np.ceil(half_width / mesh_size)) + 1)
    n_along_y = max(2, int(np.ceil(height / mesh_size)) + 1)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Corners. _b = bottom row (y=0), _t = top row (y=height);
        # _l = x=0, _mid = x=0.5, _r = x=1. mesh_size is supplied as a
        # default characteristic length, but transfinite overrides it.
        p_b_l = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, mesh_size)
        p_b_mid = gmsh.model.geo.addPoint(_X_INTERFACE, 0.0, 0.0, mesh_size)
        p_b_r = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, mesh_size)
        p_t_r = gmsh.model.geo.addPoint(1.0, height, 0.0, mesh_size)
        p_t_mid = gmsh.model.geo.addPoint(_X_INTERFACE, height, 0.0, mesh_size)
        p_t_l = gmsh.model.geo.addPoint(0.0, height, 0.0, mesh_size)

        # Edges. The interface edge is shared between the two surfaces and is
        # the ADR-0003 feature edge that pins triangles to one side or the
        # other.
        e_b_left = gmsh.model.geo.addLine(p_b_l, p_b_mid)
        e_b_right = gmsh.model.geo.addLine(p_b_mid, p_b_r)
        e_right = gmsh.model.geo.addLine(p_b_r, p_t_r)
        e_t_right = gmsh.model.geo.addLine(p_t_r, p_t_mid)
        e_t_left = gmsh.model.geo.addLine(p_t_mid, p_t_l)
        e_left = gmsh.model.geo.addLine(p_t_l, p_b_l)
        e_interface = gmsh.model.geo.addLine(p_b_mid, p_t_mid)

        # Left surface: counterclockwise around the left subdomain.
        loop_left = gmsh.model.geo.addCurveLoop(
            [e_b_left, e_interface, e_t_left, e_left]
        )
        # Right surface: the interface edge is traversed in the *opposite*
        # direction (signed -e_interface) to keep the right loop CCW.
        loop_right = gmsh.model.geo.addCurveLoop(
            [e_b_right, e_right, e_t_right, -e_interface]
        )
        surf_left = gmsh.model.geo.addPlaneSurface([loop_left])
        surf_right = gmsh.model.geo.addPlaneSurface([loop_right])

        # Transfinite (structured) meshing. The two horizontal segments of
        # each subdomain get matching node counts so the interface gets
        # conforming triangles on both sides.
        for line_id in (e_b_left, e_t_left):
            gmsh.model.geo.mesh.setTransfiniteCurve(line_id, n_along_x_half)
        for line_id in (e_b_right, e_t_right):
            gmsh.model.geo.mesh.setTransfiniteCurve(line_id, n_along_x_half)
        for line_id in (e_left, e_interface, e_right):
            gmsh.model.geo.mesh.setTransfiniteCurve(line_id, n_along_y)
        gmsh.model.geo.mesh.setTransfiniteSurface(surf_left)
        gmsh.model.geo.mesh.setTransfiniteSurface(surf_right)

        gmsh.model.geo.synchronize()

        # Physical groups - names are the keys scikit-fem will see in
        # mesh.subdomains / mesh.boundaries.
        gmsh.model.addPhysicalGroup(
            2, [surf_left], tag=_LEFT_SUBDOMAIN_TAG, name=LEFT_SUBDOMAIN_NAME
        )
        gmsh.model.addPhysicalGroup(
            2, [surf_right], tag=_RIGHT_SUBDOMAIN_TAG, name=RIGHT_SUBDOMAIN_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [e_left], tag=_LEFT_EDGE_TAG, name=LEFT_EDGE_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [e_right], tag=_RIGHT_EDGE_TAG, name=RIGHT_EDGE_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [e_b_left, e_b_right], tag=_BOTTOM_EDGE_TAG, name=BOTTOM_EDGE_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [e_t_left, e_t_right], tag=_TOP_EDGE_TAG, name=TOP_EDGE_NAME
        )

        gmsh.model.mesh.generate(2)
        gmsh.write(str(out_path))
    finally:
        gmsh.finalize()


def build_rectangle_split(mesh_size: float) -> RectangleSplitSpec:
    """Public builder callable returned by ``Problem02Slab.geometry()``."""
    return RectangleSplitSpec(mesh_size=mesh_size)


def materialise(spec: RectangleSplitSpec, cache_dir: Path) -> Path:
    """Return the path to a (possibly cached) `.msh` for ``spec``."""
    return cached_mesh(
        cache_dir=cache_dir,
        geometry_name=spec.geometry_name,
        params=spec.params(),
        build=lambda out: _build_msh(spec.mesh_size, spec.height, out),
    )
