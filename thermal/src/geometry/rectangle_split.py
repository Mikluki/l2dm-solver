# ABOUTME: gmsh OCC-kernel builder for the Problem-2 rectangle split at x=1/2.
# Two addRectangle calls are joined by occ.fragment so the interface at x=1/2
# is a conforming shared edge (ADR-0003); post-fragment surface and curve tags
# are recovered via bounding-box queries — never hard-coded through the boolean.
# Transfinite (structured) meshing is applied after occ.synchronize() to keep
# convergence rate inside the acceptance window. Builder is cache-backed (ADR-0007).

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


def _curve_bbox(curve_tag: int) -> tuple[float, float, float, float]:
    """Return (xmin, ymin, xmax, ymax) for a 1-D entity."""
    xmin, ymin, _zmin, xmax, ymax, _zmax = gmsh.model.getBoundingBox(1, curve_tag)
    return xmin, ymin, xmax, ymax


def _build_msh(mesh_size: float, height: float, out_path: Path) -> None:
    """Run gmsh (OCC kernel) to produce a `.msh` with the x=1/2 interface as a
    conforming shared edge (ADR-0003).

    Two addRectangle halves are joined by occ.fragment; post-boolean surface
    and curve tags are recovered via bounding-box queries so they remain correct
    even if fragment renumbers them. Transfinite (structured) meshing is then
    applied to both surfaces to keep the L^2 convergence rate inside the
    acceptance window (same rationale as the original geo-kernel version).
    """
    # Node count along each transfinite curve: segments + 1, rounded up so the
    # effective edge length stays at or below mesh_size.
    half_width = _X_INTERFACE
    n_along_x_half = max(2, int(np.ceil(half_width / mesh_size)) + 1)
    n_along_y = max(2, int(np.ceil(height / mesh_size)) + 1)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Two axis-aligned rectangles sharing the edge at x = _X_INTERFACE.
        surf_left_init = gmsh.model.occ.addRectangle(
            0.0, 0.0, 0.0, _X_INTERFACE, height
        )
        surf_right_init = gmsh.model.occ.addRectangle(
            _X_INTERFACE, 0.0, 0.0, 1.0 - _X_INTERFACE, height
        )

        # fragment merges the shared edge so it belongs to both surfaces
        # (ADR-0003). outDimTagsMap[0] = surfaces from surf_left_init;
        # outDimTagsMap[1] = surfaces from surf_right_init.
        _out, out_map = gmsh.model.occ.fragment(
            [(2, surf_left_init)], [(2, surf_right_init)]
        )
        surf_left = out_map[0][0][1]
        surf_right = out_map[1][0][1]

        # synchronize before any mesh or physical-group call (OCC requirement).
        gmsh.model.occ.synchronize()

        # Identify boundary curves by bounding box. The two rectangles have
        # axis-aligned edges, so xmin==xmax identifies a vertical line and
        # ymin==ymax a horizontal one.
        EPS = 1e-9
        left_curves = [
            t for _, t in gmsh.model.getBoundary([(2, surf_left)], oriented=False)
        ]
        right_curves = [
            t for _, t in gmsh.model.getBoundary([(2, surf_right)], oriented=False)
        ]

        def _vert_at(tag: int, x: float) -> bool:
            xmin, _ym, xmax, _yM = _curve_bbox(tag)
            return abs(xmin - x) < EPS and abs(xmax - x) < EPS

        def _horiz_at(tag: int, y: float) -> bool:
            _xm, ymin, _xM, ymax = _curve_bbox(tag)
            return abs(ymin - y) < EPS and abs(ymax - y) < EPS

        e_left = next(t for t in left_curves if _vert_at(t, 0.0))
        e_interface = next(t for t in left_curves if _vert_at(t, _X_INTERFACE))
        e_b_left = next(t for t in left_curves if _horiz_at(t, 0.0))
        e_t_left = next(t for t in left_curves if _horiz_at(t, height))

        e_right = next(t for t in right_curves if _vert_at(t, 1.0))
        e_b_right = next(t for t in right_curves if _horiz_at(t, 0.0))
        e_t_right = next(t for t in right_curves if _horiz_at(t, height))

        # Transfinite (structured) meshing — kernel-agnostic API after synchronize.
        for line_id in (e_b_left, e_t_left, e_b_right, e_t_right):
            gmsh.model.mesh.setTransfiniteCurve(line_id, n_along_x_half)
        for line_id in (e_left, e_interface, e_right):
            gmsh.model.mesh.setTransfiniteCurve(line_id, n_along_y)
        gmsh.model.mesh.setTransfiniteSurface(surf_left)
        gmsh.model.mesh.setTransfiniteSurface(surf_right)

        # Physical groups — string names are the keys scikit-fem sees in
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
