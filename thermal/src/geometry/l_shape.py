# ABOUTME: gmsh OCC-kernel builder for the Problem-5 L-shape geometry. The
# domain is [0,1]^2 with the upper-right quadrant removed; the reentrant corner
# sits at (1/2, 1/2). Six explicit OCC points -> six lines -> curve loop ->
# plane surface guarantees the reentrant corner is a mesh node at every
# refinement (verification submission 0007 § Acceptance 6). One named subdomain
# ("interior") and six named boundary edges; the two cut edges meeting at the
# reentrant corner are tagged distinctly so callers can apply T=0 there while
# the four outer edges carry the singular Dirichlet callable (submission 0007
# § Decisions 4-5). Builder is cache-backed (ADR-0007).

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

GEOMETRY_NAME = "l_shape"

# Physical-group names: strings are the keys scikit-fem sees in
# mesh.subdomains / mesh.boundaries.
INTERIOR_NAME = "interior"

# Boundary names. Cut edges meet at the reentrant corner (1/2, 1/2); on these
# the singular exact solution T = r^{2/3} sin(2 theta/3) is identically zero,
# so Problem 5 uses scalar 0.0 there. The other four edges form the outer
# perimeter of the L and carry the callable form of the BC.
CUT_EAST_NAME = "cut_east"  # y=1/2, x in [1/2, 1]; theta=0 in polar convention
CUT_NORTH_NAME = "cut_north"  # x=1/2, y in [1/2, 1]; theta=3pi/2 in polar convention
SOUTH_NAME = "south"  # y=0, x in [0, 1]
EAST_LOWER_NAME = "east_lower"  # x=1, y in [0, 1/2]
NORTH_LEFT_NAME = "north_left"  # y=1, x in [0, 1/2]
WEST_NAME = "west"  # x=0, y in [0, 1]

# Geometry constants (the L-shape is fixed; only mesh_size varies).
_X_CUT = 0.5
_Y_CUT = 0.5

# Numeric tags (gmsh requires them; scikit-fem keys on the *name*).
_INTERIOR_TAG = 1
_SOUTH_TAG = 11
_EAST_LOWER_TAG = 12
_CUT_EAST_TAG = 13
_CUT_NORTH_TAG = 14
_NORTH_LEFT_TAG = 15
_WEST_TAG = 16


# ============================================================================
# DATA
# ============================================================================


@dataclass(frozen=True)
class LShapeSpec:
    """Parameter handle the solver uses to materialise the cached `.msh`.

    The L-shape is geometry-fixed (one subdomain, six edges, reentrant corner
    at (1/2, 1/2)); only mesh_size enters the cache key.
    """

    mesh_size: float
    geometry_name: str = GEOMETRY_NAME

    def params(self) -> dict[str, object]:
        return {"mesh_size": float(self.mesh_size)}


# ============================================================================
# FUNCTIONS
# ============================================================================


def _build_msh(mesh_size: float, out_path: Path) -> None:
    """Build the L-shape `.msh` via the OCC kernel.

    Six OCC points -> six OCC lines -> curve loop -> plane surface, traversed
    counterclockwise so the interior is on the left. The reentrant corner
    (1/2, 1/2) is an explicit OCC point and therefore guaranteed a mesh node
    at every refinement (submission 0007 § Acceptance 6).
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(GEOMETRY_NAME)

        # Six L-shape vertices (counterclockwise: interior on the left).
        # Reentrant corner is p_cut at (1/2, 1/2).
        p_sw = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
        p_se = gmsh.model.occ.addPoint(1.0, 0.0, 0.0)
        p_e_mid = gmsh.model.occ.addPoint(1.0, _Y_CUT, 0.0)
        p_cut = gmsh.model.occ.addPoint(_X_CUT, _Y_CUT, 0.0)
        p_n_mid = gmsh.model.occ.addPoint(_X_CUT, 1.0, 0.0)
        p_nw = gmsh.model.occ.addPoint(0.0, 1.0, 0.0)

        l_south = gmsh.model.occ.addLine(p_sw, p_se)
        l_east_lower = gmsh.model.occ.addLine(p_se, p_e_mid)
        l_cut_east = gmsh.model.occ.addLine(p_e_mid, p_cut)
        l_cut_north = gmsh.model.occ.addLine(p_cut, p_n_mid)
        l_north_left = gmsh.model.occ.addLine(p_n_mid, p_nw)
        l_west = gmsh.model.occ.addLine(p_nw, p_sw)

        loop = gmsh.model.occ.addCurveLoop(
            [l_south, l_east_lower, l_cut_east, l_cut_north, l_north_left, l_west]
        )
        surf = gmsh.model.occ.addPlaneSurface([loop])

        # synchronize() before any mesh / physical-group call (OCC rule).
        gmsh.model.occ.synchronize()

        # Uniform target mesh size on every OCC point. No grading toward the
        # reentrant corner: submission 0007 § Out of scope binds (graded meshes
        # would recover rate 2 and defeat the inverted assertion).
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

        # Physical groups - string names are the keys scikit-fem sees.
        gmsh.model.addPhysicalGroup(2, [surf], tag=_INTERIOR_TAG, name=INTERIOR_NAME)
        gmsh.model.addPhysicalGroup(
            1, [l_south], tag=_SOUTH_TAG, name=SOUTH_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [l_east_lower], tag=_EAST_LOWER_TAG, name=EAST_LOWER_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [l_cut_east], tag=_CUT_EAST_TAG, name=CUT_EAST_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [l_cut_north], tag=_CUT_NORTH_TAG, name=CUT_NORTH_NAME
        )
        gmsh.model.addPhysicalGroup(
            1, [l_north_left], tag=_NORTH_LEFT_TAG, name=NORTH_LEFT_NAME
        )
        gmsh.model.addPhysicalGroup(1, [l_west], tag=_WEST_TAG, name=WEST_NAME)

        gmsh.model.mesh.generate(2)
        gmsh.write(str(out_path))
        logger.debug(
            "l_shape mesh written: h=%.4f -> %s",
            mesh_size,
            out_path,
        )
    finally:
        gmsh.finalize()


def build_l_shape(mesh_size: float) -> LShapeSpec:
    """Public builder callable returned by ``Problem05LShape.geometry()``."""
    return LShapeSpec(mesh_size=mesh_size)


def materialise(spec: LShapeSpec, cache_dir: Path) -> Path:
    """Return the path to a (possibly cached) `.msh` for ``spec``."""
    return cached_mesh(
        cache_dir=cache_dir,
        geometry_name=spec.geometry_name,
        params=spec.params(),
        build=lambda out: _build_msh(spec.mesh_size, out),
    )
