# ABOUTME: SHA-256-keyed on-disk cache for generated `.msh` files (ADR-0007).
# Builders register a deterministic key (geometry name + parameter dict); the
# cache returns the cached path on a hit, otherwise invokes the builder and
# stores the result. Re-meshing is the dominant cost of refinement studies, so
# this layer must exist even when individual mesh builds (e.g. unit square)
# are cheap.

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Callable, Mapping

logger = logging.getLogger(__name__)


# ============================================================================
# FUNCTIONS
# ============================================================================


def cache_key(geometry_name: str, params: Mapping[str, object]) -> str:
    """Deterministic SHA-256 key for a (geometry_name, params) pair.

    Parameters are serialised with sorted keys so dict ordering doesn't perturb
    the hash. Floats are stored with full ``repr`` precision.
    """
    payload = json.dumps(
        {"geometry": geometry_name, "params": dict(params)},
        sort_keys=True,
        default=repr,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def cached_mesh(
    cache_dir: Path,
    geometry_name: str,
    params: Mapping[str, object],
    build: Callable[[Path], None],
) -> Path:
    """Return the path to the cached ``.msh`` file, building it if absent.

    The builder is invoked with the target output path. The cache is hit when
    the file already exists; mtimes are not updated on hit (acceptance #6 in
    the submission brief).
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = cache_key(geometry_name, params)
    out = cache_dir / f"{geometry_name}_{key[:16]}.msh"
    if out.exists():
        logger.debug("mesh cache hit: %s", out)
        return out
    logger.debug("mesh cache miss; building: %s", out)
    build(out)
    if not out.exists():
        raise RuntimeError(f"mesh builder for {geometry_name} did not produce {out}")
    return out
