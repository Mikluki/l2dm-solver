# ABOUTME: Pytest configuration: puts the project root on sys.path so the
# tests can import the `src.*` namespace package, and provides a per-test
# artifact-directory fixture rooted at tests/_artifacts/<test_id>/.

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest


# ============================================================================
# CONFIG
# ============================================================================

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

ARTIFACT_ROOT = _PROJECT_ROOT / "tests" / "_artifacts"
MESH_CACHE_ROOT = _PROJECT_ROOT / "tests" / "_mesh_cache"

# Silence skfem's per-assemble chatter at INFO. Our own modules (src.*) stay
# at INFO so the per-level error table and fitted rates remain visible.
# Raise to DEBUG via `--log-cli-level=DEBUG` when you need the skfem internals.
logging.getLogger("skfem").setLevel(logging.WARNING)


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def artifact_dir(request: pytest.FixtureRequest) -> Path:
    """Path where failure-only diagnostics for the current test would land.

    The directory is **not** created here. Creation is the artifact emitter's
    responsibility; on a passing test no directory should appear.
    """
    test_id = request.node.name
    return ARTIFACT_ROOT / test_id


@pytest.fixture
def mesh_cache_dir() -> Path:
    """Resolved tests/_mesh_cache/ path. Created on first use by the cache."""
    return MESH_CACHE_ROOT
