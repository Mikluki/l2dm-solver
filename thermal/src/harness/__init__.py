# ABOUTME: Package marker for src.harness. The harness owns the refinement-
# study driver, the L^2 / H^1 error norms (via scikit-fem Functional so the
# norms share quadrature with assembly), and failure-only artifact emission.
# It must not bake in PDE-specific assumptions; Part 2's integral form must
# plug in unchanged.
