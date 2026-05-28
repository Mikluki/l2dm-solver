# Submission 0002 — Algebraic verification of exact solutions

**Status:** done
**Predecessors:** none. Ran independently of 0001.
**Successors:** every subsequent verification-problem submission depended on these conclusions.

## Goal

Independently re-derive each exact solution in `verification.md` Problems 1–5 and confirm it satisfies the stated PDE, BCs, source compatibility, and (where applicable) interface and flux continuity — before any of those solutions got baked into a `Problem` definition. Pure derivation; no `src/` or `tests/` edits.

## Acceptance — outcome

1. **All five problems have a recorded derivation** in `docs/derivations/algebraic-verification.md`. ✓
2. **Every check listed in `verification.md`** (PDE, BCs, source compatibility, interface continuity, flux balance, structural properties like κ₂-independence in Problem 2, rate-degradation in Problem 5) is explicitly addressed. ✓
3. **No edits to `verification.md`.** Discrepancies logged in the derivations doc; resolution by human update to `verification.md` where needed. ✓
4. **All discrepancies resolved.** ✓ — Problems 1, 2, 3, 5 verified cleanly; Problem 4's missing $O((d_{\text{sep}}/2R_{\text{out}}))^2$ outer-boundary residual was added to `verification.md` § Problem 4 (now documents both finite-separation and boundary/image error mechanisms), which 0008's acceptance threshold accommodates.

## What shipped

```
docs/derivations/algebraic-verification.md — one section per Problem 1–5
```

The discipline this submission established and that successor briefs followed: the worker re-derives without editing `verification.md` or `physics.md`; discrepancies land in `docs/derivations/` and route through human resolution, not silent patching.
