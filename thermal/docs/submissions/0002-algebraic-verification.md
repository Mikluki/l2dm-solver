# Submission 0002 — Algebraic verification of exact solutions

**Status:** in-progress (Problem 4 discrepancy pending adjudication; see `docs/derivations/algebraic-verification.md` § Problem 4)
**Predecessors:** none. Runs independently of 0001.
**Successors:** every subsequent verification-problem submission depends on the conclusions being signed off.

## Goal

Independently re-derive each exact solution in `verification.md` Problems 1–5 and confirm it satisfies the stated PDE, BCs, source compatibility, and (where applicable) interface and flux continuity — before any of those solutions get baked into a `Problem` definition.

The work is **pure derivation**. No `src/` or `tests/` edits.

## Relevant core-doc sections

- `verification.md` Problems 1–5 — solutions, BCs, sources, acceptance criteria to be re-derived.
- `physics.md` § Equations → "What Part 1 solves" — the operator under test, $-\nabla\cdot(\kappa\nabla T) = Q$.
- `CLAUDE.md` § Doc-editing rules — `verification.md` is **human-owned**; discrepancies are logged, not patched.

## Decisions resolved before this submission

1. **`verification.md` is human-owned and read-only for this worker.** Discrepancies are logged in `docs/derivations/algebraic-verification.md` and surfaced here, not patched in the source.
2. **Derivations live in `docs/derivations/algebraic-verification.md`.** One file, one section per problem. The brief points to it; the brief does not duplicate the math.
3. **Symbolic computation tools (SymPy, pen and paper) are scratch.** Only the human-readable LaTeX/math in the derivations doc lands in the repo.

## Acceptance

1. **All five problems have a recorded derivation** in `docs/derivations/algebraic-verification.md`. Bare "verified" without showing the math is not acceptance.
2. **Every check listed in `verification.md`** (PDE, BCs, source compatibility, interface continuity, flux balance, structural properties like κ₂-independence in Problem 2 or rate-degradation in Problem 5) is explicitly addressed.
3. **No edits to `verification.md`.** Discrepancies are logged in the derivations doc; the submission stays `in-progress` until a human resolves them.
4. **Status transitions:** if all 5 verify cleanly → `done`. If any discrepancy is logged → stays `in-progress`.

## Verification status

- **Problems 1, 2, 3, 5:** verify cleanly against the doc.
- **Problem 4:** discrepancy in the asymptotic justification. The doc's $O(R_0/d)$ scaling is a real finite-separation/interface mechanism, but the documented bounded-domain superposition also has a common outer-boundary residual controlled by $a/R_{\text{out}}$ with leading symmetric residual $O((a/R_{\text{out}})^2)$. Full derivation and proposed resolution in `docs/derivations/algebraic-verification.md` § Problem 4.

This submission stays at `in-progress` until `verification.md` Problem 4 is clarified to either define the approximation using the common outer boundary so the test isolates $O(R_0/d)$ interface interaction, or state that the measured discrepancy includes both error mechanisms.

## Out of scope

- Any Python code, `Problem` objects, harness code, tests.
- Recommending changes to the verification-problem set itself.
- Editing `verification.md` or `physics.md` (human-owned).
- Harness machinery unit tests (separate submission).
