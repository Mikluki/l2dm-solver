# Submission conventions

Meta-doc for `docs/submissions/`. The minimum brief structure lives in `CLAUDE.md` § Workflow → Submissions; this file records what's emerged in practice and is not load-bearing enough to push into `CLAUDE.md`.

## Head matter

Every brief opens with:

```
**Status:** <proposed | accepted | in-progress | done>
**Predecessors:** <list of NNNN, or `none`>
**Successors:** <list of NNNN, or `none`>
```

Status lifecycle:
- `proposed` — drafted, awaiting planner sign-off.
- `accepted` — planner approved; worker may implement.
- `in-progress` — implementation started but acceptance held by a discrepancy needing human adjudication. See 0002.
- `done` — acceptance passes; brief is historical from here.

Predecessor/successor links are append-only.

## Planner / worker boundary

The planner sets the physics and the verification contract; the worker decides how to code it.

- Planner: physics correctness, acceptance thresholds, structural assertions, cross-doc references, decisions with project-wide consequences (anything the doc-set encodes).
- Worker: file and module layout, API shapes, parameter names, dataclass design, local implementation patterns inside one submission's surface.

If a worker style choice starts looking like a physics or contract decision (e.g. a "convenience" coordinate-based subdomain assignment that ADR-0003 forbids), the worker stops and surfaces.

## Don't restate the docs

Briefs are pointers, not summaries. If the spec lives in `verification.md`, link to it; do not paraphrase. Restated spec drifts from the source and pollutes the worker's reading context.

A brief earns its lines by saying things *not* already in the docs the worker reads anyway: pre-resolved decisions, acceptance thresholds specific to this submission, forced-failure checks, predicted first-run failure modes, cross-references to algebraic verification.

## Optional sections (use when load-bearing, omit otherwise)

Beyond `CLAUDE.md`'s mandated goal / doc-refs / deliverable / acceptance / out-of-scope:

- **Decisions resolved before implementation** — items pre-decided in planning. The worker does not re-litigate silently. Short, self-contained, rationale > completeness.
- **Decisions left to the worker** — coding choices the planner explicitly does *not* pre-resolve. Listing them prevents the worker from inferring planner intent where none exists.
- **Pre-implementation checkpoint** — "before you start, confirm these" items.
- **Done definition** — what it takes to flip status to `accepted` or `done`, when non-trivial.

Empty or perfunctory sections are noise; omit them.

## Recurring acceptance patterns

Two patterns recur in Part 1 verification briefs. Use these by name rather than inventing new forms.

- **Convergence-rate acceptance.** ≥ 3 mesh sizes; fitted rate ≥ a floor *and* within a window of the theoretical value. The "within a window" half catches threshold-skimming (a rate that scrapes the floor is suspect, not success).
- **Forced-failure (loud-fail) check.** Run once before commit: deliberately break the thing under test (swap κ tags, perturb a Dirichlet value, multiply κ by 2 in the assembler). Confirm the relevant test fails and the failure-artifact bundle is emitted. A clean rerun produces no artifact directory. Catches "test passes for the wrong reason" without leaving the broken code in the tree.

## Naming

`NNNN-slug.md`, four-digit append-only numbering, hyphenated slug describing the unit of work (not the date). This file is `_conventions.md` — the underscore prefix marks it as meta, not a numbered submission.
