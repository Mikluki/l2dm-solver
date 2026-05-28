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
- `done` — acceptance passes; brief is historical from here.

Predecessor/successor links are append-only.

## Planner / worker boundary

The planner sets the physics and the verification contract; the worker decides how to code it.

- Planner: physics correctness, acceptance thresholds, structural assertions, cross-doc references, decisions with project-wide consequences (anything the doc-set encodes).
- Worker: file and module layout, API shapes, parameter names, dataclass design, local implementation patterns inside one submission's surface.

If a worker style choice starts looking like a physics or contract decision (e.g. a "convenience" coordinate-based subdomain assignment that `architecture.md` § Key decisions — material interfaces forbids), the worker stops and surfaces.

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

## Post-accept compaction

Briefs accumulate detail while in flight: "Decisions resolved", "Decisions deferred", "Implementation notes (worker)", forced-failure logs, convergence tables. That detail is load-bearing only until status flips to `accepted`. After that it becomes archival noise for every agent that reads the directory.

On status flip to `accepted` (or `done`), compact the brief to roughly:

- Head matter (Status / Predecessors / Successors).
- Goal (one paragraph).
- Acceptance — the criteria as they were, with each line marked ✓ or noted as superseded.
- A short "What shipped" pointer: file paths, key commit hash if non-obvious, the one or two facts a future reader actually needs (e.g. "final L² rate 1.97 on transfinite mesh; convergence table in commit `<sha>`"). Brief-specific anti-creep items (decisions that exist *only* in this brief and aren't pinned in `architecture.md` § Out of scope) fold in here as one-liners so the signal survives compaction.

The "Out of scope" section is **dropped** on accept. Project-wide exclusions are owned by `architecture.md` § Out of scope, which carries the motives and revisit triggers. Per-brief restatement was rotting on Part-1↔Part-2 boundaries and duplicating the canonical list.

Move forced-failure logs and convergence tables into the commit message that flipped the status. Briefs are not changelogs.

Worker decisions that surface a lesson worth keeping (e.g. 0003's "transfinite mesh was needed because L² error sat at the noise floor on the unstructured mesh") go into `docs/open-questions.md` if still open, or into `docs/inspector.md` / `architecture.md` § Key decisions if they generalize. They do not live in an accepted brief.

Compaction is the author's job, applied as part of the same commit that flips status. The brief that compacts itself is doing its job; the brief that grows after acceptance is not.
