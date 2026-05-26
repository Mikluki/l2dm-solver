# CLAUDE.md

## Read these first

The four core docs in `docs/` are the source of truth for this project. Do not work from summaries in this file.

- **`docs/physics.md`** — what the code computes, the equations, the symbol glossary, the regime assumptions.
- **`docs/verification.md`** — the test problems, the harness requirements, the Problem-definition interface.
- **`docs/architecture.md`** — stack, module layout, key abstractions, decisions, out-of-scope.
- **`docs/open-questions.md`** — live questions with resolution criteria.

Decisions live in `docs/decisions/` as numbered ADRs. Submissions live in `docs/submissions/`.

---

## Project disposition

This project is dual-natured. The harness — the test infrastructure built around scikit-fem — is the deliverable, and **it must be rigorous**. The solver code that the harness exercises is a means to that end.

The asymmetry: the harness gets tight tests and careful interfaces; the solver gets the smallest reasonable implementation that passes the harness.

When in doubt, ask: does this make the harness more trustworthy, or does it make the solver more clever? Prefer the first; resist the second.

---

## Coding standards

### General

- Prioritize simple, clean, maintainable solutions over clever ones.
- Make the smallest reasonable changes to achieve the desired outcome.
- Never make unrelated code changes; note tech debt inline as a comment instead.
- Preserve existing comments unless you can prove they are actively false.
- Avoid temporal context in comments; make them evergreen.
- Don't take shortcuts that create future burden. Note tech debt inline rather than silently working around it.
- Too many features or too much info is useless if it's noise. Be respectful of attention span.
- Use consistent explicit argument names across functions, especially when they share the same physical meaning (`mesh_size`, `kappa`, `R_inner`, etc.).

### Python

- Formatter: `black`, 88 char line length.
- Function args: explicit names, `Path` objects only (no raw strings for paths), dataclass when config exceeds 3 params.
- Use `logging`, not `print()`. Tests use `caplog` if they need to assert on log content. Pytest handles log capture; no per-script file handlers needed.
- Every code file starts with a module-level header comment block:
  ```python
  # ABOUTME: <one-paragraph description of file purpose, I/O, side effects>
  ```
- Caching is mandatory for heavy computation. Mesh generation is the dominant cost; cache `.msh` files by SHA-256 of geometry parameters per ADR-0007.

### Section separators

In files with multiple logical sections, use these headers:

```python
# ============================================================================
# CONFIG
# ============================================================================

# ============================================================================
# LOGGING
# ============================================================================

# ============================================================================
# FUNCTIONS
# ============================================================================

# ============================================================================
# CLASSES
# ============================================================================
```

Use only the sections that apply. Small files with one or two sections don't need all five.

### Units

- SI throughout (see `physics.md` symbol glossary).
- The 2D-vs-3D conductivity confusion is the project's most likely unit footgun. $\sigma_{2d}$ has units W/K; $\kappa$ has units W/(m·K). Name variables accordingly (`sigma_2d` vs `kappa`), never reuse a name across the two.

---

## Environment

- Python via `uv`, project-local virtual environment at `./.venv/`.
- All dependencies pinned in `pyproject.toml` (ADR-0001).
- Do not add runtime dependencies without an ADR.

---

## Workflow

### Submissions

- Each unit of work is a submission brief in `docs/submissions/NNNN-slug.md`.
- A submission brief states: goal (one sentence), relevant core-doc sections, deliverable (files created/modified, artifact produced), acceptance (concrete check), out-of-scope-for-this-submission.
- A submission is done when its acceptance check passes and the relevant core docs are updated if the submission changed anything they specify.

### Test-first cadence

The harness is the deliverable. The cadence for any submission that touches a verification problem:

1. Write the `Problem` definition (geometry, $\kappa$, $Q$, BCs, exact solution).
2. Write the pytest test wiring the problem through the harness.
3. Run the test and watch it fail in a **specific, predicted way** (e.g., "fails because the solver module doesn't exist yet").
4. Implement the missing piece.
5. Watch the test pass for the right reasons (check the convergence rate, not just pass/fail).

A test that passes on first run is suspicious — either it was already passing for unrelated reasons or it's not actually checking what you think.

### When a decision needs making mid-task

Stop. Surface the decision. Do not silently choose.

- If the decision is small and reversible, propose two options in the submission's response and let the human pick.
- If the decision is significant, draft an ADR in `docs/decisions/` (proposed status) and surface it.
- Never amend `architecture.md` to record a decision without a corresponding ADR.

### When scope creep tempts

Each "out of scope" item in `architecture.md` has a stated motive. If the current task seems to require violating one, stop and surface — do not extend scope silently. The motive usually points at why extending is worse than it looks.

---

## Documentation discipline

Never create bloat files. Every markdown file is a maintenance burden.

**Before creating any .md file, ask:**

1. Does this content already exist elsewhere?
2. Will someone actually read this, or is it info noise?
3. Should this live in code comments instead?
4. Could this be a 5-line section in an existing file instead of a new file?

**Allowed structure:**

- `README.md` at root — overview + quick start only.
- `CLAUDE.md` at root — this file.
- `docs/physics.md`, `docs/verification.md`, `docs/architecture.md`, `docs/open-questions.md` — the four core docs. Fixed set.
- `docs/decisions/NNNN-slug.md` — ADRs, one per decision.
- `docs/submissions/NNNN-slug.md` — submission briefs.

**Do not create:**

- New `.md` files at root level. Anything project-wide goes in the four core docs or an ADR.
- Per-script documentation. Use ABOUTME comments instead.
- Findings or results docs as ongoing artifacts. If a one-off investigation needs writing up, name it `FINDINGS_{timestamp}.md` and treat it as disposable.
- Narrative walkthroughs. Same.

**Single-task constraint:** if a submission produces multiple new `.md` files at root, the submission was too broad or the files should be consolidated.

### Doc-editing rules

- **Agent may edit:** `docs/open-questions.md` (append questions, mark resolved), `docs/decisions/` (propose new ADRs as `proposed` status), `docs/submissions/` (create new briefs, update status of in-flight ones).
- **Agent must propose, human approves:** `docs/architecture.md` (any structural change), ADR transitions from `proposed` to `accepted`.
- **Human-owned, agent must not silently edit:** `docs/physics.md`, `docs/verification.md`. If the agent believes either is wrong or incomplete, raise it in `docs/open-questions.md` and stop.

---

## Out-of-scope reflex

The "Out of scope" section in `architecture.md` (no parallelism, no AMR, no time dependence, no GPU, no alternative element orders mid-Part-1, no integral form yet) is binding. Each item has a stated motive and a "when to revisit" clause. If a task seems to require violating one, the response is to stop and surface, not to extend.

Additional behavioral rules not in `architecture.md`:

- Do not add configuration knobs that no test exercises.
- Do not refactor working code in unrelated submissions. Note the concern in `open-questions.md` instead.
