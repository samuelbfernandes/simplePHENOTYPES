# DECISIONS.md — simplePHENOTYPES Architectural Decisions
> Fill these in yourself — these are your calls, not Claude's.
> Each decision here unlocks the next phase of work.

---

## DECISION-001: Package Decomposition Strategy

**Question:** Do we split simplePHENOTYPES into separate packages, and if so how?

**Options:**
- A) Two packages: `genoUtils` (format converters) + `simplePHENOTYPES` (simulation)
- B) One package, internal namespace separation only
- C) Three packages: `genoUtils` + `simplePHENOTYPES` + `simBreed` (multi-gen)

**Constraints to consider:**
- Users who only want format conversion shouldn't need to install a phenotype simulator
- CRAN package interdependencies add maintenance burden
- Multi-generation simulation may attract a completely different user base

**Decision:** [ YOUR ANSWER HERE ]

**Rationale:** [ WHY ]

**Date:** [ DATE ]

---

## DECISION-002: isqg Integration Strategy

**Question:** How do we integrate multi-generation simulation given isqg is archived on CRAN?

**Options:**
- A) Depend on archived isqg — simple but fragile; CRAN may reject; no future maintenance
- B) Fork isqg and maintain it ourselves — controllable but adds a whole package to maintain
- C) Absorb isqg's C++ core into our Rust port — we own it; no archived dep; best long-term

**Recommendation:** Option C. isqg's bitwise chromosome representation maps naturally
to Rust's type system. The recombination algorithms (already in C++) port well to Rust.
This eliminates the archived dependency risk and gives you full control.

**Decision:** [ YOUR ANSWER HERE ]

**Rationale:** [ WHY ]

**Date:** [ DATE ]

---

## DECISION-003: Backward Compatibility

**Question:** Does the refactored package maintain the current `create_phenotypes()` API?

**Options:**
- A) Full backward compatibility — existing user scripts keep working; harder to clean up
- B) Breaking v2 with migration guide — clean API; forces users to update scripts
- C) Compatibility shim — new API internally, old API calls map to new ones with warnings

**Constraints:**
- Package is published; breaking changes need a major version bump and CRAN notice
- If the API has bugs, backward compat means keeping the buggy interface

**Decision:** [ YOUR ANSWER HERE ]

**Rationale:** [ WHY ]

**Date:** [ DATE ]

---

## DECISION-004: Multi-Generation Scope

**Question:** Is multi-generation simulation part of simplePHENOTYPES v2 or a separate package?

**Options:**
- A) Part of simplePHENOTYPES — keeps everything together; familiar to current users
- B) Separate `simBreed` package — clear separation; different audience; standalone utility

**Decision:** [ YOUR ANSWER HERE ]

**Rationale:** [ WHY ]

**Date:** [ DATE ]

---

## DECISION-005: Expression-Based Simulation Timeline

**Question:** When do we implement expression-based simulation?

**Options:**
- A) Part of this redesign (v2) — tackle now while architecture is being rethought
- B) Deferred to v3 — finish the marker-based refactor and multi-gen first; expression is different enough to warrant its own design cycle

**Recommendation:** Option B. Expression-based simulation has different enough input
semantics (continuous vs discrete, transcript vs marker) that designing it alongside the
marker-based redesign risks muddying both. Finish one cleanly, then spec the other.

**Decision:** [ YOUR ANSWER HERE ]

**Rationale:** [ WHY ]

**Date:** [ DATE ]

---

## DECISION-006: Utility Package Name

**Question:** What is the standalone format-converter package called?

**Options:** `genoUtils`, `markerTools`, `genomicIO`, `genoConvert`

**Constraints:** Check CRAN for name conflicts before deciding.
```bash
# Check if name is taken:
available::available("genoUtils")
```

**Decision:** [ YOUR ANSWER HERE ]

**Date:** [ DATE ]

---

## Decision Log Summary

| ID | Question | Decision | Date |
|----|----------|----------|------|
| 001 | Package decomposition | | |
| 002 | isqg integration | | |
| 003 | Backward compatibility | | |
| 004 | Multi-gen scope | | |
| 005 | Expression simulation timeline | | |
| 006 | Utility package name | | |
