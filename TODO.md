# Later work

## Closed-form infinitesimal `inCircle` (degree ? 2)

**Status:** deferred. Current mitigation is `Polynomial::trimNearZero` in `findVirtualEvents` and after infinitesimal flip `build_trigger` (see plan `infinitesimal_incircle_trim_150331d4`).

**Why:** With sites of the form `p` or `p + D·?` (one shared separation direction `D`), the Guibas–Stolfi `inCircle` polynomial has algebraic degree ? 2 and the `?³` coefficient is identically 0. Generic polynomial expansion still leaves FP residue on the cubic lead; Eigen then solves a near-singular cubic and invents bad roots.

**Done already**

- [x] `Polynomial::trimNearZero(abs=1e-12, rel=1e-10)`
- [x] Call it from `findVirtualEvents` (not primary `findEvents`)
- [x] Call it after infinitesimal flip trigger build (before monitor / enqueue)
- [x] Monitor dumps for trigger-input site polys + tiny-coeff warnings

**TODO later (only if trim is not enough)**

- [ ] Derive explicit `I(?) = c0 + c1·? + c2·?²` for shared-`D` / `?_i ? {0,1}` sites (differences `d,e,f` relative to query site `p`)
- [ ] Implement a special-case infinitesimal `inCircle` (or gate in `build_trigger` when all linear parts are parallel to one `D`) that never allocates degree 3/4
- [ ] Optionally same idea for `ccw` under shared `D` (usually constant or linear)
- [ ] Re-run monitored flip (e.g. sites 39/337 @ t?15, infinitesimal pass) and confirm degree ? 2 without relying on trim
- [ ] Add a unit test that builds shared-`D` trajectories, compares generic `inCircle`+`trimNearZero` vs closed form coeffs

**Pointers**

- Predicates: `kinDS/KineticDelaunayEventPredicates.hpp` (`inCircle`, `ccw`)
- Trajectories: `KineticDelaunay::buildInfinitesimalSiteTrajectory`
- Trim: `kinDS/Polynomial.hpp` (`trimNearZero`)
- Math sketch: plan `infinitesimal_incircle_trim_150331d4` (`[?³] ? 0` identity)
