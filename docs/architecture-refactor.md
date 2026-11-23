# kinDS architecture refactor plan

Status: planned, not started (after the oak 10 boundary-mesh fixes are committed). Written 2026-11-23 against the state after the boundary-mesh fixes on branch
`324-bifurcation-for-kinetic-voronoi-meshing`.

## Goals

- Modular functionality: independent **workloads** (segment meshlets, shell mesh, debug SVG, ...) that
  `KineticDelaunay` schedules without knowing what they do.
- Each workload declares which event types it needs; events nobody needs are not computed.
- Adding a workload takes at most one registration line.
- Smaller, single-purpose classes and files; diagnostics separated from production code paths.
- Subfolders for sources and headers.

## Current state (measured)

- One `EventCallback*` per event manager and one `CallbackManager*` (lifecycle hooks: init, finalize, retriangulated,
  graph cut, finished-branch retriangulation, failsafe). `SegmentBuilder` is that manager and owns everything:
  volumetric segment meshlets, bark meshlets, the boundary mesh (`boundary_mesh_`, ~60 uses across
  `SegmentBuilder`, flip and radius callbacks) and the debug SVGs. Its six event callbacks are `friend` classes.
- Every event type is always computed. Crossing and subdivision events only matter for segment meshlets, yet the
  crossing data (117 uses in `KineticDelaunay.cpp`) is always maintained.
- Debug SVG export is scattered over the callbacks and ~27 places in `KineticDelaunay.cpp`.
- Shared between segment and boundary mesh: mesh placement (`getPointInMeshSpace`: hinge, side planes, split blend)
  and component boundary polygons (`updateBoundary`).
- 83 files flat in `kinDS/kinDS/`; `SegmentBuilder.cpp` 16k lines, `KineticDelaunay.cpp` 10.6k lines.
- By method name, ~1.7k lines of `KineticDelaunay.cpp` and ~950 of `SegmentBuilder.cpp` are validation, checks,
  debug, logging, SVG or dumps (plus inline diagnostics).
- `KineticDelaunay.cpp` topics by method name: runtime branches ~1.7k lines, components/splits ~1.5k,
  sections/separation ~1.05k, crossings ~300.
- `TreeMesher::Settings` has 58 fields; the CLI has 37 flags; `main.cpp` is 1.8k lines with 19 helpers, many of
  them developer-only exports.
- 354 `static_cast<size_t>(-1)` sentinels; every id is a plain `size_t`.
- EvoEngine includes 12 kinDS headers, some internal (`SegmentBuilder.hpp`, `KineticDelaunayHelpers.hpp`,
  `MeshingBuffer.hpp`).

## Target architecture

### Events

- **Core events** (always run; they change triangulation, inside flags or components): flip, radius, section,
  separation, failsafe.
- **Optional events** (computed, and their data maintained, only if a registered workload needs them): crossing,
  subdivision.

### Workloads

- `Workload` interface with optional typed handlers, e.g. `onFlip(const FlipEvent&, Phase)`, plus optional lifecycle
  hooks (init, finalize, retriangulated, graph cut, branch finished, failsafe). Typed handlers replace
  `beforeEvent(Event&)` + `dynamic_cast` (12 casts today) and the `friend` callback classes.
- Dispatch order: registration order for "before", reverse for "after", so a debug workload registered first sees
  the state before and after all others.
- Event needs are deduced at compile time: a template adapter detects the implemented handlers with C++20
  `requires` and builds the event mask. Explicit `static constexpr` override for special cases.
- Registration: `kinetic_delaunay.addWorkload<ShellMeshWorkload>(services, settings)` at the composition root
  (`TreeMesher`), driven by settings. No static self-registration: kinDS is a static library and the linker may drop
  unreferenced registrars. Optional name -> factory table for selecting workloads from the CLI / EvoEngine (one
  table line per workload).
- Threading contract: handlers may be called concurrently for different runtime branches.

### Initial workloads

- `SegmentMeshWorkload`: what remains of `SegmentBuilder` (segment and bark meshlets, strips, crossings).
  Needs all events. Further split into `StripTracker`, `BarkMeshlets` (today's "intersection meshes"),
  `ClosingCaps` (~1k lines), `PlaceholderResolver` (flexible placeholders, ~570 lines), `MeshletStore` (storage,
  UVs, export, collapse post-processing) and a small crossing / Voronoi-vertex buffer cache.
- `ShellMeshWorkload` (today's `BoundaryMesh`): the shell mesh plus the boundary parts of the flip, radius and section handlers, the cap
  pass and `advanceBoundaryMesh`. Needs flip, radius, section.
- `DebugSvgWorkload`: before/after SVGs (`SegmentBuilderVisualDebug`). Queries segment-workload state only if that
  workload is registered. Internal error dumps of `KineticDelaunay` become diagnostics (see below).
- Diagnostics as observer workloads: `ValidationWorkload` (`--check-intersections`, sites in hull, meshlet
  validation), error-dump writer, optionally `DebugExportWorkload` (profile planes, raw meshlets, scope exports).

### Split `KineticDelaunay`

Thin coordinator owning:

- kinetic triangulation: graph, inside flags, event managers, predicates, scheduler;
- `ComponentTracker`: components, outline polygons, splits/merges, graph cuts, finished-branch retriangulation;
- `RuntimeBranchManager`: `BranchScheduler`, pending splits, parallel scopes;
- `SeparationManager`: infinitesimal separation, virtual offsets;
- `CrossingTracker`: crossing data and crossing events (optional, enabled on demand).

### Shared services

- `StrandPositionModel`: one answer to "where is strand s at t". It must keep three coordinate spaces explicitly
  distinct (separate functions / types, never implicit conversions):
  1. **profile / kinetic coordinates** used for the triangulation and event predicates (reference-frame policy:
     hold, blend, hinge unfold);
  2. **shifted coordinates** for infinitesimal time (profile coordinates plus the virtual separation offset);
  3. **world / mesh coordinates** used for mesh placement (hinge, side planes, split blend).
  It unifies
  today's `StrandTree` frame transforms, `KineticDelaunay` `getPointAt` variants / post-split frame transitions
  (hold, blend) and `SegmentBuilder` placement (hinge, side planes, split blend). Recent hinge and frame bugs all came
  from the interaction of these three.
- Component outlines (from `ComponentTracker`).
- Workloads receive services plus a read-only `KineticDelaunay` view through one context object (no friend access).

### Boundary front conventions

`BoundaryFront` helper owning the implicit rules that caused all boundary-mesh bugs fixed on 2026-11-23: fronts live
on the outside-facing half-edge, left corner = origin, corners swap when a triangle joins vs leaves. Operations:
front for a face edge (joining / leaving), split a front at an apex, merge two fronts with a corner bridge. Radius
cases 0-3 and both flip paths use it instead of re-deriving twins and corners.

### Folder layout

Include paths stay under `kinDS/...`.

```
kinDS/core/         Polynomial, Trajectory, EventTime/KineticAlgorithm, schedulers, Statistics, Logger, SlotTable
kinDS/geometry/     HalfEdgeDelaunayGraph, Delaunator2D, Polygon2D, PlaneProjector, SimilarityBlend, GeometryUtils
kinDS/kinetic/      KineticDelaunay (+ collaborators), events, predicates, StrandTree, Workload interface/dispatcher
kinDS/services/     StrandPositionModel, component outlines
kinDS/mesh/         VoronoiMesh, ObjExporter, MeshIntersection, MeshCGAL, Validator, MeshingBuffer
kinDS/workloads/    segment_mesh/, boundary_mesh/, debug_svg/, validation/
kinDS/diagnostics/  dumps, consistency checks, formatting helpers
kinDS/pipeline/     TreeMesher
include/kinDS/      public API (TreeMesher, workload registration, settings, mesh types, exporter)
tools/              developer exports, regression harness
```

### Further cleanups

- **Typed ids** (`StrandId`, `HalfEdgeId`, `FaceId`, `VoronoiVertexId`, `ComponentId`, `InputBranchId`,
  `RuntimeBranchId`) with `kInvalid` / `std::optional` instead of `size_t(-1)`. Introduce gradually, starting with
  the workload interfaces.
- **Naming** (decided): "boundary" currently means bark meshlets, the `BoundaryMesh` shell and component outlines.
  Rename to `BarkMeshlets`, `ShellMesh`, `ComponentOutline`; `SegmentBuilder` -> `SegmentMeshWorkload`
  (`BoundaryMeshWorkload` accordingly becomes `ShellMeshWorkload`).
- **Settings per workload**; debug options in one `DebugSettings`; `TreeMesher::Settings` keeps only the workload
  list and kinetic settings. Table-driven CLI parsing (flag, help, target field).
- **Developer tools out of the library CLI**: profile-plane / raw-meshlet / scope exports into `tools/` or a
  `DebugExportWorkload`.
- **Regression harness as a CTest target** (`tools/regress`): 9 oaks + oak 10 alpha sweep, byte-identical comparison
  against a baseline, mesh checks (watertight, manifold, orientation, long edges, zero-area faces). Replaces ad-hoc
  scripts.
- **Public API boundary**: EvoEngine includes only `include/kinDS/`.
- **No hidden global state**: replace `thread_local` flags (e.g. `in_parallel_section_finish`,
  `emit_collapsed_triangles`) and static globals with explicit context parameters.
- **Legacy and failsafe paths** (decided): keep `--legacy-hinge` (`legacy_mirrored_hinge` in `PlaneProjector`,
  `StrandTree`), the failsafe paths and the "shouldn't occur but does" branches for now. They may be moved behind the
  diagnostics layer but not removed. Commented-out code can go.
- **Metadata**: structured fields serialized only when metadata is enabled, instead of JSON strings built in hot
  paths.
- **File size guideline**: no source file above ~2k lines.

## Migration steps

Every step keeps outputs byte-identical to the baseline (combined and boundary meshes), except step 8, which may only
change performance.

0. Baseline: record meshes and timings for the 9 oaks and the oak 10 alpha sweep.
1. Regression harness as CTest target (`tools/regress`).
2. Extract diagnostics (validation, dumps, debug formatting) into `kinDS/diagnostics/`; production paths keep only
   hooks.
3. Folder restructure: `git mv`, includes, CMake globs, EvoEngine includes. No code changes.
4. `Workload` interface and dispatcher; wrap today's `SegmentBuilder` + six callbacks as one workload; fold
   `CallbackManager` into workload hooks.
5. Extract shared services (`StrandPositionModel`, component outlines); split `KineticDelaunay` collaborators.
   Largest step; several commits.
6. `ShellMeshWorkload` with the `BoundaryFront` helper; `boundary_mesh.obj` byte-identical.
7. `DebugSvgWorkload` and `ValidationWorkload`.
8. Event demand: skip crossing / subdivision events and crossing data when no workload needs them; verify a
   boundary-only run gives the same `boundary_mesh.obj`, faster.
9. Ergonomics and docs: compile-time event-mask deduction, `addWorkload<T>`, optional name table, per-workload
   settings, table-driven CLI, public `include/kinDS/`, README example for a new workload, EvoEngine integration.

Typed ids and renames are applied file by file as they are touched.

## Risks

- Core events may implicitly depend on crossing data (graph cuts, radius boundary transitions); step 8 must verify
  it can be disabled without changing core results.
- Boundary and segment code currently run interleaved inside one callback; splitting may change order within an
  event. Byte-identical checks guard every step.
- `SegmentBuilder`'s size makes step 5 the biggest; split it into several commits.

## Decisions (2026-11-23)

- All suggestions accepted, including the `StrandPositionModel` (with the three coordinate spaces kept distinct) and
  the renames.
- Legacy (`--legacy-hinge`) and failsafe code stay; they may be hidden behind the diagnostics layer.

## Open decisions

1. Bark meshlets belong to `SegmentMeshWorkload`; `ShellMeshWorkload` is only the shell mesh (`boundary_mesh.obj`).
   Confirm.
2. One-line `addWorkload<T>()` vs automatic static registration (linker caveat).
3. May `TreeMesher`'s settings/API change for EvoEngine (workload list instead of individual flags)?
4. Debug SVG: keep drawing segment-workload state, or only the kinetic structure?
5. Should the exported file name change with the rename (`boundary_mesh.obj` -> `shell_mesh.obj`)?
