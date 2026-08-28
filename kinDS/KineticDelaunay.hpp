#pragma once
#include "HalfEdgeDelaunayGraph.hpp"
#include "HalfEdgeDelaunayGraphToSVG.hpp"
#include "KineticAlgorithm.hpp"
#include "Logger.hpp"
#include "Polynomial.hpp"
#include "ProgressBar.hpp"
#include "StrandTree.hpp"
#include <filesystem>
#include <format>
#include <functional>
#include <glm/gtx/exterior_product.hpp>
#include <glm/gtx/string_cast.hpp>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace kinDS
{
struct BoundaryPoint
{
  size_t vertex_id;
  size_t he_id;
  glm::dvec2 p;
};

/// Data recorded when a component split awaits graph cut / retriangulation.
struct PendingBranchSplit
{
  size_t parent_component_id = static_cast<size_t>(-1);
  /// Runtime branch that is splitting (the cut / finalize key). One runtime branch may span multiple kinetic components.
  size_t parent_runtime_branch = static_cast<size_t>(-1);
  /// Reference strand from the parent component used for the frozen frame.
  size_t reference_vertex = static_cast<size_t>(-1);
  /// Full pre-split parent-component strand list.
  std::vector<size_t> frozen_parent_strands;
  /// Component ids from the split: retained parent first, then split-off components.
  std::vector<size_t> split_component_ids;
  glm::dvec2 old_branch_centroid { 0.0, 0.0 };
  glm::dvec2 new_branch_centroid { 0.0, 0.0 };
  glm::dvec2 separation_direction { 0.0, 0.0 };
  /// Kinetic time at which sites are frozen for virtual / infinitesimal separation.
  double infinitesimal_t_event = 0.0;
  /// True while cross-piece events are driven by virtual time along @ref separation_direction.
  bool infinitesimal_active = false;
  /// Bumped on activate-replace and finalize so stale infinitesimal events skip.
  uint64_t infinitesimal_epoch = 0;
  /// When true, finalize checks (convex seams → graph cut, section
  /// boundary force-cut) are skipped: the pending pieces are no longer topologically distinct (e.g. after a
  /// radius outside→inside component consolidate). Reactivation clears this flag when the same input-branch
  /// split is noted again.
  bool on_hiatus = false;
  /// True after a same-t @ref SeparationEvent with @c after_radius_dispatch_order has been enqueued to apply the
  /// graph cut once radius (and other lower-order) handlers at this time finish. Dedups re-activate calls.
  bool cut_event_queued = false;
  /// Sorted unique input-branch ids spanned by the pending pieces (identity for hiatus reactivation).
  std::vector<size_t> input_branch_ids;
};

/// Registry of pending branch splits and per-strand lookup into the frozen parent frame.
struct PendingBranchSplitState
{
  static constexpr size_t no_parent_component = static_cast<size_t>(-1);

  void clear();
  void resetStrandLookup(size_t strand_count);
  PendingBranchSplit& getOrCreate(size_t parent_component_id);
  const PendingBranchSplit* findByParent(size_t parent_component_id) const;
  void registerStrandsForSplit(
    size_t parent_component_id, const std::vector<std::vector<size_t>>& new_components);
  const std::vector<size_t>* frozenStrandsForStrand(size_t strand_id) const;

  std::unordered_map<size_t, PendingBranchSplit> by_parent_;
  std::vector<size_t> strand_parent_component_;
};

/// Survives @ref clearPendingBranchSplits: hold shared frame through @c hold_end_height, then blend on @c blend_section.
struct PostSplitFrameTransition
{
  size_t parent_component_id = static_cast<size_t>(-1);
  /// Shared input-branch frame held until @ref hold_end_height.
  size_t common_reference_branch = 0;
  /// Last height still fully expressed in the common frame (S+1 when the split is in section S).
  size_t hold_end_height = 0;
  /// Section whose piece lerps common(@c hold_end_height) → native(@c hold_end_height + 1).
  size_t blend_section = 0;
  std::vector<size_t> strand_ids;
  bool rotation_warning_emitted = false;
};

inline constexpr double kPostSplitFrameBlendRotationWarnDegrees = 20.0;

struct PostSplitFrameTransitionState
{
  static constexpr size_t no_transition = static_cast<size_t>(-1);

  void clear();
  void resetStrandLookup(size_t strand_count);
  PostSplitFrameTransition& add(PostSplitFrameTransition transition);
  const PostSplitFrameTransition* findForStrand(size_t strand_id) const;
  void expireBeforeHeight(size_t height);

  std::vector<PostSplitFrameTransition> transitions_;
  /// Parallel to strand id → index into @ref transitions_, or @ref no_transition.
  std::vector<size_t> strand_transition_index_;
};

static glm::dvec2 polygonCentroid(const std::vector<BoundaryPoint>& polygon)
{
  double A = 0.0;
  glm::dvec2 C { 0.0, 0.0 };

  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& p = polygon[i].p;
    const glm::dvec2& q = polygon[(i + 1) % n].p;

    double cross = glm::cross(p, q);
    A += cross;
    C += (p + q) * cross;
  }

  A *= 0.5;

  if (std::abs(A) < 1e-12)
    return C; // degenerate polygon

  return C / (6.0 * A);
}

static std::vector<size_t> buildComponentMap(const std::vector<std::vector<size_t>>& components, size_t vertex_count)
{
  std::vector<size_t> component_map(vertex_count);

  for (size_t i = 0; i < components.size(); i++)
  {
    for (const auto v : components[i])
    {
      component_map[v] = i;
    }
  }

  return component_map;
}

/**
 * \brief Class for computing the Delaunay triangulation of a set of cubic Hermite splines.
 *
 * This follows Guibas, L.J., Mitchell, J.S.B., Roos, T. (1992).
 * Voronoi diagrams of moving points in the plane. In:
 * Schmidt, G., Berghammer, R. (eds) Graph-Theoretic Concepts in Computer Science. WG 1991.
 * Lecture Notes in Computer Science, vol 570. Springer, Berlin, Heidelberg.
 * https://doi.org/10.1007/3-540-55121-2_11
 */
class KineticDelaunay
{
 public:
  using Event = KineticAlgorithm::Event;
  using EventCallback = KineticAlgorithm::EventCallback;
  using EventManager = KineticAlgorithm::EventManager;

  enum class ComponentSplitPolicy
  {
    Retriangulate,
    InPlaceCut,
  };

  class FlipEvent;
  class RadiusEvent;
  class CrossingEvent;

  class CallbackManager
  {
   public:
    virtual ~CallbackManager() = default;

    virtual void init() { }
    virtual void finalize(double t) { }
    virtual void onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots) { (void)t; (void)prev_face_slots; (void)prev_he_slots; }
    virtual void onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots) { (void)t; (void)prev_face_slots; (void)prev_he_slots; }
    /// Invoked immediately before a targeted pending runtime-branch graph cut in @ref applyPendingRuntimeBranchSplit.
    virtual void onBeforeComponentGraphSplit(double t) { (void)t; }
    /// Invoked when frozen-site infinitesimal separation becomes active (before virtual event recompute).
    virtual void onInfinitesimalSeparationActivated(size_t parent_component_id, double t)
    {
      (void)parent_component_id;
      (void)t;
    }
  };

  class FlipEvent;
  class RadiusEvent;
  class CrossingEvent;
  class SectionEvent;

  class FlipEventManager;
  class RadiusEventManager;
  class CrossingEventManager;
  class SectionEventManager;
  class SubdivisionEvent;
  class SubdivisionEventManager;
  class SeparationEvent;
  class SeparationEventManager;

  /// Forward declaration; full definition in `KineticDelaunayCrossingEvent.hpp`.
  struct CrossingData;

  struct ComponentData
  {
    std::vector<std::vector<size_t>> components;
    std::vector<size_t> component_map;
    // [component_index][boundary_no][point_no] - the first boundary is the outer one, any additional ones are holes in
    // the polygon
    std::vector<std::vector<std::vector<BoundaryPoint>>> component_boundaries;
    std::vector<glm::dvec2> component_centroids;
    std::vector<double> component_last_updated;
  };

  ComponentData component_data;

  /// Runtime-branch bookkeeping. Mirrors @ref ComponentData (a per-strand id lookup plus a per-branch strand list) but
  /// is keyed by runtime branch ids, which are allocated monotonically and never reused. Also tracks split lineage:
  /// the parent each branch split off from, and, for a branch whose split is scheduled but not yet executed, the child
  /// branch ids that split will produce.
  struct RuntimeBranchData
  {
    static constexpr size_t no_branch = static_cast<size_t>(-1);

    /// Per-strand runtime branch id (index = strand id); @ref no_branch when unassigned/retired.
    std::vector<size_t> branch_map;
    /// Strand ids per runtime branch id (index = branch id); inverse of @ref branch_map over live strands. Rebuilt via
    /// @ref rebuildBranchStrandLists.
    std::vector<std::vector<size_t>> branches;
    /// Parent runtime branch each branch split from (index = branch id); @ref no_branch for original branches.
    std::vector<size_t> parent;
    /// Liveness per runtime branch id (index = branch id).
    std::vector<bool> alive;
    /// For a parent branch mid-split, the child branch ids whose split has not been executed yet. Populated when the
    /// pending split is first noted and cleared when the split is executed.
    std::unordered_map<size_t, std::vector<size_t>> pending_splits;
    /// Monotonic runtime-branch id allocator; ids are never reused.
    size_t next_branch_id = 0;

    /// Reserve per-branch bookkeeping so @p branch_id is addressable.
    void ensureBranchCapacity(size_t branch_id)
    {
      if (branch_id + 1 > branches.size())
      {
        branches.resize(branch_id + 1);
      }
      if (branch_id + 1 > parent.size())
      {
        parent.resize(branch_id + 1, no_branch);
      }
      if (branch_id + 1 > alive.size())
      {
        alive.resize(branch_id + 1, false);
      }
    }

    /// Allocate a fresh, live runtime branch id, optionally recording its @p parent_branch.
    size_t allocate(size_t parent_branch = no_branch)
    {
      const size_t branch_id = next_branch_id++;
      ensureBranchCapacity(branch_id);
      alive[branch_id] = true;
      parent[branch_id] = parent_branch;
      return branch_id;
    }

    /// Runtime branch id for @p strand_id, or @ref no_branch when out of range.
    size_t branchForStrand(size_t strand_id) const
    {
      return strand_id < branch_map.size() ? branch_map[strand_id] : no_branch;
    }

    /// Rebuild @ref branches from the current @ref branch_map (skipping @ref no_branch entries).
    void rebuildBranchStrandLists()
    {
      for (std::vector<size_t>& strands : branches)
      {
        strands.clear();
      }
      for (size_t strand_id = 0; strand_id < branch_map.size(); ++strand_id)
      {
        const size_t branch_id = branch_map[strand_id];
        if (branch_id == no_branch)
        {
          continue;
        }
        ensureBranchCapacity(branch_id);
        branches[branch_id].push_back(strand_id);
      }
    }

    /// Record that @p parent_branch will split into @p child_branches (not yet executed).
    void noteSplit(size_t parent_branch, std::vector<size_t> child_branches)
    {
      pending_splits[parent_branch] = std::move(child_branches);
    }

    /// Clear the pending-split record for @p parent_branch once its split has been executed.
    void completeSplit(size_t parent_branch) { pending_splits.erase(parent_branch); }

    /// Child branch ids of @p parent_branch whose split has not been executed, or nullptr if none pending.
    const std::vector<size_t>* pendingChildBranches(size_t parent_branch) const
    {
      const auto it = pending_splits.find(parent_branch);
      return it == pending_splits.end() ? nullptr : &it->second;
    }

    /// True iff @p branch_id is a child produced by a split that has not been executed yet (so it is still part of a
    /// larger, not-yet-separated "unsplit" branch).
    bool isPendingSplitChild(size_t branch_id) const
    {
      if (branch_id >= parent.size() || parent[branch_id] == no_branch)
      {
        return false;
      }
      const std::vector<size_t>* children = pendingChildBranches(parent[branch_id]);
      if (children == nullptr)
      {
        return false;
      }
      for (size_t child : *children)
      {
        if (child == branch_id)
        {
          return true;
        }
      }
      return false;
    }

    /// The enclosing unsplit branch id for @p branch_id: walks up pending-split parents until reaching a branch that is
    /// not itself a pending-split child. Returns @p branch_id when it is already unsplit.
    size_t unsplitBranchId(size_t branch_id) const
    {
      size_t current = branch_id;
      while (isPendingSplitChild(current))
      {
        current = parent[current];
      }
      return current;
    }

    /// Strand ids of the unsplit branch rooted at @p branch_id: its own strands plus those of all (recursively) pending
    /// split-off children.
    std::vector<size_t> unsplitBranchStrands(size_t branch_id) const
    {
      std::vector<size_t> strands;
      std::vector<size_t> pending = { branch_id };
      while (!pending.empty())
      {
        const size_t current = pending.back();
        pending.pop_back();
        if (current < branches.size())
        {
          strands.insert(strands.end(), branches[current].begin(), branches[current].end());
        }
        if (const std::vector<size_t>* children = pendingChildBranches(current))
        {
          pending.insert(pending.end(), children->begin(), children->end());
        }
      }
      return strands;
    }
  };

  RuntimeBranchData runtime_branch_data_;

  std::pair<glm::dvec2, glm::dvec2> computeAngularBisector(size_t he_id, double t) const;

  std::pair<double, double> delaunayVoronoiEdgeIntersection(
    size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const;

  /**
   * \brief Compute the half-edges crossed by the Voronoi edge between the given start point and destination, starting
   * from the given face.
   *
   * Currently assumes that the start face is finite and that the destination lies outside of only one edge of the start
   * triangle.
   */
  std::pair<std::vector<size_t>, std::vector<double>> computeCrossedHalfEdges(
    size_t start_face_id, const glm::dvec2& destination, const glm::dvec2& start_point, double t) const;

 private:
  friend class SectionEventManager;
  friend class SubdivisionEvent;
  friend class SeparationEvent;

  StrandTree branch_trajs;
  HalfEdgeDelaunayGraph graph;
  std::unique_ptr<KineticAlgorithm> kinetic_algorithm_;
  Statistics statistics_;
  bool collect_statistics_ = false;
  // Reused managers: one per event type.
  // Kept as pointers to avoid forcing complete manager types in this header.
  std::unique_ptr<FlipEventManager> flip_event_manager_;
  std::unique_ptr<RadiusEventManager> radius_event_manager_;
  std::unique_ptr<CrossingEventManager> crossing_event_manager_;
  std::unique_ptr<SectionEventManager> section_event_manager_;
  std::unique_ptr<SubdivisionEventManager> subdivision_event_manager_;
  std::unique_ptr<SeparationEventManager> separation_event_manager_;
  /// Strand mesh subdivision schedule; consumed on the first @ref enqueueScheduledSubdivisionEvents in @ref compute.
  std::vector<std::pair<size_t, double>> subdivision_schedule_;
  CallbackManager* callback_manager_ = nullptr;
  size_t sections_advanced = 0; // Counter for the number of sections advanced; starts at @ref start_section_
  /// First section index to initialize and process (inclusive). Default 0.
  size_t start_section_ = 0;
  /// Last section index to process (inclusive). Empty means @ref getSectionCount() - 1.
  std::optional<size_t> end_section_;
  double cutoff; // Cutoff radius for boundary events
  std::vector<bool> face_inside; // Tracks whether faces are inside or outside the boundary

  std::vector<std::vector<size_t>> branches; // track which vertices/splines belong to which branch
  std::vector<glm::dvec2> dummy_boundary;
  bool add_dummy_boundary;
  size_t prev_component_count = 1;
  PendingBranchSplitState pending_branch_splits_;
  /// Post-cut frame hold/blend; survives @ref clearPendingBranchSplits.
  PostSplitFrameTransitionState post_split_frame_transitions_;
  ComponentSplitPolicy component_split_policy_ = ComponentSplitPolicy::InPlaceCut;
  double separation_offset_scale_ = 500.0;
  /// Virtual time of the infinitesimal event currently being handled (0 outside that path).
  mutable double current_infinitesimal_t_ = 0.0;
  /// Set while @ref EventManager::computeEvents runs with an @ref InfinitesimalComputeContext
  /// (via @ref ScopedInfinitesimalEventCompute). Prefer the computeEvents parameter at call sites.
  mutable bool computing_infinitesimal_events_ = false;
  mutable double infinitesimal_recompute_min_x_ = 0.0;
  mutable double infinitesimal_schedule_t_ = 0.0;
  mutable size_t infinitesimal_schedule_parent_ = static_cast<size_t>(-1);
  mutable uint64_t infinitesimal_schedule_epoch_ = 0;
  std::vector<EventTime> quadrilateral_last_updated;
  std::vector<EventTime> face_last_updated;
  bool on_the_fly_boundary = true;
  std::optional<std::filesystem::path> visual_debug_output_root_;
  bool visual_debug_enabled_ = false;
  bool error_files_enabled_ = false;
  /// When true, visual-debug SVGs use pending split-off child runtime branch folders (and own strand sets) as soon as
  /// a radius event notes the pending split — not only after the graph cut.
  bool visual_debug_separate_pending_splits_ = false;
  std::optional<double> flip_polynomial_dump_target_time_;
  std::optional<size_t> flip_polynomial_dump_target_half_edge_;
  bool diagnostics_enabled_ = false;
  /// When true, after each kinetic event verify every live site lies in its component's graph convex hull.
  /// Default off; CLI @c --check-sites-in-hull.
  bool sites_inside_convex_hull_check_enabled_ = false;

  // crossing-related data (see public `CrossingData` forward declaration above).

  // Owned by `CrossingEventManager`. This is only an alias reference so the existing
  // code can keep using the name `crossing_data`.
  CrossingData& crossing_data;

  glm::dvec3 computeVoronoiVertexHomogenous(size_t voronoi_vertex_id, double t,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  void reassignVoronoiVerticesOnBoundary(size_t he_id, double t,
    std::optional<InfinitesimalComputeContext> infinitesimal = std::nullopt);

  void reassignVoronoiVerticesInQuadrilateral(size_t quad_index, double t,
    const std::map<size_t, size_t>& pre_flip_quad_faces,
    std::optional<InfinitesimalComputeContext> infinitesimal = std::nullopt);

  void precomputeStep(double t);

  /**
   * Piecewise-linear site motion on [@p section, @p section + 1] using @ref StrandTree::getPiecePolynomial with the
   * same per-strand reference branch as @ref getPointAt for times in that interval.
   */
  Trajectory<2> getSitePiecePolynomial(size_t strand_id, size_t section, double schedule_time) const;

  void growGraphSlotArrays();
  size_t findContainingTriForVoronoiVertex(size_t voronoi_vertex_id, double t) const;
  void initializeFaceState(size_t face_index, double t);
  void initializeNewFacesAfterGraphUpdate(double t, size_t first_new_face_slot);
  void clearPendingBranchSplits();
  /// Drop pending kinetic entries for @p parent_runtime_branch_id and mark its runtime children as final.
  void completePendingRuntimeBranchSplit(size_t parent_runtime_branch_id, double t);
  /// Build a runtime-branch cut map that severs only @p target_parent_runtime_branch from its pending children.
  /// Other pending-split children are collapsed to their unsplit parent so they are not cut yet.
  std::vector<size_t> buildRuntimeBranchCutMapForParent(size_t target_parent_runtime_branch) const;
  /// Set @ref prev_component_count from live components minus still-uncut pending split extras.
  void syncPrevComponentCountWithPendingSplits();
  /// Clear strand list and boundary/centroid/last_updated for an emptied (merged-away) component slot.
  void clearComponentSupportData(size_t component_id);
  /// Section index used to resolve input-branch ids when noting / reactivating a pending split at @p split_time.
  size_t pendingSplitBranchSection(double split_time) const;
  /// Sorted unique input-branch ids among non-dummy strands in @p strand_groups at @p branch_section.
  std::vector<size_t> collectInputBranchIdsForStrandGroups(
    const std::vector<std::vector<size_t>>& strand_groups, size_t branch_section) const;
  /// On-hiatus pending split with the same @c input_branch_ids set, if any.
  PendingBranchSplit* findHiatusPendingSplitWithInputBranches(const std::vector<size_t>& input_branch_ids);
  /// Assign split-off strands to pending child runtime branches (allocate only when none exist yet).
  void assignPendingSplitChildRuntimeBranches(PendingBranchSplit& split, size_t parent_component_id,
    const std::vector<size_t>& split_component_ids, size_t branch_section, double split_time, bool allow_allocate);
  /// At section @p section_index, look ahead to height @p section_index + 1 and register hold/blend when a
  /// live component's strands already occupy multiple input branches there (upcoming split).
  void registerUpcomingPostSplitFrameTransitions(size_t section_index);
  void registerPostSplitFrameTransition(size_t parent_component_id, double split_time,
    const std::vector<size_t>& affected_strands);
  void maybeWarnPostSplitFrameBlendRotation(PostSplitFrameTransition& transition) const;
  const PostSplitFrameTransition* postSplitFrameTransitionForStrand(size_t strand_id) const;
  bool pendingSplitSeamsAreConvex(size_t parent_component_id, double t) const;
  void handleSeparationEventAtTime(size_t parent_component_id, double t);
  /// Apply the in-place graph cut (or retriangulation) for one pending parent *runtime branch* only.
  void applyPendingRuntimeBranchSplit(double t, size_t parent_runtime_branch_id);
  /// Activate frozen-site virtual separation, or handle convex seams: enqueue a same-t @ref SeparationEvent when
  /// @p apply_cut_now is false, otherwise apply the graph cut immediately.
  void activateInfinitesimalSeparationOrApplyCut(size_t parent_component_id, double t, bool apply_cut_now = false);
  /// Recompute flip/radius/crossing for all cross-piece simplices using virtual site trajectories.
  /// Used to seed the infinitesimal event queue on activate; post-event refresh uses the same local
  /// neighbor paradigm as regular events (under @ref ScopedInfinitesimalEventCompute).
  /// Roots are scheduled at kinetic @p t with @c infinitesimal_t in (@p min_virtual_x, +inf).
  void recomputeEventsAfterInfinitesimalSeparation(size_t parent_component_id, double t, double min_virtual_x);
  /// If seams are convex and not on hiatus: bump epoch, clear active, apply graph cut.
  bool maybeFinalizeInfinitesimalSeparation(size_t parent_component_id, double t);
  void collectSeparationRecomputeTargets(size_t parent_component_id, std::unordered_set<size_t>& affected_quads,
    std::unordered_set<size_t>& affected_faces) const;
  const PendingBranchSplit* activeSeparationForStrand(size_t strand_id) const;
  const PendingBranchSplit* infinitesimalSplitForStrand(size_t strand_id) const;
  Trajectory<2> buildInfinitesimalSiteTrajectory(size_t strand_id, double t_event,
    const std::vector<size_t>& event_strand_ids) const;
  Trajectory<2> addSeparationOffsetToPiecePolynomial(
    const Trajectory<2>& base, size_t strand_id, size_t section, double schedule_time) const;
  void debugSeparationTrackedFlipProbe(
    size_t parent_component_id, double t, const char* phase, size_t even_half_edge_id = 22) const;
  /// Diagnostic: logs per-input-branch profile-plane transform column metrics (lengths, in-plane orthogonality,
  /// handedness) for the branches involved in a pending split, to check whether each transform is decomposable into
  /// scalar * orthogonal (a prerequisite for skew-free frame switching via @ref PlaneProjector).
  void logSplitTransformOrthonormalityDiagnostic(size_t parent_component_id, double t) const;
  /// Strand ids whose min input branch defines the motion frame for @p strand_id.
  void collectReferenceBranchStrandPool(size_t strand_id, std::vector<size_t>& pool) const;
  void updateRuntimeBranchMapFromInputBranches(double t);
  void updateRuntimeBranchMapFromLiveGraph(double t, const char* trigger);
  size_t allocateRuntimeBranch();
  void markRuntimeBranchDead(size_t runtime_branch_id, double t, const char* reason,
    std::optional<size_t> input_branch_id = std::nullopt);
  void logRuntimeBranchEvent(double t, const std::string& event_line);
  bool runtime_branch_log_header_written_ = false;
  void retireFinishedInputBranches(double t);
  void validateFinishedInputBranchMatchesRuntime(size_t section, size_t input_branch_id) const;
  std::vector<std::vector<size_t>> extractGraphConnectedComponents() const;
  std::vector<size_t> extractGraphConnectedComponent(size_t u, std::vector<bool>& visited) const;
  size_t getRuntimeBranchIdForFace(size_t face_id) const;
  void onGraphRetriangulated(double t, size_t prev_face_slots, size_t prev_he_slots);
  /// Post-cut bookkeeping. Does not clear pending splits or remap runtime branches — the caller finalizes the
  /// targeted pending parent via @ref completePendingRuntimeBranchSplit. Optional live-graph remap is off by default.
  void onGraphCutApplied(double t, size_t prev_face_slots, size_t prev_he_slots, bool update_runtime_branch_map = false,
    const HalfEdgeDelaunayGraph::RuntimeBranchSplitResult* split_result = nullptr);
  /// After a graph cut, recompute flip events for (1) infinite edges that bordered live+tombstoned faces and
  /// (2) finite capped outer edges. Stamps @ref quadrilateral_last_updated so stale schedules are ignored.
  void refreshEventsAfterGraphCut(double t, const HalfEdgeDelaunayGraph::RuntimeBranchSplitResult& split_result);

  void handleEvents();

  /// Enqueues one @ref SubdivisionEvent per entry in @ref subdivision_schedule_ (called once from @ref compute).
  void enqueueScheduledSubdivisionEvents();

  size_t getBranchIndex(size_t strand_id, size_t t) const;

  const std::vector<std::vector<size_t>>& getBranches(size_t t) const;

  const std::vector<size_t>& getBranchStrands(size_t t, size_t branch_id);

  std::vector<double> findEvents(
    Polynomial& event_trigger, double min_fraction, bool only_positive_to_negative = false);
  /// Like @ref findEvents but accepts any finite root strictly greater than @p min_x (no section upper bound).
  std::vector<double> findVirtualEvents(
    Polynomial& event_trigger, double min_x, bool only_positive_to_negative = false);

 public:
  KineticDelaunay(const StrandTree& branch_trajs, double cutoff, bool add_dummy_splines);
  ~KineticDelaunay();

  bool isDummyBoundary(size_t v) const;

  /// True when @p strand_id has at least one incident half-edge in the live Delaunay graph.
  bool isStrandLiveInGraph(size_t strand_id) const;

  bool computeBoundaryOnTheFly() const;

  /// Site position with independent frame / separation controls.
  /// @param apply_reference_transform Map into the component reference-branch frame (@ref StrandTree::evaluateTransformed).
  /// @param include_virtual_offset Add @ref separationOffsetAt.
  glm::dvec2 getPointAt(size_t v, double t, bool apply_reference_transform = true,
    bool include_virtual_offset = true) const;

  glm::dvec2 getPointAt(double t, size_t v, bool apply_reference_transform = true,
    bool include_virtual_offset = true) const;

  /// Kinetic Delaunay / SVG coordinates: reference-branch frame with virtual separation (never object-space).
  glm::dvec2 getPointInDelaunaySpace(size_t v, double t) const;
  glm::dvec2 getPointInDelaunaySpace(size_t v, double t, size_t reference_branch) const;

  /** Branch frame used by getPointAt for sites in the same component as @p strand_id. */
  size_t getReferenceBranch(size_t strand_id, double t) const;

  /// Input-branch section index at the upper bound of an event computation interval.
  size_t inputBranchSectionIndexAtIntervalUpperBound(double event_interval_upper_bound) const;

  /// Smallest live strand id sharing @p strand_id's runtime branch; used for transform frame lookup.
  size_t representativeStrandIdForRuntimeBranch(size_t strand_id) const;

  /// One shared reference branch for all @p strand_ids (e.g. flip quadrilateral vertices).
  size_t getSharedReferenceBranchForStrands(const std::vector<size_t>& strand_ids, double branch_lookup_time) const;

  /// Distinct input-branch ids among @p event_strand_ids at @p event_interval_upper_bound.
  std::vector<size_t> collectDistinctInputBranchesForEventTrigger(
    const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const;

  /// True when @p event_strand_ids for one event trigger involve more than one input branch at @p event_interval_upper_bound.
  bool eventTriggerUsesSharedTransformedFrame(
    const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const;

  /// Shared reference input branch for a multi-branch event trigger at @p event_interval_upper_bound.
  size_t sharedReferenceBranchForEventTrigger(
    const std::vector<size_t>& event_strand_ids, double event_interval_upper_bound) const;

  /// Piece polynomial for one strand in an event trigger (frame choice uses eventIntervalUpperBound(@p schedule_time)).
  Trajectory<2> getSitePiecePolynomialForEventStrands(size_t strand_id, size_t section, double schedule_time,
    const std::vector<size_t>& event_strand_ids) const;

  /// Parameter of the Voronoi/Delaunay crossing along the Delaunay edge, computed in Delaunay profile space.
  double delaunayVoronoiEdgeIntersectionParameter(size_t delaunay_edge_id, size_t voronoi_edge_id, double t) const;

  /// Recompute stored @c delaunay_edge_param for all intersections on one Delaunay edge at @p t (kinetic space).
  void refreshDelaunayEdgeIntersectionParams(size_t delaunay_edge_id, double t);
  /// Like @ref refreshDelaunayEdgeIntersectionParams, then sort the CrossingData list by recomputed param.
  /// Safe after flips when all mesh-pair edge links on the edge are unset (@c -1).
  void refreshAndSortDelaunayEdgeIntersectionParams(size_t delaunay_edge_id, double t);

  /// Recompute @c delaunay_edge_param for intersections on all three edges of Delaunay face @p face_id.
  void refreshTriangleDelaunayEdgeIntersectionParams(size_t face_id, double t);

  glm::dvec2 getPointAtWithReferenceBranch(size_t v, double t, size_t reference_branch,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  Trajectory<2> getSitePiecePolynomialWithReferenceBranch(
    size_t strand_id, size_t section, size_t reference_branch, double schedule_time) const;

  std::vector<glm::dvec2> getPointsAt(double t, bool apply_reference_transform = true,
    bool include_virtual_offset = true) const;

  glm::dvec3 getPointInObjectSpace(size_t v, double t) const;

  const StrandTree& getStrandTree() const;

  /// Input branch ids whose last real-strand section is @p t (same criterion as @ref retireFinishedInputBranches).
  std::vector<size_t> inputBranchesFinishingAtSection(double t) const;

  void setVisualDebugOutputRoot(const std::filesystem::path& root);
  const std::optional<std::filesystem::path>& getVisualDebugOutputRoot() const;
  /// When true, emit full visual-debug artifacts (segmentbuilder-driven SVGs, branch-split dumps, …).
  void setVisualDebugEnabled(bool enabled);
  bool isVisualDebugEnabled() const;
  /// When true, emit failure SVG/TXT dumps (@c --error-files). Also implied by @ref isVisualDebugEnabled.
  void setErrorFilesEnabled(bool enabled);
  bool isErrorFilesEnabled() const;
  bool shouldDumpErrorFiles() const;
  /// When true, write visual-debug SVGs into pending child branch folders (with that branch's strands only) while a
  /// split is still pending. Default off; CLI @c --svg-separate-pending-splits.
  void setVisualDebugSeparatePendingSplits(bool enabled);
  bool visualDebugSeparatePendingSplits() const;
  /// Path to @c branches.txt under @ref getVisualDebugOutputRoot when visual debug is enabled.
  /// Truncated whenever @ref setVisualDebugOutputRoot is called so each run starts with a fresh log.
  std::optional<std::filesystem::path> getRuntimeBranchLogPath() const;

  void setFlipPolynomialDumpTargetTime(std::optional<double> target_time);
  const std::optional<double>& getFlipPolynomialDumpTargetTime() const;

  void setFlipPolynomialDumpTargetHalfEdge(std::optional<size_t> half_edge_id);
  const std::optional<size_t>& getFlipPolynomialDumpTargetHalfEdge() const;

  void computeComponentData(double t);

  /// On a radius outside→inside transition, if the triangle's vertices belong to different kinetic components,
  /// fold every higher id into the smallest id: copy strands (updating @c component_map while iterating), clear
  /// emptied slots, then recompute the survivor's boundary / centroid / last_updated. Slots are never compacted.
  /// Pending splits that touch any absorbed / survivor id are put @c on_hiatus (finalize suppressed).
  void consolidateComponentsAtTriangle(const std::array<int, 3>& triangle_vertices, double t);

  /// Put every pending split whose @c split_component_ids intersects @p component_ids onto hiatus.
  void putPendingSplitsOnHiatusTouchingComponents(const std::vector<size_t>& component_ids, double t);

  /// True if any kinetic pending split for @p parent_runtime_branch_id is currently on hiatus.
  bool isPendingRuntimeBranchOnHiatus(size_t parent_runtime_branch_id) const;

  /// True once @ref component_data reflects the current @ref HalfEdgeDelaunayGraph after the latest section
  /// retriangulation. While false, connected components may already be split in component data but the mesh
  /// is still the pre-split triangulation (see @ref prev_component_count).
  bool isGraphRetriangulatedForComponents() const;

  void setComponentSplitPolicy(ComponentSplitPolicy policy) { component_split_policy_ = policy; }
  ComponentSplitPolicy getComponentSplitPolicy() const { return component_split_policy_; }
  void setSeparationOffsetScale(double scale) { separation_offset_scale_ = scale; }
  double getSeparationOffsetScale() const { return separation_offset_scale_; }
  /// Section fraction in [0,1] where the active separation ramp ends, if it ends inside the section.
  std::optional<double> separationRampEndSectionFraction(size_t strand_id, size_t section) const;

  /// Virtual separation offset for @p strand_id at @p t (added when @c include_virtual_offset is true).
  /// Uses @ref current_infinitesimal_t_ (see @ref ScopedCurrentInfinitesimalTime).
  glm::dvec2 separationOffsetAt(size_t strand_id, double t) const;

  /// Kinetic @p real_time paired with the infinitesimal coordinate currently in effect for handle/export.
  EventTime eventTimeAt(double real_time) const { return EventTime(real_time, current_infinitesimal_t_); }

  /// Infinitesimal coordinate bound during event handling / SVG export (0 outside those paths).
  double currentInfinitesimalTime() const { return current_infinitesimal_t_; }

  /// Temporarily sets @ref current_infinitesimal_t_ so @ref getPointAt / SVG exports apply the matching virtual shift.
  class ScopedCurrentInfinitesimalTime
  {
   public:
    ScopedCurrentInfinitesimalTime(KineticDelaunay& kd, double infinitesimal_t)
      : kd_(kd)
      , previous_(kd.current_infinitesimal_t_)
    {
      kd_.current_infinitesimal_t_ = infinitesimal_t;
    }
    ScopedCurrentInfinitesimalTime(KineticDelaunay& kd, EventTime t)
      : ScopedCurrentInfinitesimalTime(kd, t.infinitesimal_time)
    {
    }
    ~ScopedCurrentInfinitesimalTime() { kd_.current_infinitesimal_t_ = previous_; }

    ScopedCurrentInfinitesimalTime(const ScopedCurrentInfinitesimalTime&) = delete;
    ScopedCurrentInfinitesimalTime& operator=(const ScopedCurrentInfinitesimalTime&) = delete;

   private:
    KineticDelaunay& kd_;
    double previous_;
  };

  /// Enables virtual / infinitesimal polynomial scheduling for @ref findVirtualEvents and site piece polys.
  /// Restores prior flags on destruction. No-op (inactive) if the pending split is missing or not active.
  class ScopedInfinitesimalEventCompute
  {
   public:
    ScopedInfinitesimalEventCompute(KineticDelaunay& kd, size_t parent_component_id, double t, double min_virtual_x);
    ~ScopedInfinitesimalEventCompute();

    ScopedInfinitesimalEventCompute(const ScopedInfinitesimalEventCompute&) = delete;
    ScopedInfinitesimalEventCompute& operator=(const ScopedInfinitesimalEventCompute&) = delete;

    bool active() const { return active_; }

   private:
    KineticDelaunay& kd_;
    bool active_ = false;
    bool previous_computing_ = false;
    double previous_min_x_ = 0.0;
    double previous_schedule_t_ = 0.0;
    size_t previous_parent_ = static_cast<size_t>(-1);
    uint64_t previous_epoch_ = 0;
  };

  /// Record a pending branch split until the next graph cut/retriangulation.
  /// If a matching on-hiatus pending split already exists for the same input-branch set, reactivates it
  /// (clears @c on_hiatus, refreshes piece ids / runtime child assignment) without re-scheduling separation.
  void notePendingBranchSplit(size_t parent_component_id, double split_time,
    const std::vector<size_t>& pre_split_parent_strands, const std::vector<std::vector<size_t>>& new_components,
    const std::vector<size_t>& split_component_ids);
  /// Recompute @ref PendingBranchSplit separation centroids from current @ref component_data centroids.
  void refreshPendingSplitSeparationCentroids(size_t parent_component_id, double t);
  /// Centroid-delta direction for @p split at @p t in the same point frame used by event polys / offsets.
  /// When @p shared_reference_branch is set, all sites are evaluated in that shared transformed frame;
  /// otherwise each site uses @p apply_reference_transform with its native reference branch (or local support).
  glm::dvec2 computeSeparationDirection(const PendingBranchSplit& split, double t, bool apply_reference_transform,
    std::optional<size_t> shared_reference_branch = std::nullopt) const;
  /// After a pending split is noted: if seams are already convex, enqueue a same-t @ref SeparationEvent
  /// (@ref SeparationEvent::after_radius_dispatch_order) so the graph cut runs after radius meshing; otherwise
  /// start infinitesimal virtual separation immediately (graph stays connected).
  void maybeScheduleSeparationOrApplyPendingSplit(size_t parent_component_id, double split_time);

  /// Pending split data keyed by the parent component id, if recorded.
  std::optional<PendingBranchSplit> getPendingBranchSplit(size_t parent_component_id) const;

  /// Calls @p visitor for each pending branch split (parent component id, split data).
  void visitPendingBranchSplits(
    const std::function<void(size_t parent_component_id, const PendingBranchSplit& split)>& visitor) const;

  /**
   * Best (largest-area) seam outline loop walked from the precomputed @p start_edges (see
   * @ref collectFutureBranchSeamStartEdges). The walk follows runtime branch ids alone; it performs no seam/partner
   * lookups of its own.
   */
  std::vector<BoundaryPoint> collectFutureRuntimeBranchOutline(
    const std::vector<size_t>& start_edges, double t) const;

  /// One outer outline per future runtime branch listed in the pending split metadata.
  std::vector<std::vector<BoundaryPoint>> collectPendingSplitBranchOutlines(
    size_t parent_component_id, double t) const;

  /// Delaunay simplices whose flip/radius/crossing events are recomputed when another separation
  /// iteration is scheduled. Only mixed-shift (retained vs separated) targets.
  VisualDebugHighlight buildSeparationRecomputeHighlight(size_t parent_component_id) const;

  /// Base-to-offset segments for separated strands in an active pending split.
  std::vector<HalfEdgeDelaunayGraphToSVG::SeparationOffsetSegment> collectSeparationOffsetSegments(
    size_t parent_component_id, double t) const;

  /// Smallest @ref StrandTree::getBranchIndex among live strands in @p component_id at @p branch_lookup_height.
  size_t minInputBranchForComponent(size_t component_id, size_t branch_lookup_height) const;

  /// Smallest @ref StrandTree::getBranchIndex among @p strand_ids at @p branch_lookup_height.
  size_t minInputBranchForStrands(const std::vector<size_t>& strand_ids, size_t branch_lookup_height) const;

  CrossingData& getCrossingDataMutable();
  const CrossingData& getCrossingData() const;

  /** Per-edge intersection list consistency; intended after @ref EventCallback::afterEvent (e.g. debug SVG export). */
  void validateCrossingIntersectionInvariants(const char* context, double t) const;

  /** Cached Voronoi-vertex list iterators must reference their own id in the containing triangle list. */
  void validateVoronoiVertexIteratorInvariants(const char* context, double t) const;

  std::vector<HalfEdgeDelaunayGraphToSVG::IntersectionDebugInfo> getCrossingIntersectionDebugData() const;

  const HalfEdgeDelaunayGraph& init(CallbackManager* callback_manager = nullptr);
  void registerSectionEventCallback(EventCallback* callback);
  void registerFlipEventCallback(EventCallback* callback);
  void registerRadiusEventCallback(EventCallback* callback);
  void registerCrossingEventCallback(EventCallback* callback);
  void registerSubdivisionEventCallback(EventCallback* callback);
  void registerSeparationEventCallback(EventCallback* callback);
  /// Pairs `(strand_id, t)` sorted by non-decreasing `t`; enqueued once when @ref compute builds the initial queue.
  void setSubdivisionSchedule(std::vector<std::pair<size_t, double>> schedule);
  void registerEventCallbacks(EventCallback* section_callback, EventCallback* flip_callback,
    EventCallback* radius_callback, EventCallback* crossing_callback);

  const HalfEdgeDelaunayGraph& getGraph() const;

  size_t getSectionCount() const;

  /// Kinetic section window for @ref init / @ref compute. Clamped to start in `[0, getSectionCount())` and
  /// end in `[start, getSectionCount()]`. @p end_section is the exclusive stop / finalize time: section events
  /// run for `[start, end)`, and @ref KineticAlgorithm::processEvents discards anything with
  /// `occurrence_time >= end`. When @p start_section > 0, bootstrap uses that height and input branches there
  /// become runtime branches directly.
  void setSectionRange(size_t start_section, std::optional<size_t> end_section = std::nullopt);
  size_t getStartSection() const { return start_section_; }
  /// Exclusive kinetic stop / finalize time (defaults to @c getSectionCount(), i.e. tree height).
  size_t getEndSection() const;

  Statistics& statistics() { return statistics_; }
  const Statistics& statistics() const { return statistics_; }
  void setCollectStatistics(bool enabled) { collect_statistics_ = enabled; }
  bool collectStatistics() const { return collect_statistics_; }
  /// Live non-dummy strands and alive runtime branches (excludes retired / phased-out).
  std::pair<size_t, size_t> countLiveStrandsAndBranches() const;

  // Computes the Delaunay triangulation of the given splines
  void compute();

  std::vector<size_t> extractConnectedComponent(size_t u, std::vector<bool>& visited) const;

  const std::vector<glm::dvec2>& getDummyBoundary() const;

  std::vector<std::vector<size_t>> checkForSplit(const std::array<int, 3>& tri_vertices, double t) const;
  /// Same split check using an explicit face-inside state, allowing callbacks to predict a radius event's target state.
  std::vector<std::vector<size_t>> checkForSplit(
    const std::array<int, 3>& tri_vertices, const std::vector<bool>& inside_state, double t) const;

  std::vector<std::vector<size_t>> extractConnectedComponents() const;

  std::vector<BoundaryPoint> traverseBoundary(size_t start_he_id, double t,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  std::vector<std::vector<BoundaryPoint>> extractComponentBoundaries(const std::vector<size_t>& component, double t,
    std::vector<bool>& he_visited, bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  std::vector<BoundaryPoint> extractComponentBoundary(const std::vector<size_t>& component, double t,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  const std::vector<bool>& getFacesInside() const;

  bool getFaceInside(size_t face_index) const;

  /** Alpha / radius cutoff used for inside-outside classification (see radius events). */
  double getCutoff() const { return cutoff; }

  void setFaceInside(size_t face_index, bool value, double t);

  /** All three vertices share one input branch with exactly three strands at @p t. */
  bool isMinimalInputBranchTriangle(const std::array<int, 3>& vertices, double t) const;

  /**
   * Runtime branch of @p strand_id (see @ref RuntimeBranchData::branch_map).
   * Not the inside-face component id and not the input @ref StrandTree branch index.
   */
  size_t getRuntimeBranchIdForStrand(size_t strand_id) const;

  /** Inside-face connected component id for @p strand_id. */
  size_t getInsideFaceComponentId(size_t strand_id) const;

  /** Lowest strand id in inside-face component @p component_id (BFS seed order). */
  size_t getComponentLowestStrandId(size_t component_id) const;

  /**
   * Runtime branch of a live half-edge, resolved from its incident face(s).
   * Returns 0 when the half-edge is dead or the branch cannot be resolved.
   */
  size_t getRuntimeBranchIdForHalfEdge(size_t half_edge_id) const;

  /** True iff the live Delaunay graph contains exactly one finite triangle in @p runtime_branch_id. */
  bool runtimeBranchHasSingleFiniteTriangle(size_t runtime_branch_id) const;

  const std::vector<size_t>& getRuntimeBranchMap() const { return runtime_branch_data_.branch_map; }
  const RuntimeBranchData& getRuntimeBranchData() const { return runtime_branch_data_; }

  /// Strand ids of the unsplit runtime branch @p branch_id: the branch's own strands plus those of any pending (not yet
  /// executed) split-off children. Use this to treat a mid-split branch as a single, not-yet-separated branch.
  std::vector<size_t> collectUnsplitRuntimeBranchStrands(size_t branch_id) const
  {
    return runtime_branch_data_.unsplitBranchStrands(branch_id);
  }
  /// Maps a (possibly pending split-off child) runtime branch id to its enclosing unsplit branch id.
  size_t unsplitRuntimeBranchId(size_t branch_id) const { return runtime_branch_data_.unsplitBranchId(branch_id); }
  /// True iff @p branch_id is an unsplit branch (not a pending split-off child).
  bool isUnsplitRuntimeBranch(size_t branch_id) const { return !runtime_branch_data_.isPendingSplitChild(branch_id); }

  bool mustRemainInside(size_t face_index, double t) const;

  bool isOnComponentBoundary(size_t he_id) const;

  bool isOnComponentBoundaryOutside(size_t he_id) const;

  bool isOnFutureBranchSeamForComponent(
    size_t he_id, size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const;
  /// Cross-component (tombstoning-criterion) seam half-edges pointing inward into @p component_id (destination in the
  /// component, origin in a partner). Each is an entry seed for the outline walk and is not itself part of the outline.
  std::vector<size_t> collectFutureBranchSeamStartEdges(
    size_t component_id, const std::unordered_set<size_t>& partner_component_ids) const;
  /// Next seam half-edge after @p he_id: the first edge, rotating around the pivot vertex from @p he_id's `next`, whose
  /// origin and destination lie in the same runtime branch. Uses runtime branch ids only.
  size_t nextOnFutureBranchSeamId(size_t he_id) const;
  /// Walk one seam outline loop. @p start_he_id is a start seed (see @ref collectFutureBranchSeamStartEdges); the walk
  /// advances via @ref nextOnFutureBranchSeamId and closes when a half-edge's destination returns to the seed's
  /// destination.
  std::vector<BoundaryPoint> traverseFutureBranchSeam(size_t start_he_id, double t) const;

  size_t nextOnComponentBoundaryId(size_t he_id) const;

  // CrossingData accessors. Parameter `voronoi_vertex_id` is a Delaunay face index (dual circumcenter),
  // NOT a Voronoi cell / site id (Delaunay vertex index).
  bool isCrossingDataVoronoiVertexRegistered(size_t voronoi_vertex_id) const;
  void requireLiveRegisteredVoronoiVertex(size_t voronoi_vertex_id, const char* context) const;
  size_t getCrossingDataContainingTriId(size_t voronoi_vertex_id) const;
  std::vector<size_t> getCrossingDataVoronoiVerticesInTri(size_t tri_id) const;
  glm::dvec3 getVoronoiVertexHomogeneous(size_t voronoi_vertex_id, double t,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  /**
   * \brief Compute the (possibly clamped) Voronoi vertex position for the Delaunay edge represented by half_edge_id.
   *
   * For finite triangles this is the circumcenter; for infinite / hull cases, this returns a finite point obtained
   * by moving a neighboring circumcenter along a perpendicular direction so it can be used for meshing and
   * intersection computations.
   */
  glm::dvec3 computeVoronoiVertexClampedInfinity(size_t half_edge_id, double t,
    bool apply_reference_transform = true, bool include_virtual_offset = true) const;
  glm::dvec3 computeVoronoiVertexClampedInfinityWithReferenceBranch(size_t half_edge_id, double t,
    size_t reference_branch, bool apply_reference_transform = true, bool include_virtual_offset = true) const;

  /**
   * Debug sanity checks: compare @ref face_inside against circumradius and @ref cutoff (and @ref mustRemainInside).
   * Intended to be called only when @ref SegmentBuilder::diagnostics is enabled.
   */
  /// Sentinel for disabled diagnostic targets. Never matches unset / invalid / infinite entity ids.
  static constexpr size_t kDiagnosticsMonitorDisabledId = static_cast<size_t>(-1);
  static constexpr bool isDiagnosticsMonitorIdEnabled(size_t monitor_id)
  {
    return monitor_id != kDiagnosticsMonitorDisabledId;
  }
  /// True only when both ids are enabled (not @ref kDiagnosticsMonitorDisabledId) and equal.
  /// Disabled / unset / invalid sentinels (-1) never match.
  static constexpr bool matchesDiagnosticsMonitorId(size_t candidate_id, size_t monitor_id)
  {
    return isDiagnosticsMonitorIdEnabled(monitor_id) && isDiagnosticsMonitorIdEnabled(candidate_id)
      && candidate_id == monitor_id;
  }

  /// Set to @ref kDiagnosticsMonitorDisabledId to disable face monitoring.
  static constexpr size_t kDiagnosticsMonitoredFaceId = kDiagnosticsMonitorDisabledId;
  /// Debug target: Voronoi vertex whose crossing-event trigger roots are traced.
  /// Set to @ref kDiagnosticsMonitorDisabledId to disable.
  static constexpr size_t kDiagnosticsMonitoredCrossingVoronoiVertexId = 1734;
  /// Debug target: undirected Delaunay edge id highlighted in crossing trigger logs (optional; not a filter when disabled).
  /// Set to @ref kDiagnosticsMonitorDisabledId to disable.
  static constexpr size_t kDiagnosticsMonitoredCrossingDelaunayEdgeId = kDiagnosticsMonitorDisabledId;
  /// Suspected missed crossing time; crossing diagnostics use inclusive window [@ref kDiagnosticsMonitoredCrossingTime,
  /// @ref kDiagnosticsMonitoredCrossingTime + 1] so section-event recomputes at the lower bound are included.
  static constexpr double kDiagnosticsMonitoredCrossingTime = 10.0;
  static constexpr double kDiagnosticsMonitoredCrossingTimeEpsilon = 0.05;
  /// Debug target: undirected Delaunay edge id for flip-event trigger / handle diagnostics.
  /// Directed half-edge 1158 ⇒ undirected edge 579 (also matches twin 1159).
  static constexpr size_t kDiagnosticsMonitoredFlipDelaunayEdgeId = 3644 / 2;
  /// Flip create/discard / trigger-root logging is constrained to [floor(t), floor(t)+1).
  static constexpr double kDiagnosticsMonitoredFlipTime = 30.0;
  void setDiagnosticsEnabled(bool enabled);
  bool diagnosticsEnabled() const;
  /// Optional per-event sanity check: all live sites lie inside the graph convex hull (same topology as SVG).
  /// Off by default; enable via CLI @c --check-sites-in-hull. Failures log @c KINDS_WARNING (possible bad flip).
  void setSitesInsideConvexHullCheckEnabled(bool enabled);
  bool sitesInsideConvexHullCheckEnabled() const;
  void validateSitesInsideConvexHull(const char* context, EventTime t) const;
  /// Bounds-checked diagnostic id queries; invalid ids are ignored by monitor logging.
  bool isDiagnosticsStrandIdValid(size_t strand_id) const;
  bool isDiagnosticsFaceIdValid(size_t face_id) const;
  bool isDiagnosticsHalfEdgeIdValid(size_t half_edge_id) const;
  bool isDiagnosticsMonitoredFaceValid() const;
  bool isDiagnosticsMonitoredCrossingValid() const;
  /// When diagnostics + monitored VV are enabled, log stored containing triangle vs
  /// @ref findContainingTriForVoronoiVertex (same as initialization) at @p t.
  void logDiagnosticsMonitoredCrossingContainingTriangle(double t, const char* context) const;
  void validateFlipAdjacentFaceInsideConsistency(size_t half_edge_id, double t) const;
  void validateAllFaceInsideStatesAtTime(double t, const char* context) const;
  void logFaceInsideStateAtTime(size_t face_id, double t, const char* context) const;
  void logRadiusEventTriggerRoots(size_t face_id, size_t he_id, double t, double min_fraction,
    Polynomial event_trigger, const std::array<size_t, 3>& strand_ids,
    const std::array<Trajectory<2>, 3>& trajectories) const;
  /// Log all real roots, sign changes, and findEvents filter/enqueue decisions for one crossing trigger.
  void logCrossingEventTriggerRoots(size_t voronoi_vertex_id, size_t he_id, size_t edge_index, double t,
    double min_fraction, const Polynomial& event_trigger, bool only_positive_to_negative) const;
  /// Log all real roots, sign changes, and findEvents filter/enqueue decisions for one flip trigger.
  /// @p trigger_predicate is @c "ccw" (convex-boundary) or @c "inCircle" (interior).
  void logFlipEventTriggerRoots(size_t he_id, double t, double min_fraction, const Polynomial& event_trigger,
    const char* trigger_pass, const char* trigger_predicate) const;
};
} // namespace kinDS