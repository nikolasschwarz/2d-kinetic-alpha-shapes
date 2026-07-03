#pragma once
#include "Delaunator2D.hpp"
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <glm/glm.hpp>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

namespace kinDS
{
// A graph that can represent the delaunay triangulation such that edges are explicitly stored and can be flipped in its
// quadrilateral.
class HalfEdgeDelaunayGraph
{
 public:
  struct HalfEdge
  {
    int origin = -1; // index into vertices
    int next = -1; // index into half_edges, always represents the next half-edge in the face
    int face = -1; // index into faces
    // twin = index ^ 1
  };

  struct Triangle
  {
    std::array<size_t, 3> half_edges;
  };

  struct TriangleKeyHash
  {
    size_t operator()(const std::array<size_t, 3>& k) const noexcept
    {
      size_t h = 1469598103934665603ull; // FNV-1a offset basis
      for (size_t v : k)
      {
        h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      }
      return h;
    }
  };

  static bool isTombstone(const HalfEdge& he) { return he.origin == -1 && he.next == -1 && he.face == -1; }

 private:
  size_t vertex_count = 0; // Number of vertices in the triangulation
  std::vector<Triangle> triangles;
  std::vector<HalfEdge> half_edges; // List of half-edges in the triangulation
  std::vector<size_t> vertex_to_half_edge; // Maps each vertex to one of its outgoing half-edges for easy access
  std::vector<size_t> live_even_half_edges_;
  std::vector<size_t> live_faces_;
  std::vector<uint8_t> half_edge_live_;

  /**
   * @brief Build the data structure from an index buffer of triangles.
   *
   * Constructs the half-edge data structure from a list of triangle indices. Construction is done in O(n) time.
   * Twin edges are stored implicitly by storing them next to each other. The twin index can be computed as `index ^ 1`.
   *
   * @param index_buffer index buffer of triangles, size must be multiple of 3 to be valid and each index must be in
   * range [0, vertex_count).
   */
  void build(const std::vector<size_t>& index_buffer);

  void rebuildLiveIndices();
  void applyFreshHalfEdgeFaceReferences(
    const std::vector<HalfEdge>& fresh_half_edges, const std::vector<int>& tri_remap);
  bool halfEdgeInFaceCycle(size_t he_id, size_t face_id) const;
  bool faceHasInfiniteVertex(size_t face_id) const;
  bool isDeadFaceSlot(const Triangle& tri) const;
  void tombstoneHalfEdge(size_t he_id);
  void killFaceSlot(Triangle& tri, size_t invalid_he);

  struct OuterBoundaryMaps
  {
    std::vector<int> new_outer_outgoing;
    std::vector<int> new_outer_incoming;
    std::vector<int> existing_outer_outgoing;
    std::vector<int> existing_incoming_inf;
  };

  struct CreateInfiniteFacesDebugContext
  {
    const std::function<glm::dvec2(size_t)>* vertex_at = nullptr;
    std::string output_prefix;
  };

  struct AngularBisectorRay
  {
    int apex = -1;
    glm::dvec2 origin {};
    glm::dvec2 direction {};
  };

  OuterBoundaryMaps buildOuterBoundaryMaps(const std::unordered_set<size_t>& new_outer_half_edge_set) const;
  bool infiniteHalfEdgeBordersTombstonedTriangle(size_t he_id) const;
  void clearTombstonedFaceRef(size_t he_id);
  void spliceExistingInfiniteEdgesToNewOuter(
    const OuterBoundaryMaps& maps, std::vector<int>& incoming_inf_override);
  int resolveNextOuterBoundaryHalfEdge(int vertex, const OuterBoundaryMaps& maps) const;
  void createInfiniteFacesFromBoundary(const std::vector<size_t>& outer_half_edges,
    const std::function<glm::dvec2(size_t)>& vertex_at,
    const CreateInfiniteFacesDebugContext* debug_context = nullptr);
  void writeCreateInfiniteFacesPhaseDebugSvg(const std::string& filepath, const std::vector<glm::dvec2>& points,
    const std::string& phase_title, const std::unordered_set<size_t>* highlight_outer_half_edges = nullptr) const;
  void maybeWriteCreateInfiniteFacesPhaseDebug(const CreateInfiniteFacesDebugContext* debug_context,
    const std::string& phase_suffix, const std::string& phase_title,
    const std::unordered_set<size_t>* highlight_outer_half_edges = nullptr) const;
  std::optional<AngularBisectorRay> computeAngularBisectorRayFromMesh(
    size_t he_id, const std::function<const glm::dvec2*(int)>& try_point) const;
  static bool isInfiniteHalfEdgePair(const HalfEdge& he, const HalfEdge& twin);
  static std::string formatHalfEdgeTopologyLabel(
    size_t he_id, const HalfEdge& he, const HalfEdge& twin);
  size_t prevAroundFace(size_t he_id) const;
  size_t nextOnOuterBoundaryRaw(size_t he_id) const;
  size_t prevOnOuterBoundaryRaw(size_t he_id) const;

 public:
  HalfEdgeDelaunayGraph() = default;

  void init(const std::vector<glm::dvec2>& site_positions);

  void init(const std::vector<std::vector<glm::dvec2>>& splines);

  void update(
    size_t vertex_count,
    const std::vector<std::vector<size_t>>& components,
    const std::function<glm::dvec2(size_t)>& vertex_position);

  /**
   * Apply a pending runtime-branch split without retriangulating: tombstone cross-component
   * finite triangles/edges, preserve surviving infinite topology (infinity is a wildcard, not
   * branch-owned), and add infinite faces only along newly opened outer edges.
   *
   * Here "outer" means: this directed half-edge borders an infinite or tombstoned triangle on
   * its own side and a regular finite triangle on its twin side. This is distinct from any
   * alpha-shape boundary classification derived from inside/outside face flags.
   */
  void applyRuntimeBranchSplit(
    const std::vector<size_t>& component_map,
    const std::function<glm::dvec2(size_t)>& vertex_at, std::optional<double> debug_time = std::nullopt);

  /** Tombstone every face/edge incident to a vertex flagged in @p dead_vertex_mask. */
  void killVertices(const std::vector<bool>& dead_vertex_mask);

  // Flips an edge between two triangles by rotating it counter-clockwise in its quadrilateral
  void flipEdge(size_t he_id);
  // Other methods to manipulate and query the triangulation can be added here.
  void printDebug(std::ostream* out = nullptr) const;

  static glm::dvec2 circumcenter(const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c);

  std::vector<std::pair<glm::dvec2, bool>> computeCircumcenters(const std::vector<glm::dvec2>& vertices) const;

  /**
   * Perpendicular ray direction for the Voronoi vertex dual to an infinite Delaunay face.
   * Uses the canonical outside-directed hull half-edge, so the result is invariant to how a
   * triangle stores its three half-edge slots after flips/reordering.
   */
  std::optional<glm::dvec2> infiniteVoronoiRayDirection(
    size_t face_id, const std::function<glm::dvec2(int vertex_index)>& vertex_at) const;

  bool isLiveHalfEdge(size_t he_id) const;
  bool isLiveFace(size_t face_id) const;

  /** @throws std::runtime_error if any live half-edge references a dead or out-of-range face. */
  void validateLiveHalfEdgeFaceReferences() const;

  const HalfEdge& halfEdge(size_t he_id) const;
  const Triangle& face(size_t face_id) const;

  size_t halfEdgeSlotCount() const;
  size_t faceSlotCount() const;
  size_t getVertexCount() const;

  size_t liveDelaunayEdgeCount() const { return live_even_half_edges_.size(); }
  size_t liveFaceCount() const { return live_faces_.size(); }
  size_t liveDelaunayEdgeId(size_t live_index) const { return live_even_half_edges_.at(live_index); }
  size_t liveFaceId(size_t live_index) const { return live_faces_.at(live_index); }

  // utils

  /**
   * Determines whether a half-edge is outside the boundary. This includes all edges from/to infinity as well as the
   * outer half-edges along the boundary.
   */
  bool isOutsideConvexBoundary(size_t he_id) const;

  /**
   * Determines whether a half-edge is on the boundary. This includes both, the inner and the outer half-edges along the
   * boundary.
   */
  bool isOnConvexBoundary(size_t he_id) const;

  /**
   * Determines whether a half-edge is on the outer side along the boundary.
   */
  bool isOnConvexBoundaryOutside(size_t he_id) const;

  /**
   * Determines whether a half-edge is connected to the vertex at infinity.
   */
  bool isInfinite(size_t he_id) const;
  int destination(size_t he_id) const;
  int triangleOppositeVertex(size_t he_id) const;
  std::array<int, 3> adjacentTriangleVertices(size_t he_id) const;
  size_t neighborEdgeId(size_t he_id) const;
  size_t nextOnConvexBoundaryId(size_t he_id) const;

  size_t prevOnConvexBoundaryId(size_t he_id) const;

  std::vector<size_t> neighbors(size_t v);
  std::vector<size_t> inducedNeighbors(size_t v, const std::vector<bool>& face_inside) const;

  /** Neighbors connected by a live half-edge with at least one live incident face. */
  std::vector<size_t> inducedNeighborsFromLiveGraph(size_t v) const;

  /**
   * Returns the half-edge ID of the twin half-edge. The twin half-edge is the half-edge that shares the same edge but has opposite orientation and belongs to the adjacent face.
   */
  static size_t twin(size_t he_id);

  /**
   * Returns the half-edge ID of the half-edge that comes before the given half-edge in its face.
   */
  size_t prev(size_t he_id) const;

  std::array<int, 3> getTriangleVertexIndices(size_t face_id) const;
  std::array<size_t, 3> getTriangleHalfEdgeIndices(size_t start_he_id) const;
  std::array<size_t, 4> getQuadBoundaryHalfEdgeIndices(size_t quad_id) const;

  void reorder_from_old(const std::vector<Triangle>& old_triangles, const std::vector<HalfEdge>& old_half_edges);

  // ---------------- Live Delaunay edge iterator ----------------
  class LiveDelaunayEdgeIterator
  {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const size_t*;
    using reference = const size_t&;

    LiveDelaunayEdgeIterator(const std::vector<size_t>* live, size_t idx)
      : live_(live)
      , idx_(idx)
    {
    }

    value_type operator*() const { return (*live_)[idx_]; }

    LiveDelaunayEdgeIterator& operator++()
    {
      ++idx_;
      return *this;
    }

    LiveDelaunayEdgeIterator operator++(int)
    {
      LiveDelaunayEdgeIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const LiveDelaunayEdgeIterator& other) const { return live_ == other.live_ && idx_ == other.idx_; }

    bool operator!=(const LiveDelaunayEdgeIterator& other) const { return !(*this == other); }

   private:
    const std::vector<size_t>* live_ = nullptr;
    size_t idx_ = 0;
  };

  class LiveDelaunayEdgeRange
  {
   public:
    explicit LiveDelaunayEdgeRange(const HalfEdgeDelaunayGraph* graph)
      : graph_(graph)
    {
    }

    LiveDelaunayEdgeIterator begin() const
    {
      return LiveDelaunayEdgeIterator(&graph_->live_even_half_edges_, 0);
    }

    LiveDelaunayEdgeIterator end() const
    {
      return LiveDelaunayEdgeIterator(&graph_->live_even_half_edges_, graph_->live_even_half_edges_.size());
    }

   private:
    const HalfEdgeDelaunayGraph* graph_ = nullptr;
  };

  LiveDelaunayEdgeRange liveDelaunayEdges() const { return LiveDelaunayEdgeRange(this); }

  // ---------------- Live face iterator ----------------
  class LiveFaceIterator
  {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const size_t*;
    using reference = const size_t&;

    LiveFaceIterator(const std::vector<size_t>* live, size_t idx)
      : live_(live)
      , idx_(idx)
    {
    }

    value_type operator*() const { return (*live_)[idx_]; }

    LiveFaceIterator& operator++()
    {
      ++idx_;
      return *this;
    }

    LiveFaceIterator operator++(int)
    {
      LiveFaceIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const LiveFaceIterator& other) const { return live_ == other.live_ && idx_ == other.idx_; }

    bool operator!=(const LiveFaceIterator& other) const { return !(*this == other); }

   private:
    const std::vector<size_t>* live_ = nullptr;
    size_t idx_ = 0;
  };

  class LiveFaceRange
  {
   public:
    explicit LiveFaceRange(const HalfEdgeDelaunayGraph* graph)
      : graph_(graph)
    {
    }

    LiveFaceIterator begin() const { return LiveFaceIterator(&graph_->live_faces_, 0); }

    LiveFaceIterator end() const { return LiveFaceIterator(&graph_->live_faces_, graph_->live_faces_.size()); }

   private:
    const HalfEdgeDelaunayGraph* graph_ = nullptr;
  };

  LiveFaceRange liveFaces() const { return LiveFaceRange(this); }

  // ---------------- Incident edge iterator ----------------
  class IncidentEdgeIterator
  {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t; // returning halfedge indices
    using difference_type = std::ptrdiff_t;
    using pointer = const size_t*;
    using reference = const size_t&;

    IncidentEdgeIterator(const HalfEdgeDelaunayGraph* g, size_t v, bool end = false)
      : g_(g)
      , v_(v)
      , start_he_(g ? g->vertex_to_half_edge[v] : npos)
      , curr_he_(end ? npos : start_he_)
      , first_(true)
    {
      if (!end && g_ && curr_he_ != npos && !g_->isLiveHalfEdge(curr_he_))
      {
        advance();
      }
    }

    value_type operator*() const { return curr_he_; }

    IncidentEdgeIterator& operator++()
    {
      advance();
      return *this;
    }

    IncidentEdgeIterator operator++(int)
    {
      IncidentEdgeIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const IncidentEdgeIterator& other) const
    {
      return curr_he_ == other.curr_he_ && v_ == other.v_ && g_ == other.g_;
    }

    bool operator!=(const IncidentEdgeIterator& other) const { return !(*this == other); }

   private:
    void advance()
    {
      if (curr_he_ == npos)
      {
        return;
      }

      do
      {
        curr_he_ = g_->neighborEdgeId(curr_he_);
        if (curr_he_ == start_he_ && !first_)
        {
          curr_he_ = npos;
          return;
        }
        first_ = false;
      } while (!g_->isLiveHalfEdge(curr_he_));
    }

    const HalfEdgeDelaunayGraph* g_;
    size_t v_;
    size_t start_he_;
    size_t curr_he_;
    bool first_;
    static constexpr size_t npos = static_cast<size_t>(-1);
  };

  IncidentEdgeIterator incidentEdgesBegin(size_t v) const { return IncidentEdgeIterator(this, v, false); }

  IncidentEdgeIterator incidentEdgesEnd(size_t v) const { return IncidentEdgeIterator(this, v, true); }

  // ---------------- Convex hull edge iterator ----------------
  class ConvexHullEdgeIterator
  {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t; // returning halfedge indices
    using difference_type = std::ptrdiff_t;
    using pointer = const size_t*;
    using reference = const size_t&;

    ConvexHullEdgeIterator(const HalfEdgeDelaunayGraph* g, size_t he_id, bool end = false)
      : g_(g)
      , start_he_(he_id)
      , curr_he_(end ? npos : he_id)
      , first_(true)
    {
      if (!end)
      {
        assert(g_->isLiveHalfEdge(curr_he_) && "Iterator started on tombstone half-edge!");
        assert(g_->triangleOppositeVertex(curr_he_) == -1 && "Iterator started on non-boundary half-edge!");
      }
    }

    value_type operator*() const { return curr_he_; }

    ConvexHullEdgeIterator& operator++()
    {
      if (curr_he_ == npos)
      {
        return *this;
      }

      do
      {
        curr_he_ = g_->nextOnConvexBoundaryId(curr_he_);
        assert(g_->triangleOppositeVertex(curr_he_) == -1 && "Iterator moved to non-boundary half-edge!");
        if (curr_he_ == start_he_ && !first_)
        {
          curr_he_ = npos;
          return *this;
        }
        first_ = false;
      } while (!g_->isLiveHalfEdge(curr_he_));

      return *this;
    }

    ConvexHullEdgeIterator operator++(int)
    {
      ConvexHullEdgeIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const ConvexHullEdgeIterator& other) const { return curr_he_ == other.curr_he_ && g_ == other.g_; }

    bool operator!=(const ConvexHullEdgeIterator& other) const { return !(*this == other); }

   private:
    const HalfEdgeDelaunayGraph* g_;
    size_t start_he_;
    size_t curr_he_;
    bool first_;
    static constexpr size_t npos = static_cast<size_t>(-1);
  };

  ConvexHullEdgeIterator boundaryEdgesBegin() const
  {
    for (size_t he_id : live_even_half_edges_)
    {
      if (isOnConvexBoundaryOutside(he_id))
      {
        return ConvexHullEdgeIterator(this, he_id, false);
      }
    }
    return boundaryEdgesEnd();
  }

  ConvexHullEdgeIterator boundaryEdgesEnd() const { return ConvexHullEdgeIterator(this, 0, true); }
};
} // namespace kinDS
