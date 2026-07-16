#pragma once

#include "KineticDelaunay.hpp"

#include <cassert>
#include <glm/glm.hpp>
#include <list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace kinDS
{
struct KineticDelaunay::CrossingData
{
  struct VoronoiDelaunayEdgeIntersection;
  typedef std::list<VoronoiDelaunayEdgeIntersection>::iterator EdgeIntersectionRef;

  // Per-edge intersection ref lists are heap-allocated so vector growth does not move list objects and
  // invalidate cached iterators (delaunay_ref / voronoi_ref) on MSVC debug builds.
  struct EdgeIntersectionRefListSlots
  {
    using RefList = std::list<EdgeIntersectionRef>;

    RefList& operator[](size_t edge_id) { return *slots_.at(edge_id); }
    const RefList& operator[](size_t edge_id) const { return *slots_.at(edge_id); }

    size_t size() const { return slots_.size(); }

    void growTo(size_t count)
    {
      if (count <= slots_.size())
      {
        return;
      }

      const size_t old_count = slots_.size();
      slots_.resize(count);
      for (size_t edge_id = old_count; edge_id < count; ++edge_id)
      {
        slots_[edge_id] = std::make_unique<RefList>();
      }
    }

    void resize(size_t count)
    {
      const size_t old_count = slots_.size();
      slots_.resize(count);
      for (size_t edge_id = old_count; edge_id < count; ++edge_id)
      {
        slots_[edge_id] = std::make_unique<RefList>();
      }
    }

    void clearAllLists()
    {
      for (const auto& slot : slots_)
      {
        if (slot)
        {
          slot->clear();
        }
      }
    }

   private:
    std::vector<std::unique_ptr<RefList>> slots_;
  };

  static constexpr size_t invalid_containing_tri_id = static_cast<size_t>(-1);

 private:
  using VoronoiVertexList = std::list<size_t>;
  using VoronoiVertexListIterator = VoronoiVertexList::iterator;

  std::vector<size_t> voronoi_vertex_to_containing_tri_id;
  std::vector<std::unique_ptr<VoronoiVertexList>> tri_id_to_voronoi_vertices;
  std::vector<std::optional<VoronoiVertexListIterator>> voronoi_vertex_to_iterator;

  void ensureTriListSlots(size_t face_count)
  {
    const size_t old_count = tri_id_to_voronoi_vertices.size();
    tri_id_to_voronoi_vertices.resize(face_count);
    for (size_t tri_id = old_count; tri_id < face_count; ++tri_id)
    {
      tri_id_to_voronoi_vertices[tri_id] = std::make_unique<VoronoiVertexList>();
    }
  }

  VoronoiVertexList& triList(size_t tri_id)
  {
    return *tri_id_to_voronoi_vertices.at(tri_id);
  }

  const VoronoiVertexList& triList(size_t tri_id) const
  {
    return *tri_id_to_voronoi_vertices.at(tri_id);
  }

 public:
  void assertVoronoiVertexIteratorMatchesId(size_t voronoi_vertex_id) const
  {
    assert(voronoi_vertex_id < voronoi_vertex_to_iterator.size());
    assert(voronoi_vertex_to_iterator[voronoi_vertex_id].has_value());
    assert(**voronoi_vertex_to_iterator[voronoi_vertex_id] == voronoi_vertex_id);
  }

  std::list<VoronoiDelaunayEdgeIntersection> edge_intersections;
  EdgeIntersectionRefListSlots voronoi_edge_intersections;
  EdgeIntersectionRefListSlots delaunay_edge_intersections;

  std::vector<double> last_crossing;

  struct VoronoiDelaunayEdgeIntersection
  {
    std::optional<std::list<EdgeIntersectionRef>::iterator> delaunay_ref;
    std::optional<std::list<EdgeIntersectionRef>::iterator> voronoi_ref;
    size_t delaunay_edge_id;
    size_t voronoi_edge_id;
    double delaunay_edge_param;
    // SegmentBuilder boundary-interval mesh linkage along one Delaunay edge.
    size_t prev_segment_mesh_pair_index = static_cast<size_t>(-1);
    size_t next_segment_mesh_pair_index = static_cast<size_t>(-1);
  };

  void growTo(size_t face_count)
  {
    if (face_count <= voronoi_vertex_to_containing_tri_id.size())
    {
      return;
    }

    voronoi_vertex_to_containing_tri_id.resize(face_count, invalid_containing_tri_id);
    ensureTriListSlots(face_count);
    voronoi_vertex_to_iterator.resize(face_count);
    last_crossing.resize(face_count, 0.0);
  }

  void growEdgeSlotsTo(size_t delaunay_edge_count)
  {
    voronoi_edge_intersections.growTo(delaunay_edge_count);
    delaunay_edge_intersections.growTo(delaunay_edge_count);
  }

  void init(size_t face_count)
  {
    voronoi_vertex_to_containing_tri_id.clear();
    voronoi_vertex_to_containing_tri_id.resize(face_count, invalid_containing_tri_id);

    tri_id_to_voronoi_vertices.clear();
    ensureTriListSlots(face_count);

    voronoi_vertex_to_iterator.clear();
    voronoi_vertex_to_iterator.resize(face_count);

    last_crossing.clear();
    last_crossing.resize(face_count, 0.0);
  }

  bool isVoronoiVertexRegistered(size_t voronoi_vertex_id) const
  {
    return voronoi_vertex_id < voronoi_vertex_to_containing_tri_id.size()
      && voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] != invalid_containing_tri_id
      && voronoi_vertex_to_iterator[voronoi_vertex_id].has_value();
  }

  std::optional<size_t> peekContainingTriId(size_t voronoi_vertex_id) const
  {
    if (!isVoronoiVertexRegistered(voronoi_vertex_id))
    {
      return std::nullopt;
    }
    return voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
  }

  void setVoronoiVertexTriId(size_t voronoi_vertex_id, size_t tri_id)
  {
    if (voronoi_vertex_id >= voronoi_vertex_to_containing_tri_id.size() || tri_id >= tri_id_to_voronoi_vertices.size())
    {
      throw std::runtime_error("setVoronoiVertexTriId: voronoi vertex or triangle index out of range");
    }

    if (voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] != invalid_containing_tri_id
      && !voronoi_vertex_to_iterator[voronoi_vertex_id].has_value())
    {
      throw std::runtime_error("setVoronoiVertexTriId: Voronoi vertex " + std::to_string(voronoi_vertex_id)
        + " has containing triangle " + std::to_string(voronoi_vertex_to_containing_tri_id[voronoi_vertex_id])
        + " but no valid cached list iterator");
    }

    if (isVoronoiVertexRegistered(voronoi_vertex_id))
    {
      const size_t current_tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
      if (current_tri_id == tri_id)
      {
        assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
        return;
      }

      assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
      triList(current_tri_id).erase(*voronoi_vertex_to_iterator[voronoi_vertex_id]);
      voronoi_vertex_to_iterator[voronoi_vertex_id].reset();
      voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] = invalid_containing_tri_id;
    }

    voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] = tri_id;
    voronoi_vertex_to_iterator[voronoi_vertex_id]
      = triList(tri_id).emplace(triList(tri_id).end(), voronoi_vertex_id);
    assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
  }

  void unsetVoronoiVertex(size_t voronoi_vertex_id)
  {
    if (voronoi_vertex_id >= voronoi_vertex_to_containing_tri_id.size())
    {
      return;
    }

    if (!isVoronoiVertexRegistered(voronoi_vertex_id))
    {
      assert(!voronoi_vertex_to_iterator[voronoi_vertex_id].has_value());
      return;
    }

    assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
    const size_t tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
    triList(tri_id).erase(*voronoi_vertex_to_iterator[voronoi_vertex_id]);
    voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] = invalid_containing_tri_id;
    voronoi_vertex_to_iterator[voronoi_vertex_id].reset();
  }

  void moveVertex(size_t voronoi_vertex_id, size_t target_tri_id, double t)
  {
    if (voronoi_vertex_id >= voronoi_vertex_to_containing_tri_id.size() || target_tri_id >= tri_id_to_voronoi_vertices.size())
    {
      throw std::runtime_error("moveVertex: voronoi vertex or target triangle index out of range");
    }

    const size_t current_tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
    if (current_tri_id == target_tri_id)
    {
      if (!isVoronoiVertexRegistered(voronoi_vertex_id))
      {
        setVoronoiVertexTriId(voronoi_vertex_id, target_tri_id);
      }
      else
      {
        assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
      }
      return;
    }

    if (!isVoronoiVertexRegistered(voronoi_vertex_id))
    {
      setVoronoiVertexTriId(voronoi_vertex_id, target_tri_id);
      return;
    }

    assertVoronoiVertexIteratorMatchesId(voronoi_vertex_id);
    triList(current_tri_id).erase(*voronoi_vertex_to_iterator[voronoi_vertex_id]);
    voronoi_vertex_to_iterator[voronoi_vertex_id].reset();
    voronoi_vertex_to_containing_tri_id[voronoi_vertex_id] = invalid_containing_tri_id;
    setVoronoiVertexTriId(voronoi_vertex_id, target_tri_id);

    KINDS_DEBUG("Voronoi vertex " << voronoi_vertex_id << " moved from triangle " << current_tri_id << " to "
                                  << target_tri_id << " at t = " << t);
  }

  size_t getContainingTriId(size_t voronoi_vertex_id) const
  {
    if (voronoi_vertex_id >= voronoi_vertex_to_containing_tri_id.size())
    {
      throw std::runtime_error("getContainingTriId: voronoi vertex index out of range");
    }

    if (!isVoronoiVertexRegistered(voronoi_vertex_id))
    {
      const size_t containing_tri_id = voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
      std::string msg = "getContainingTriId: Voronoi vertex " + std::to_string(voronoi_vertex_id) + " is not registered";
      if (containing_tri_id != invalid_containing_tri_id)
      {
        msg += " (stale containing_tri_id=" + std::to_string(containing_tri_id) + " without list iterator)";
      }
      throw std::runtime_error(msg);
    }
    return voronoi_vertex_to_containing_tri_id[voronoi_vertex_id];
  }

  const std::vector<size_t>& getContainingTriIds() const { return voronoi_vertex_to_containing_tri_id; }

  // Note: we copy this into a vector because we need this for reassigning Voronoi vertices in a quadrilateral and we
  // will be modifying the underlying list while iterating
  std::vector<size_t> getVoronoiVerticesInTri(size_t tri_id) const
  {
    const VoronoiVertexList& vertices = triList(tri_id);
    return std::vector<size_t>(vertices.begin(), vertices.end());
  }

  /**
   * Throws if cached iterators, containing-triangle ids, or per-triangle lists disagree.
   * When @p kd is non-null, may write a debug SVG before throwing.
   */
  void validateVoronoiVertexIteratorInvariants(
    const char* context, const HalfEdgeDelaunayGraph& graph, const KineticDelaunay* kd = nullptr, double t = 0.0) const;

  void computeEdgeIntersections(const KineticDelaunay& kd, double t);

  // Update Voronoi–Delaunay edge intersections after a single crossing event.
  void updateAfterCrossingEvent(const KineticDelaunay& kd, const KineticDelaunay::CrossingEvent& e);

  // Remove a single intersection from all three data structures (global list,
  // per-Voronoi-edge list, and per-Delaunay-edge list).
  void removeIntersection(EdgeIntersectionRef intersection_ref);

  // Remove all intersections whose Delaunay edge no longer exists in the current graph.
  void removeIntersectionsOnDeadDelaunayEdges(const HalfEdgeDelaunayGraph& graph);

  /**
   * Throws if per-edge lists, cached refs, or Delaunay-edge param ordering are inconsistent.
   * When @p kd is non-null, writes a highlighted debug SVG before throwing.
   */
  void validateIntersectionInvariants(const char* context, const KineticDelaunay* kd = nullptr, double t = 0.0) const;
};

class KineticDelaunay::CrossingEvent final : public KineticDelaunay::Event
{
 public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;
  size_t voronoi_vertex_id;

  CrossingEvent(
    KineticDelaunay* kd, double t, size_t he_id, double creation_time, glm::dvec2 position, size_t voronoi_vertex_id)
    : KineticDelaunay::Event(t, creation_time, 40u)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
    , voronoi_vertex_id(voronoi_vertex_id)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::CrossingEventManager final : public KineticDelaunay::EventManager
{
 public:
  explicit CrossingEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

  // Expose crossing data to the KineticDelaunay orchestration code and to events.
  CrossingData& getCrossingDataMutable() { return crossing_data_; }
  const CrossingData& getCrossingData() const { return crossing_data_; }

 private:
  KineticDelaunay* kd_;
  CrossingData crossing_data_;
};

inline void KineticDelaunay::CrossingEvent::handleEvent()
{
  auto* kd = kd_;
  if (!kd)
  {
    throw std::runtime_error("CrossingEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;

  // Check if the event is still valid
  if (creation_time < kd->crossing_data.last_crossing[voronoi_vertex_id])
  {
    return;
  }

  if (voronoi_vertex_id >= graph.faceSlotCount() || !graph.isLiveFace(voronoi_vertex_id)
    || !kd->crossing_data.isVoronoiVertexRegistered(voronoi_vertex_id))
  {
    return;
  }

  size_t containing_tri_id = kd->crossing_data.getContainingTriId(voronoi_vertex_id);

  // Outdated if the containing triangle or the dual triangle (same index as the Voronoi vertex) was
  // updated by a flip after this event was scheduled.
  if (creation_time < kd->face_last_updated[containing_tri_id])
  {
    return;
  }
  if (creation_time < kd->face_last_updated[voronoi_vertex_id])
  {
    return;
  }

  auto* event_handler = kd->crossing_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  kd->crossing_data.last_crossing[voronoi_vertex_id] = occurrence_time;

  KINDS_DEBUG("Processing crossing event at time " << occurrence_time << " for Voronoi vertex ID " << voronoi_vertex_id
                                                   << " crossing half-edge ID " << half_edge_id);

  // move to neighboring triangle
  KINDS_DEBUG("Moving Voronoi vertex " << voronoi_vertex_id << " from triangle " << containing_tri_id << " to triangle "
                                       << graph.halfEdge(half_edge_id ^ 1).face);
  kd->crossing_data.moveVertex(voronoi_vertex_id, graph.halfEdge(half_edge_id ^ 1).face, occurrence_time);

  // Update Voronoi–Delaunay edge intersections stored in crossing_data in response to this crossing.
  kd->crossing_data.updateAfterCrossingEvent(*kd, *this);

  // Recompute params on all three Delaunay edges of this Voronoi vertex (Delaunay face) before callbacks mesh.
  kd->refreshTriangleDelaunayEdgeIntersectionParams(voronoi_vertex_id, occurrence_time);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  // After callbacks (e.g. debug SVG export); intersection lists must be consistent.
  kd->validateVoronoiVertexIteratorInvariants("CrossingEvent:afterEvent", occurrence_time);
  kd->validateCrossingIntersectionInvariants("CrossingEvent:afterEvent", occurrence_time);

  // Re-compute crossing events for this Voronoi vertex
  kd->crossing_event_manager_->computeEvents(occurrence_time, voronoi_vertex_id);

  if (kd->diagnosticsEnabled() && kd->isDiagnosticsMonitoredFaceValid())
  {
    kd->logFaceInsideStateAtTime(kDiagnosticsMonitoredFaceId, occurrence_time, "crossing_event");
  }
}

/** Empty optional formats as "null". Log form is V{idx}xD{idx} for list positions along voronoi/delaunay edge lists. */
std::string formatCrossingIntersectionForLog(const KineticDelaunay& kd,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection);

bool tryComputeCrossingIntersectionPosition2D(const KineticDelaunay& kd,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, double t, glm::dvec2& out_xy,
  bool apply_reference_transform = true, bool include_virtual_offset = true);

/** No-op if @p intersection is empty; otherwise log KINDS_ERROR on mismatch with expected dual edge / half-edge. */
void validateClosingCapCrossingRef(const KineticDelaunay& kd, const char* context_msg,
  std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef> intersection, size_t expected_voronoi_edge_id,
  int delaunay_half_edge_id);

} // namespace kinDS
