#include "SegmentBuilderSubdivisionCallback.hpp"

#include "SegmentBuilderVisualDebug.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunaySubdivisionEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilder.hpp"

#include <algorithm>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace kinDS
{
void SegmentBuilderSubdivisionCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* sub = dynamic_cast<KineticDelaunay::SubdivisionEvent*>(&e);
  if (!sub)
  {
    return;
  }

  const size_t strand_id = sub->strand_id;
  const double t = sub->occurrence_time;

  if (segment_builder_.visual_debug)
  {
    auto& debug_graph = segment_builder_.kin_del.getGraph();
    writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, debug_graph, t, "before",
      "subdivision_strand" + std::to_string(strand_id) + "_seq" + std::to_string(sub->queue_sequence_),
      VisualDebugHighlight::forSubdivisionStrand(debug_graph, strand_id));
  }

  KINDS_DEBUG("Inserting subdivision for strand " << strand_id << " at t = " << t);
  //   Traverse all half-edges around this strand and insert a new vertex into the corresponding segment meshes
  auto& graph = segment_builder_.kin_del.getGraph();

  // compute boundary polygon for component at time t
  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
  size_t component_id = segment_builder_.kin_del.component_data.component_map[strand_id];
  segment_builder_.updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = polygonCentroid(boundary_polygon);
  std::vector<size_t> cell_boundary_he_ids;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    cell_boundary_he_ids.push_back(*it);
  }

  // The authoritative boundary-interval list comes from createClosingMesh() tracing.
  // We use that exact list for finish/start calls.
  auto format_intersection_ref
    = [](const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& ref) -> std::string
  {
    if (!ref.has_value())
    {
      return "null";
    }
    std::ostringstream oss;
    oss << "{de=" << ref.value()->delaunay_edge_id << ", ve=" << ref.value()->voronoi_edge_id
        << ", de_param=" << ref.value()->delaunay_edge_param << "}";
    return oss.str();
  };

  auto interval_delaunay_edge
    = [](const SegmentBuilder::BoundaryIntersectionInterval& interval) -> std::optional<size_t>
  {
    if (interval.start_intersection.has_value())
    {
      return interval.start_intersection.value()->delaunay_edge_id;
    }
    if (interval.end_intersection.has_value())
    {
      return interval.end_intersection.value()->delaunay_edge_id;
    }
    return std::nullopt;
  };

  auto format_refs_vector = [&](const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef>& refs) -> std::string
  {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < refs.size(); ++i)
    {
      if (i > 0)
      {
        oss << ", ";
      }
      oss << format_intersection_ref(std::make_optional(refs[i]));
    }
    oss << "]";
    return oss.str();
  };

  std::unordered_map<size_t, std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef>> boundary_refs_cache;
  auto boundary_refs_for_edge = [&](size_t d_edge_id)
    -> const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef>&
  {
    auto it = boundary_refs_cache.find(d_edge_id);
    if (it != boundary_refs_cache.end())
    {
      return it->second;
    }
    const auto inserted
      = boundary_refs_cache.emplace(d_edge_id, segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_id));
    return inserted.first->second;
  };

  auto find_ref_index = [](const std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef>& refs,
                         KineticDelaunay::CrossingData::EdgeIntersectionRef needle) -> size_t
  {
    for (size_t i = 0; i < refs.size(); ++i)
    {
      if (refs[i] == needle)
      {
        return i;
      }
    }
    return static_cast<size_t>(-1);
  };

  auto validate_traced_interval = [&](const SegmentBuilder::BoundaryIntersectionInterval& interval, size_t interval_idx)
  {
    const auto d_edge_opt = interval_delaunay_edge(interval);
    if (!d_edge_opt.has_value())
    {
      std::ostringstream oss;
      oss << "Subdivision traced interval has no intersection endpoints (idx=" << interval_idx << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }
    const size_t d_edge_id = d_edge_opt.value();
    const auto& refs = boundary_refs_for_edge(d_edge_id);
    const size_t count = refs.size();
    if (count == 0)
    {
      std::ostringstream oss;
      oss << "Subdivision traced interval references d_edge " << d_edge_id
          << " but boundary-order intersection list is empty.";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    const bool has_start = interval.start_intersection.has_value();
    const bool has_end = interval.end_intersection.has_value();
    if (!has_start && !has_end)
    {
      std::ostringstream oss;
      oss << "Subdivision traced interval has null start and end (idx=" << interval_idx << ", d_edge=" << d_edge_id
          << ", refs=" << format_refs_vector(refs) << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    if (!has_start && has_end)
    {
      if (refs.front() != interval.end_intersection.value())
      {
        std::ostringstream oss;
        oss << "Subdivision traced interval [null,ref] is not boundary-first (idx=" << interval_idx
            << ", d_edge=" << d_edge_id << ", end=" << format_intersection_ref(interval.end_intersection)
            << ", expected_first=" << format_intersection_ref(std::make_optional(refs.front()))
            << ", refs=" << format_refs_vector(refs) << ").";
        KINDS_ERROR(oss.str());
        throw std::runtime_error(oss.str());
      }
      return;
    }

    if (has_start && !has_end)
    {
      if (refs.back() != interval.start_intersection.value())
      {
        std::ostringstream oss;
        oss << "Subdivision traced interval [ref,null] is not boundary-last (idx=" << interval_idx
            << ", d_edge=" << d_edge_id << ", start=" << format_intersection_ref(interval.start_intersection)
            << ", expected_last=" << format_intersection_ref(std::make_optional(refs.back()))
            << ", refs=" << format_refs_vector(refs) << ").";
        KINDS_ERROR(oss.str());
        //throw std::runtime_error(oss.str());
      }
      return;
    }

    const size_t start_idx = find_ref_index(refs, interval.start_intersection.value());
    const size_t end_idx = find_ref_index(refs, interval.end_intersection.value());
    if (start_idx == static_cast<size_t>(-1) || end_idx == static_cast<size_t>(-1))
    {
      std::ostringstream oss;
      oss << "Subdivision traced interval [ref,ref] endpoints not found on boundary-ordered list (idx=" << interval_idx
          << ", d_edge=" << d_edge_id << ", start=" << format_intersection_ref(interval.start_intersection)
          << ", end=" << format_intersection_ref(interval.end_intersection) << ", refs=" << format_refs_vector(refs)
          << ").";
      KINDS_ERROR(oss.str());
      throw std::runtime_error(oss.str());
    }

    if (end_idx == start_idx + 1)
    {
      return;
    }

    // Validation-only normalization for edges with >1 crossings:
    // allow reversed adjacent pairs to be diagnosed as swapped-order intervals.
    if (count > 1 && start_idx == end_idx + 1)
    {
      KINDS_WARNING("Subdivision traced interval appears reversed (normalizable) on d_edge "
                    << d_edge_id << " at idx " << interval_idx << " start_idx=" << start_idx
                    << " end_idx=" << end_idx << ". Caller should pass boundary-direction order.");
      return;
    }

    std::ostringstream oss;
    oss << "Subdivision traced interval [ref,ref] is not adjacent in boundary order (idx=" << interval_idx
        << ", d_edge=" << d_edge_id << ", start_idx=" << start_idx << ", end_idx=" << end_idx
        << ", start=" << format_intersection_ref(interval.start_intersection)
        << ", end=" << format_intersection_ref(interval.end_intersection) << ", refs=" << format_refs_vector(refs)
        << ").";
    KINDS_ERROR(oss.str());
    throw std::runtime_error(oss.str());
  };

  auto extend_intersection_neighbor_at_shared_crossing
    = [&](size_t neighbor_pair_idx, KineticDelaunay::CrossingData::EdgeIntersectionRef shared_ref,
        bool update_start_on_neighbor, size_t cur_pair_idx)
  {
    segment_builder_.extendIntersectionMeshAtSharedCrossing(neighbor_pair_idx, shared_ref, update_start_on_neighbor, t,
      SegmentBuilder::BoundaryEventType::Subdivision, SegmentBuilder::BoundarySegmentAction::SegmentRemapped, nullptr, true,
      cur_pair_idx);
  };

  // finish old meshes
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    segment_builder_.finishMesh(*it, t, boundary_polygon, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
  }

  size_t new_segment_id = segment_builder_.segment_properties.size();

  // Create closing mesh and retrieve the exact Delaunay-boundary intervals traced for this strand.
  std::vector<SegmentBuilder::BoundaryIntersectionInterval> traced_boundary_intervals;
  size_t closing_mesh_index
    = segment_builder_.createClosingMesh(strand_id, t, boundary_polygon, centroid, &traced_boundary_intervals);

  // Finish boundary-interval meshes using the interval list produced during closing-mesh tracing.
  for (size_t interval_idx = 0; interval_idx < traced_boundary_intervals.size(); ++interval_idx)
  {
    const auto& interval = traced_boundary_intervals[interval_idx];
    validate_traced_interval(interval, interval_idx);
    const size_t mesh_id = segment_builder_.resolveIntersectionMeshPairIndex(
      interval.voronoi_cell_id, interval.start_intersection, interval.end_intersection, t);

    const auto d_edge_opt_finish = interval_delaunay_edge(interval);
    if (d_edge_opt_finish.has_value())
    {
      const auto& refs_finish = boundary_refs_for_edge(d_edge_opt_finish.value());
      if (!refs_finish.empty())
      {
        if (interval.start_intersection.has_value())
        {
          extend_intersection_neighbor_at_shared_crossing(
            interval.start_intersection.value()->prev_segment_mesh_pair_index, interval.start_intersection.value(), false,
            mesh_id);
        }
        else
        {
          extend_intersection_neighbor_at_shared_crossing(
            refs_finish.front()->prev_segment_mesh_pair_index, refs_finish.front(), false, mesh_id);
        }

        if (interval.end_intersection.has_value())
        {
          extend_intersection_neighbor_at_shared_crossing(
            interval.end_intersection.value()->next_segment_mesh_pair_index, interval.end_intersection.value(), true,
            mesh_id);
        }
        else
        {
          extend_intersection_neighbor_at_shared_crossing(
            refs_finish.back()->next_segment_mesh_pair_index, refs_finish.back(), true, mesh_id);
        }
      }
    }

    KINDS_DEBUG("Subdivision boundary interval finish: strand/cell=" << interval.voronoi_cell_id << " t=" << t
      << " interval_idx=" << interval_idx << " mesh_id=" << mesh_id
      << " start_ref=" << format_intersection_ref(interval.start_intersection)
      << " end_ref=" << format_intersection_ref(interval.end_intersection));
    segment_builder_.finishMeshFromIntersections(interval.voronoi_cell_id, t, interval.start_intersection,
      interval.end_intersection, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
  }

  MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_builder_.segment_mesh_pairs[closing_mesh_index];
  segment_mesh_pair.segment_index0 = segment_builder_.strand_to_segment_indices[strand_id].back();
  segment_mesh_pair.segment_index1 = new_segment_id;

  // Create a new segment mesh property for the new segment
  MeshStructure::SegmentProperties properties;
  segment_builder_.segment_properties.push_back(properties);
  segment_builder_.strand_to_segment_indices[strand_id].push_back(new_segment_id);

  // Start new meshes
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    segment_builder_.startNewMesh(*it, t, false, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::NewSegment);

    // insert vertices into adjacent meshes
    auto& he = graph.halfEdge(*it);

    size_t adjacent_he_id = he.next;
    size_t voronoi_vertex_id = graph.halfEdge(adjacent_he_id).face;
    size_t containing_face_id = segment_builder_.kin_del.getCrossingDataContainingTriId(voronoi_vertex_id);
    bool voronoi_vertex_inside = segment_builder_.kin_del.getFacesInside()[containing_face_id];

    // If the vertex is outside, the mesh is actually not a neighbor and we can skip it.
    if (!voronoi_vertex_inside)
    {
      continue;
    }

    size_t adjacent_segment_mesh_pair_index
      = segment_builder_.half_edge_index_to_segment_mesh_pair_index[adjacent_he_id];
    VoronoiMesh& adjacent_mesh = segment_builder_.meshes[adjacent_segment_mesh_pair_index];

    const size_t adj_even = adjacent_he_id & ~size_t(1);
    const auto& adj_he = graph.halfEdge(adj_even);
    const auto& adj_twin = graph.halfEdge(adj_even ^ size_t(1));
    const std::string adj_vertex_meta = segment_builder_.composeRegularStripVertexMetadata(t, adj_even / 2, adj_even,
      static_cast<int>(adj_he.origin), static_cast<int>(adj_twin.origin), SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentRemapped, std::nullopt, "cross", "subdivision_adjacent_incident");
    const std::string adj_face_meta = segment_builder_.composeRegularStripFaceMetadata(t, adj_even / 2, adj_even,
      static_cast<int>(adj_he.origin), static_cast<int>(adj_twin.origin), SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentRemapped, "subdivision_adjacent_extend");

    glm::dvec3 vertex = segment_builder_.computeVoronoiVertex(adjacent_he_id, t);
    size_t new_vertex_index = segment_builder_.addMeshletVertex(
      adjacent_mesh, boundary_polygon, centroid, vertex, strand_id, t, std::optional<size_t>(voronoi_vertex_id), adj_vertex_meta);
    auto& segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[adjacent_segment_mesh_pair_index];

    if (!segments.empty())
    {
      int he_id_left = segments.front().start_half_edge_id;
      int he_id_right = segments.back().end_half_edge_id;

      if (adjacent_he_id % 2 == 0 && he_id_left == -1)
      {
        size_t last_left = segments.front().mesh_start_vertex_id;
        size_t last_right = segments.front().mesh_end_vertex_id;
        {
          const size_t tris_before = adjacent_mesh.getTriangleCount();
          segment_builder_.addMeshletTriangle(adjacent_mesh, last_left, last_right, new_vertex_index, adj_face_meta);
          if (segment_builder_.diagnostics)
          {
            std::ostringstream note;
            note << "extend_subdivision_adjacent d_tris=" << (adjacent_mesh.getTriangleCount() - tris_before);
            segment_builder_.meshletDiagnosticLogLine("extend_mesh", adjacent_he_id, t, note.str().c_str());
          }
        }
        segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
      }
      else if (adjacent_he_id % 2 != 0 && he_id_right == -1)
      {
        size_t last_left = segments.back().mesh_start_vertex_id;
        size_t last_right = segments.back().mesh_end_vertex_id;
        {
          const size_t tris_before = adjacent_mesh.getTriangleCount();
          segment_builder_.addMeshletTriangle(adjacent_mesh, last_left, last_right, new_vertex_index, adj_face_meta);
          if (segment_builder_.diagnostics)
          {
            std::ostringstream note;
            note << "extend_subdivision_adjacent d_tris=" << (adjacent_mesh.getTriangleCount() - tris_before);
            segment_builder_.meshletDiagnosticLogLine("extend_mesh", adjacent_he_id, t, note.str().c_str());
          }
        }
        segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
      }
    }
  }

  // Restart exactly the same traced intervals so the event-time boundary handoff is consistent.
  for (size_t interval_idx = 0; interval_idx < traced_boundary_intervals.size(); ++interval_idx)
  {
    const auto& interval = traced_boundary_intervals[interval_idx];
    validate_traced_interval(interval, interval_idx);
    const size_t mesh_id_before = segment_builder_.resolveIntersectionMeshPairIndex(
      interval.voronoi_cell_id, interval.start_intersection, interval.end_intersection, t);
    KINDS_DEBUG("Subdivision boundary interval start: strand/cell=" << interval.voronoi_cell_id << " t=" << t
      << " interval_idx=" << interval_idx << " mesh_id_before=" << mesh_id_before
      << " start_ref=" << format_intersection_ref(interval.start_intersection)
      << " end_ref=" << format_intersection_ref(interval.end_intersection));
    const size_t mesh_id = segment_builder_.startNewMeshFromIntersections(
      interval.voronoi_cell_id, t, interval.start_intersection, interval.end_intersection, false,
      SegmentBuilder::BoundaryEventType::Subdivision, SegmentBuilder::BoundarySegmentAction::NewSegment);
    KINDS_DEBUG("Subdivision boundary interval started: strand/cell=" << interval.voronoi_cell_id << " t=" << t
      << " interval_idx=" << interval_idx << " mesh_id=" << mesh_id);
  }
}

void SegmentBuilderSubdivisionCallback::afterEvent(KineticDelaunay::Event& e)
{
  (void)e;
}
} // namespace kinDS
