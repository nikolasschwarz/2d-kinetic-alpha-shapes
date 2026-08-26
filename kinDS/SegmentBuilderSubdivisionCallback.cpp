#include "SegmentBuilderSubdivisionCallback.hpp"

#include "SegmentBuilderVisualDebug.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunaySubdivisionEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilder.hpp"

#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace kinDS
{
namespace
{
std::string formatIntersectionRef(const std::optional<KineticDelaunay::CrossingData::EdgeIntersectionRef>& ref)
{
  if (!ref.has_value())
  {
    return "null";
  }
  std::ostringstream oss;
  oss << "{de=" << ref.value()->delaunay_edge_id << ", ve=" << ref.value()->voronoi_edge_id
      << ", de_param=" << ref.value()->delaunay_edge_param << "}";
  return oss.str();
}
} // namespace

void SegmentBuilderSubdivisionCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* sub = dynamic_cast<KineticDelaunay::SubdivisionEvent*>(&e);
  if (!sub)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");

  const size_t strand_id = sub->strand_id;
  const double t = sub->occurrence_time;

  if (segment_builder_.visual_debug)
  {
    auto& debug_graph = segment_builder_.kin_del.getGraph();
    const size_t runtime_branch_id = segment_builder_.kin_del.getRuntimeBranchIdForStrand(strand_id);
    writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, debug_graph, t, "before",
      "subdivision_strand" + std::to_string(strand_id) + "_seq" + std::to_string(sub->queue_sequence_),
      VisualDebugHighlight::forSubdivisionStrand(debug_graph, strand_id), runtime_branch_id,
      /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr, /*explicit_runtime_branch_ids=*/nullptr,
      sub->creation_time);
  }

  KINDS_DEBUG("Inserting subdivision for strand " << strand_id << " at t = " << t);
  auto& graph = segment_builder_.kin_del.getGraph();

  std::vector<bool> he_visited(graph.halfEdgeSlotCount(), false);
  size_t component_id = segment_builder_.kin_del.component_data.component_map[strand_id];
  segment_builder_.updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = polygonCentroid(boundary_polygon);

  // Finish regular Voronoi-edge strips around the subdivided cell.
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    segment_builder_.finishMesh(*it, t, boundary_polygon, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
  }
  segment_builder_.markInteriorMeshletsCompletedForVoronoiCell(strand_id);

  size_t new_segment_id = segment_builder_.segment_properties.size();

  // Closing mesh tracing is the authoritative boundary-interval list for finish/start.
  std::vector<SegmentBuilder::BoundaryIntersectionInterval> traced_boundary_intervals;
  size_t closing_mesh_index
    = segment_builder_.createClosingMesh(strand_id, t, boundary_polygon, centroid, &traced_boundary_intervals);

  // Finish each traced boundary interval: resolve once, extend neighbors, finish.
  // Retire (flush + mark) only after all finishes so adjacent intervals can still extend shared endpoints.
  std::vector<size_t> finished_mesh_ids;
  finished_mesh_ids.reserve(traced_boundary_intervals.size());
  for (size_t interval_idx = 0; interval_idx < traced_boundary_intervals.size(); ++interval_idx)
  {
    const auto& interval = traced_boundary_intervals[interval_idx];
    const size_t mesh_id = segment_builder_.resolveIntersectionMeshPairIndex(
      interval.voronoi_cell_id, interval.start_intersection, interval.end_intersection, t);

    std::optional<size_t> d_edge_opt;
    if (interval.start_intersection.has_value())
    {
      d_edge_opt = interval.start_intersection.value()->delaunay_edge_id;
    }
    else if (interval.end_intersection.has_value())
    {
      d_edge_opt = interval.end_intersection.value()->delaunay_edge_id;
    }
    if (d_edge_opt.has_value() && mesh_id != static_cast<size_t>(-1))
    {
      const auto refs = segment_builder_.getBoundaryIntersectionsInBoundaryOrder(d_edge_opt.value());
      if (!refs.empty())
      {
        if (interval.start_intersection.has_value())
        {
          segment_builder_.extendIntersectionMeshAtSharedCrossing(
            interval.start_intersection.value()->prev_segment_mesh_pair_index, interval.start_intersection.value(), false,
            t, SegmentBuilder::BoundaryEventType::Subdivision, SegmentBuilder::BoundarySegmentAction::SegmentRemapped,
            nullptr, true, mesh_id);
        }
        else
        {
          segment_builder_.extendIntersectionMeshAtSharedCrossing(refs.front()->prev_segment_mesh_pair_index,
            refs.front(), false, t, SegmentBuilder::BoundaryEventType::Subdivision,
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped, nullptr, true, mesh_id);
        }

        if (interval.end_intersection.has_value())
        {
          segment_builder_.extendIntersectionMeshAtSharedCrossing(
            interval.end_intersection.value()->next_segment_mesh_pair_index, interval.end_intersection.value(), true, t,
            SegmentBuilder::BoundaryEventType::Subdivision, SegmentBuilder::BoundarySegmentAction::SegmentRemapped,
            nullptr, true, mesh_id);
        }
        else
        {
          segment_builder_.extendIntersectionMeshAtSharedCrossing(refs.back()->next_segment_mesh_pair_index,
            refs.back(), true, t, SegmentBuilder::BoundaryEventType::Subdivision,
            SegmentBuilder::BoundarySegmentAction::SegmentRemapped, nullptr, true, mesh_id);
        }
      }
    }

    KINDS_DEBUG("Subdivision boundary interval finish: strand/cell=" << interval.voronoi_cell_id << " t=" << t
      << " interval_idx=" << interval_idx << " mesh_id=" << mesh_id
      << " start_ref=" << formatIntersectionRef(interval.start_intersection)
      << " end_ref=" << formatIntersectionRef(interval.end_intersection));
    segment_builder_.finishMeshFromIntersections(interval.voronoi_cell_id, t, interval.start_intersection,
      interval.end_intersection, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentCompleted, nullptr, mesh_id);

    if (mesh_id != static_cast<size_t>(-1))
    {
      finished_mesh_ids.push_back(mesh_id);
    }
  }

  for (size_t mesh_id : finished_mesh_ids)
  {
    segment_builder_.flushPendingRadiusComplementarySplitsForPair(mesh_id);
    segment_builder_.markBoundaryMeshletCompleted(mesh_id);
  }

  MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_builder_.segment_mesh_pairs[closing_mesh_index];
  segment_mesh_pair.segment_index0 = segment_builder_.strand_to_segment_indices[strand_id].back();
  segment_mesh_pair.segment_index1 = new_segment_id;

  MeshStructure::SegmentProperties properties;
  segment_builder_.segment_properties.push_back(properties);
  segment_builder_.strand_to_segment_indices[strand_id].push_back(new_segment_id);

  // Start new regular Voronoi-edge strips and extend adjacent inside strips.
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    segment_builder_.startNewMesh(*it, t, false, SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::NewSegment);

    auto& he = graph.halfEdge(*it);

    size_t adjacent_he_id = he.next;
    size_t voronoi_vertex_id = graph.halfEdge(adjacent_he_id).face;
    size_t containing_face_id = segment_builder_.kin_del.getCrossingDataContainingTriId(voronoi_vertex_id);
    bool voronoi_vertex_inside = segment_builder_.kin_del.getFacesInside()[containing_face_id];

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
      SegmentBuilder::BoundarySegmentAction::SegmentRemapped, std::nullopt, "cross",
      "subdivision_adjacent_incident", "Voronoi vertex");
    const std::string adj_face_meta = segment_builder_.composeRegularStripFaceMetadata(t, adj_even / 2, adj_even,
      static_cast<int>(adj_he.origin), static_cast<int>(adj_twin.origin), SegmentBuilder::BoundaryEventType::Subdivision,
      SegmentBuilder::BoundarySegmentAction::SegmentRemapped, "subdivision_adjacent_extend");

    glm::dvec3 vertex = segment_builder_.computeVoronoiVertex(adjacent_he_id, t);
    size_t new_vertex_index = segment_builder_.addMeshletVertex(
      adjacent_mesh, boundary_polygon, centroid, vertex, strand_id, t, false, std::optional<size_t>(voronoi_vertex_id), adj_vertex_meta);
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

  // Restart the same traced intervals on new pairs (reuse=false); no resolve — links were cleared on retire.
  for (size_t interval_idx = 0; interval_idx < traced_boundary_intervals.size(); ++interval_idx)
  {
    const auto& interval = traced_boundary_intervals[interval_idx];
    const size_t mesh_id = segment_builder_.startNewMeshFromIntersections(
      interval.voronoi_cell_id, t, interval.start_intersection, interval.end_intersection, false,
      SegmentBuilder::BoundaryEventType::Subdivision, SegmentBuilder::BoundarySegmentAction::NewSegment);
    KINDS_DEBUG("Subdivision boundary interval started: strand/cell=" << interval.voronoi_cell_id << " t=" << t
      << " interval_idx=" << interval_idx << " mesh_id=" << mesh_id
      << " start_ref=" << formatIntersectionRef(interval.start_intersection)
      << " end_ref=" << formatIntersectionRef(interval.end_intersection));
  }
}

void SegmentBuilderSubdivisionCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* sub = dynamic_cast<KineticDelaunay::SubdivisionEvent*>(&e);
  if (!sub)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");
  if (segment_builder_.diagnostics)
  {
    segment_builder_.logDiagnosticsMonitoredFaceInsideState(sub->occurrence_time, "subdivision_event");
  }
}
} // namespace kinDS
