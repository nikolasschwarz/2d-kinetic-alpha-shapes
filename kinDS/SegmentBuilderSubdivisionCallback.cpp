#include "SegmentBuilderSubdivisionCallback.hpp"

#include "HalfEdgeDelaunayGraphToSVG.hpp"
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
void SegmentBuilderSubdivisionCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* sub = dynamic_cast<KineticDelaunay::SubdivisionEvent*>(&e);
  if (!sub)
  {
    return;
  }

  const size_t strand_id = sub->strand_id;
  const double t = sub->occurrence_time;

  KINDS_DEBUG("Inserting subdivision for strand " << strand_id << " at t = " << t);
  //   Traverse all half-edges around this strand and insert a new vertex into the corresponding segment meshes
  auto& graph = segment_builder_.kin_del.getGraph();

  // compute boundary polygon for component at time t
  std::vector<bool> he_visited(graph.getHalfEdges().size(), false);
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

  // finish old meshes
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    segment_builder_.finishMesh(*it, t, boundary_polygon);
  }

  // Finish all boundary-interval meshes on boundary Delaunay edges of this Voronoi cell.
  for (size_t he_id : cell_boundary_he_ids)
  {
    const size_t he_even = he_id & ~1;
    if (!segment_builder_.kin_del.isOnComponentBoundary(he_even))
    {
      continue;
    }
    const size_t d_edge_id = he_even / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs;
    refs.reserve(d_intersections.size());
    for (const auto& ref : d_intersections)
    {
      refs.push_back(ref);
    }

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      segment_builder_.finishMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      segment_builder_.finishMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      segment_builder_.finishMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::SegmentCompleted);
    }
  }

  size_t new_segment_id = segment_builder_.segment_properties.size();

  // create a closing mesh
  size_t closing_mesh_index = segment_builder_.createClosingMesh(strand_id, t, boundary_polygon, centroid);
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
    segment_builder_.startNewMesh(*it, t);

    // insert vertices into adjacent meshes
    auto& he = graph.getHalfEdges()[*it];

    size_t adjacent_he_id = he.next;
    size_t voronoi_vertex_id = graph.getHalfEdges()[adjacent_he_id].face;
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

    glm::dvec3 vertex = segment_builder_.computeVoronoiVertex(adjacent_he_id, t);
    size_t new_vertex_index = segment_builder_.addMeshletVertex(
      adjacent_mesh, boundary_polygon, centroid, vertex, strand_id, t, std::optional<size_t>(voronoi_vertex_id));
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
          segment_builder_.addMeshletTriangle(adjacent_mesh, last_left, last_right, new_vertex_index);
          std::ostringstream note;
          note << "extend_subdivision_adjacent d_tris=" << (adjacent_mesh.getTriangleCount() - tris_before);
          segment_builder_.meshletDiagnosticLogLine("extend_mesh", adjacent_he_id, t, note.str().c_str());
        }
        segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
      }
      else if (adjacent_he_id % 2 != 0 && he_id_right == -1)
      {
        size_t last_left = segments.back().mesh_start_vertex_id;
        size_t last_right = segments.back().mesh_end_vertex_id;
        {
          const size_t tris_before = adjacent_mesh.getTriangleCount();
          segment_builder_.addMeshletTriangle(adjacent_mesh, last_left, last_right, new_vertex_index);
          std::ostringstream note;
          note << "extend_subdivision_adjacent d_tris=" << (adjacent_mesh.getTriangleCount() - tris_before);
          segment_builder_.meshletDiagnosticLogLine("extend_mesh", adjacent_he_id, t, note.str().c_str());
        }
        segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
      }
    }
  }

  // Restart all boundary-interval meshes for the updated strand topology after subdivision.
  for (size_t he_id : cell_boundary_he_ids)
  {
    const size_t he_even = he_id & ~1;
    if (!segment_builder_.kin_del.isOnComponentBoundary(he_even))
    {
      continue;
    }
    const size_t d_edge_id = he_even / 2;
    const auto& d_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[d_edge_id];
    if (d_intersections.empty())
    {
      continue;
    }

    std::vector<KineticDelaunay::CrossingData::EdgeIntersectionRef> refs;
    refs.reserve(d_intersections.size());
    for (const auto& ref : d_intersections)
    {
      refs.push_back(ref);
    }

    {
      const size_t first_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, std::nullopt, refs.front());
      const bool reuse_existing = (first_cell != strand_id);
      segment_builder_.startNewMeshFromIntersections(
        first_cell, t, std::nullopt, refs.front(), reuse_existing, SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::NewSegment);
    }
    for (size_t k = 0; k + 1 < refs.size(); ++k)
    {
      const size_t mid_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs[k], refs[k + 1]);
      const bool reuse_existing = (mid_cell != strand_id);
      segment_builder_.startNewMeshFromIntersections(
        mid_cell, t, refs[k], refs[k + 1], reuse_existing, SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::NewSegment);
    }
    {
      const size_t last_cell
        = segment_builder_.determineVoronoiCellForBoundaryIntersectionInterval(d_edge_id, refs.back(), std::nullopt);
      const bool reuse_existing = (last_cell != strand_id);
      segment_builder_.startNewMeshFromIntersections(
        last_cell, t, refs.back(), std::nullopt, reuse_existing, SegmentBuilder::BoundaryEventType::Subdivision,
        SegmentBuilder::BoundarySegmentAction::NewSegment);
    }
  }
}

void SegmentBuilderSubdivisionCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* sub = dynamic_cast<KineticDelaunay::SubdivisionEvent*>(&e);
  if (!sub)
  {
    return;
  }
  if (!segment_builder_.visual_debug)
  {
    return;
  }

  auto& graph = segment_builder_.kin_del.getGraph();
  std::vector<glm::dvec2> points = segment_builder_.kin_del.getPointsAt(sub->occurrence_time);
  std::string filename = "t" + std::to_string(sub->occurrence_time) + "_segmentbuilder_after_subdivision_strand"
    + std::to_string(sub->strand_id) + "_seq" + std::to_string(sub->queue_sequence_) + ".svg";
  const auto& containing_tri_ids = segment_builder_.kin_del.getCrossingData().getContainingTriIds();
  auto intersection_debug_data = segment_builder_.kin_del.getCrossingIntersectionDebugData();
  HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &segment_builder_.kin_del.getFacesInside(), true,
    &containing_tri_ids, &intersection_debug_data);
  KINDS_INFO("SegmentBuilder wrote SVG: " << filename);
}
} // namespace kinDS
