#include "SegmentBuilderFlipCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "SegmentBuilderVisualDebug.hpp"

#include <cmath>
#include <optional>
#include <sstream>
#include <string>

namespace kinDS
{
void SegmentBuilderFlipCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&e);
  if (!flip)
  {
    return;
  }
  auto& graph = segment_builder_.kin_del.getGraph();

  auto vertex = graph.halfEdge(flip->half_edge_id).origin;
  if (vertex == -1)
  {
    vertex = graph.destination(flip->half_edge_id);
  }
  const size_t component_id = segment_builder_.kin_del.component_data.component_map[static_cast<size_t>(vertex)];
  const size_t runtime_branch_id = segment_builder_.kin_del.getRuntimeBranchIdForHalfEdge(flip->half_edge_id);

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    flip->occurrence_time, "before", "flip_he" + std::to_string(flip->half_edge_id),
    VisualDebugHighlight::forFlip(graph, flip->half_edge_id), runtime_branch_id);
  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = polygonCentroid(boundary_polygon);

  // Finish the segment mesh pair of the edge being flipped
  glm::dvec3 event_point { flip->position[0], flip->position[1], flip->occurrence_time };
  size_t segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[flip->half_edge_id];

  if(graph.isInfinite(flip->half_edge_id) && segment_builder_.kin_del.computeBoundaryOnTheFly())
  {
    // Do nothing for now
  }
  else
  {
    VoronoiMesh& mesh = segment_builder_.meshes[segment_mesh_pair_index];
    const auto& last_segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    if (!last_segments.empty())
    {
      const size_t pre_even_flip_he = flip->half_edge_id & ~1;
      const size_t pre_left_voronoi_vertex_id = graph.halfEdge(pre_even_flip_he).face;
      size_t event_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid, event_point, vertex,
        flip->occurrence_time, std::optional<size_t>(pre_left_voronoi_vertex_id));
      size_t last_left = last_segments.front().mesh_start_vertex_id;
      size_t last_right = last_segments.back().mesh_end_vertex_id;
      // create one triangle to the event point
      {
        const size_t tris_before = mesh.getTriangleCount();
        segment_builder_.addMeshletTriangle(mesh, last_left, last_right, event_vertex_index);
        if (segment_builder_.diagnostics)
        {
          std::ostringstream note;
          note << "extend_flip_before d_tris=" << (mesh.getTriangleCount() - tris_before);
          segment_builder_.meshletDiagnosticLogLine(
            "extend_mesh", flip->half_edge_id, flip->occurrence_time, note.str().c_str());
        }
      }
    }
    else // TODO: I think this should't occur, but it does
    {
      const size_t pre_even_flip_he = flip->half_edge_id & ~1;
      const size_t pre_left_voronoi_vertex_id = graph.halfEdge(pre_even_flip_he).face;
      const size_t pre_left_containing_tri_id
        = segment_builder_.kin_del.getCrossingDataContainingTriId(pre_left_voronoi_vertex_id);
      const bool pre_left_inside = segment_builder_.kin_del.getFaceInside(pre_left_containing_tri_id);
      const bool event_pos_finite = std::isfinite(flip->position[0]) && std::isfinite(flip->position[1]);
      // If this segment list is empty, afterEvent skipped seeding this mesh because the Voronoi vertex was not usable.
      assert(!(pre_left_inside && event_pos_finite));
    }
  }

  // For the boundary mesh, handle the case that a boundary edge is flipped. This means the opposite vertex becomes a
  // boundary vertex
  // Only applies if the edge is on the component boundary as well
  if (segment_builder_.kin_del.isOnComponentBoundary(flip->half_edge_id))
  {
    // TODO: make sure this is equivalent to component boundary
    size_t outer_he_id = segment_builder_.kin_del.isOnComponentBoundaryOutside(flip->half_edge_id) ? flip->half_edge_id
                                                                                                     : graph.twin(flip->half_edge_id);
    size_t inner_he_id = outer_he_id ^ 1;

    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);
    const auto& boundary_last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id];

    glm::dvec2 new_boundary_vertex = segment_builder_.kin_del.getPointAt(flip->occurrence_time, opposite_vertex);

    size_t new_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    // TODO: raw UVs
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { new_boundary_vertex[0], new_boundary_vertex[1], flip->occurrence_time }, centroid, opposite_vertex, flip->occurrence_time);

    // create one triangle to the event point
    segment_builder_.addBoundaryTriangle(boundary_last_vertices.first, boundary_last_vertices.second, new_boundary_vertex_index);

    // update last left and right indices of the other two half-edges of the triangle
    size_t he1_id = graph.halfEdge(inner_he_id).next;
    size_t he2_id = graph.halfEdge(he1_id).next;

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id]
      = std::make_pair(boundary_last_vertices.first, new_boundary_vertex_index);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id]
      = std::make_pair(new_boundary_vertex_index, boundary_last_vertices.second);

    // reset last left and right vertices of the half-edge because it is not on the boundary anymore
    segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id] = std::make_pair(-1, -1);
  }
}

void SegmentBuilderFlipCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&e);
  if (!flip)
  {
    return;
  }
  auto& graph = segment_builder_.kin_del.getGraph();

  int vertex = graph.halfEdge(flip->half_edge_id).origin;
  if (vertex == -1)
  {
    vertex = graph.destination(flip->half_edge_id);
  }
  const size_t component_id = segment_builder_.kin_del.component_data.component_map[static_cast<size_t>(vertex)];
  segment_builder_.updateBoundaries(flip->occurrence_time, { component_id });

  const auto& he = graph.halfEdge(flip->half_edge_id);
  const auto& twin_he = graph.halfEdge(flip->half_edge_id ^ 1);
  // Create a new segment mesh pair for the two new edges created by the flip
  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : segment_builder_.strand_to_segment_indices[he.origin].back();
  segment_mesh_pair.segment_index1
    = twin_he.origin == -1 ? -1 : segment_builder_.strand_to_segment_indices[twin_he.origin].back();

  segment_builder_.half_edge_index_to_segment_mesh_pair_index[flip->half_edge_id] = segment_builder_.segment_mesh_pairs.size();
  segment_builder_.half_edge_index_to_segment_mesh_pair_index[flip->half_edge_id ^ 1]
    = segment_builder_.segment_mesh_pairs.size();

  segment_builder_.segment_mesh_pairs.push_back(segment_mesh_pair);

  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = polygonCentroid(boundary_polygon);

  const size_t even_flip_he = flip->half_edge_id & ~1;
  const size_t left_voronoi_vertex_id = graph.halfEdge(even_flip_he).face;
  const size_t left_containing_tri_id = segment_builder_.kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);
  const bool left_inside = segment_builder_.kin_del.getFaceInside(left_containing_tri_id);
  const bool flip_pos_finite = std::isfinite(flip->position[0]) && std::isfinite(flip->position[1]);
  const bool seed_mesh_with_flip_vertex = left_inside && flip_pos_finite;

  // For now also create a mesh, but this might be changed later
  VoronoiMesh mesh;

  // add last vertex indices
  segment_builder_.segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  if (seed_mesh_with_flip_vertex)
  {
    size_t index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid,
      glm::dvec3 { flip->position[0], flip->position[1], flip->occurrence_time }, vertex, flip->occurrence_time,
      std::optional<size_t>(left_voronoi_vertex_id));
    segment_builder_.segment_mesh_pair_last_left_and_right_vertex.back().emplace_back(
      SegmentBuilder::MeshingData { static_cast<int>(index), static_cast<int>(index), -1, -1 });
  }

  segment_builder_.registerMeshletWithSuffix(std::move(mesh), std::string("_voronoi") + std::to_string(flip->half_edge_id / 2),
    flip->occurrence_time);

  auto quad_he_ids = graph.getQuadBoundaryHalfEdgeIndices(flip->half_edge_id / 2);
  
  // Update the meshes of the neighboring segments to connect to the new vertex created at the flip position, if applicable
  for (size_t he_id : quad_he_ids)
  {
    if(segment_builder_.kin_del.getGraph().isInfinite(he_id) && segment_builder_.kin_del.computeBoundaryOnTheFly())
    {
      continue;
    }

    size_t segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[he_id];
    VoronoiMesh& mesh_ref = segment_builder_.meshes[segment_mesh_pair_index];

    auto& segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    if (segments.empty())
    {
      continue;
    }

    // Same geometric seed as above; alpha warning already emitted from addMeshletVertex for the new meshlet.
    size_t new_vertex_index = segment_builder_.addMeshletVertex(mesh_ref, boundary_polygon, centroid,
      glm::dvec3 { flip->position[0], flip->position[1], flip->occurrence_time }, vertex, flip->occurrence_time);

    int he_id_left = segments.front().start_half_edge_id;
    int he_id_right = segments.back().end_half_edge_id;

    if (he_id % 2 == 0 && he_id_left == -1)
    {
      size_t last_left = segments.front().mesh_start_vertex_id;
      size_t last_right = segments.front().mesh_end_vertex_id;
      {
        const size_t tris_before = mesh_ref.getTriangleCount();
        segment_builder_.addMeshletTriangle(mesh_ref, last_left, last_right, new_vertex_index);
        if (segment_builder_.diagnostics)
        {
          std::ostringstream note;
          note << "extend_flip_neighbor d_tris=" << (mesh_ref.getTriangleCount() - tris_before);
          segment_builder_.meshletDiagnosticLogLine(
            "extend_mesh", he_id, flip->occurrence_time, note.str().c_str());
        }
      }
      segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
    }
    else if (he_id % 2 != 0 && he_id_right == -1)
    {
      size_t last_left = segments.back().mesh_start_vertex_id;
      size_t last_right = segments.back().mesh_end_vertex_id;
      {
        const size_t tris_before = mesh_ref.getTriangleCount();
        segment_builder_.addMeshletTriangle(mesh_ref, last_left, last_right, new_vertex_index);
        if (segment_builder_.diagnostics)
        {
          std::ostringstream note;
          note << "extend_flip_neighbor d_tris=" << (mesh_ref.getTriangleCount() - tris_before);
          segment_builder_.meshletDiagnosticLogLine(
            "extend_mesh", he_id, flip->occurrence_time, note.str().c_str());
        }
      }
      segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
    }
  }

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    flip->occurrence_time, "after", "flip_he" + std::to_string(flip->half_edge_id),
    VisualDebugHighlight::forFlip(graph, flip->half_edge_id),
    segment_builder_.kin_del.getRuntimeBranchIdForHalfEdge(flip->half_edge_id));

  if (segment_builder_.kin_del.isOnComponentBoundary(flip->half_edge_id))
  {
    size_t outer_he_id = segment_builder_.kin_del.isOnComponentBoundaryOutside(flip->half_edge_id) ? flip->half_edge_id
                                                                                                     : graph.twin(flip->half_edge_id);
    size_t inner_he_id = outer_he_id ^ 1;

    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);
    const auto& boundary_last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id];

    glm::dvec2 old_boundary_vertex = segment_builder_.kin_del.getPointAt(flip->occurrence_time, opposite_vertex);

    size_t old_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { old_boundary_vertex[0], old_boundary_vertex[1], flip->occurrence_time }, centroid, opposite_vertex, flip->occurrence_time);

    size_t he1_id = graph.halfEdge(inner_he_id).next;
    size_t he2_id = graph.halfEdge(he1_id).next;

    size_t tri_index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second, old_boundary_vertex_index);

    if (tri_index == size_t(-1))
    {
      KINDS_DEBUG("\he1_id: " << he1_id << "\ntwin: " << (he1_id ^ 1)
                              << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first << ", "
                              << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second
                              << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].first
                              << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].second
                              << "\nopposite: " << opposite_vertex);
    }

    tri_index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second, old_boundary_vertex_index);

    if (tri_index == size_t(-1))
    {
      KINDS_DEBUG("\he2_id: " << he2_id << "\ntwin: " << (he2_id ^ 1)
                              << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first << ", "
                              << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second
                              << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].first
                              << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].second
                              << "\nopposite: " << opposite_vertex);
    }

    segment_builder_.half_edge_to_boundary_vertex_index[outer_he_id] = old_boundary_vertex_index;
    segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id] = std::make_pair(
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second);

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id] = std::make_pair(-1, -1);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id] = std::make_pair(-1, -1);
  }

  segment_builder_.refreshCrossingRefsForAllStrips();
}
} // namespace kinDS

