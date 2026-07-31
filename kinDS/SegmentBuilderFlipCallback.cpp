#include "SegmentBuilderFlipCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilderVisualDebug.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>

namespace kinDS
{
namespace
{
glm::dvec3 unshiftedFlipEventPoint(const KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, size_t he_id, double t)
{
  const std::vector<size_t> strand_ids = collectFlipQuadrilateralStrandIds(graph, he_id);
  if (strand_ids.empty())
  {
    return glm::dvec3(0.0, 0.0, t);
  }

  glm::dvec2 center(0.0);
  for (size_t strand_id : strand_ids)
  {
    center += kin_del.getPointAt(strand_id, t, false, false);
  }
  center /= static_cast<double>(strand_ids.size());
  return glm::dvec3(center, t);
}

size_t runtimeBranchIdForFlipEdge(const KineticDelaunay& kin_del, const HalfEdgeDelaunayGraph& graph, size_t flip_half_edge_id)
{
  int vertex = graph.halfEdge(flip_half_edge_id).origin;
  if (vertex < 0)
  {
    vertex = graph.destination(flip_half_edge_id);
  }
  if (vertex >= 0)
  {
    return kin_del.getRuntimeBranchIdForStrand(static_cast<size_t>(vertex));
  }
  return kin_del.getRuntimeBranchIdForHalfEdge(flip_half_edge_id);
}

void logFlipMonitoredEdgeDiagnostics(SegmentBuilder& segment_builder, const HalfEdgeDelaunayGraph& graph,
  const KineticDelaunay::FlipEvent& flip, const char* phase)
{
  if (!segment_builder.diagnostics)
  {
    return;
  }
  const bool in_monitored_window = std::isfinite(flip.occurrence_time)
    && flip.occurrence_time >= std::floor(KineticDelaunay::kDiagnosticsMonitoredFlipTime)
    && flip.occurrence_time < std::floor(KineticDelaunay::kDiagnosticsMonitoredFlipTime) + 1.0;
  // Match KineticDelaunay flip diagnostics: only the monitored edge inside [floor(t), floor(t)+1).
  // Disabled monitor id (-1) never matches unset/invalid edges.
  if (!in_monitored_window
    || !KineticDelaunay::matchesDiagnosticsMonitorId(
         flip.half_edge_id / 2, KineticDelaunay::kDiagnosticsMonitoredFlipDelaunayEdgeId))
  {
    return;
  }

  std::ostringstream ctx;
  ctx << "flip_" << phase << "_he" << flip.half_edge_id << "_window_t"
      << "_d" << KineticDelaunay::kDiagnosticsMonitoredFlipDelaunayEdgeId;
  segment_builder.logDiagnosticsMonitoredDelaunayEdgeState(flip.occurrence_time, ctx.str().c_str(),
    KineticDelaunay::kDiagnosticsMonitoredFlipDelaunayEdgeId);

  std::ostringstream quad_oss;
  quad_oss << "flip monitored-edge context " << phase << " flip_he=" << flip.half_edge_id << " t="
           << flip.occurrence_time << " quad_edges=[";
  const auto quad_he_ids = graph.getQuadBoundaryHalfEdgeIndices(flip.half_edge_id / 2);
  for (size_t i = 0; i < quad_he_ids.size(); ++i)
  {
    if (i > 0)
    {
      quad_oss << ", ";
    }
    quad_oss << (quad_he_ids[i] / 2);
  }
  quad_oss << "]";
  KINDS_MONITOR(quad_oss.str());
}
} // namespace

void SegmentBuilderFlipCallback::beforeEvent(KineticDelaunay::Event& e)
{
  buffered_flip_voronoi_vertex_id_.reset();
  buffered_flip_mesh_position_.reset();
  buffered_flip_delaunay_xy_.reset();

  auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&e);
  if (!flip)
  {
    return;
  }
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "before");
  auto& graph = segment_builder_.kin_del.getGraph();

  auto vertex = graph.halfEdge(flip->half_edge_id).origin;
  if (vertex == -1)
  {
    vertex = graph.destination(flip->half_edge_id);
  }
  const size_t component_id = segment_builder_.kin_del.component_data.component_map[static_cast<size_t>(vertex)];
  const size_t runtime_branch_id = runtimeBranchIdForFlipEdge(segment_builder_.kin_del, graph, flip->half_edge_id);

  // Buffer coinciding flip-edge Voronoi coordinates once from pre-flip topology/strands. afterEvent must
  // reuse these: face↔strand association changes during the flip and recomputation introduces mismatch.
  {
    const size_t even_he = flip->half_edge_id & ~size_t { 1 };
    const int left_face = graph.halfEdge(even_he).face;
    const int right_face = graph.halfEdge(even_he ^ 1).face;
    if (left_face >= 0 && right_face >= 0)
    {
      const size_t left_vv = static_cast<size_t>(left_face);
      const size_t right_vv = static_cast<size_t>(right_face);
      if (!graph.faceHasInfiniteVertex(left_vv) && !graph.faceHasInfiniteVertex(right_vv))
      {
        const size_t canonical_vv = std::min(left_vv, right_vv);
        const glm::dvec3 fallback_profile
          = unshiftedFlipEventPoint(segment_builder_.kin_del, graph, flip->half_edge_id, flip->occurrence_time);
        const auto object_space = segment_builder_.computeMeshVoronoiVertexObjectSpace(
          canonical_vv, fallback_profile, static_cast<size_t>(vertex), flip->occurrence_time);
        buffered_flip_voronoi_vertex_id_ = canonical_vv;
        buffered_flip_mesh_position_ = object_space.position;
        buffered_flip_delaunay_xy_ = glm::dvec2(
          segment_builder_.computeVoronoiVertex(graph.face(canonical_vv).half_edges[0], flip->occurrence_time));
      }
    }
  }

  writeSegmentBuilderVisualDebugSvg(segment_builder_.visual_debug, segment_builder_.kin_del, graph,
    flip->occurrence_time, "before", "flip_he" + std::to_string(flip->half_edge_id),
    VisualDebugHighlight::forFlip(graph, flip->half_edge_id), runtime_branch_id,
    /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr, /*explicit_runtime_branch_ids=*/nullptr,
    flip->creation_time);
  logFlipMonitoredEdgeDiagnostics(segment_builder_, graph, *flip, "before");
  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = polygonCentroid(boundary_polygon);

  SegmentBuilder::MeshletVertexRuntimeInfo flip_voronoi_runtime;
  flip_voronoi_runtime.explicit_mesh_position = buffered_flip_mesh_position_;
  flip_voronoi_runtime.explicit_delaunay_xy = buffered_flip_delaunay_xy_;

  // Finish the segment mesh pair of the edge being flipped.
  // Prefer buffered Delaunay XY as the profile input; mesh-space placement comes from explicit_mesh_position.
  const glm::dvec3 event_point = buffered_flip_delaunay_xy_.has_value()
    ? glm::dvec3(buffered_flip_delaunay_xy_->x, buffered_flip_delaunay_xy_->y, flip->occurrence_time)
    : unshiftedFlipEventPoint(segment_builder_.kin_del, graph, flip->half_edge_id, flip->occurrence_time);
  size_t segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[flip->half_edge_id];

  if(graph.isInfinite(flip->half_edge_id) && segment_builder_.kin_del.computeBoundaryOnTheFly())
  {
    // Do nothing for now
  }
  else
  {
    VoronoiMesh& mesh = segment_builder_.meshes[segment_mesh_pair_index];
    const auto& last_segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    const auto flip_vertex_metadata = [](size_t voronoi_vertex_id)
    {
      return SegmentBuilder::MetadataBuilder()
        .addString("event_type", "flip_event")
        .addString("source", "Voronoi vertex")
        .addSize("voronoi_vertex_id", voronoi_vertex_id)
        .build();
    };
    const std::string flip_face_metadata = segment_builder_.store_mesh_metadata
      ? SegmentBuilder::MetadataBuilder()
          .addString("event_type", "flip_event")
          .addDouble("t", flip->occurrence_time)
          .build()
      : std::string {};
    if (!last_segments.empty())
    {
      const size_t pre_even_flip_he = flip->half_edge_id & ~1;
      const size_t pre_left_voronoi_vertex_id = static_cast<size_t>(graph.halfEdge(pre_even_flip_he).face);
      const size_t meshing_voronoi_vertex_id = buffered_flip_voronoi_vertex_id_.value_or(
        canonicalFlipEdgeVoronoiVertexIdForMeshing(graph, flip->half_edge_id, pre_left_voronoi_vertex_id));
      size_t event_vertex_index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid, event_point, vertex,
        flip->occurrence_time, false, std::optional<size_t>(meshing_voronoi_vertex_id),
        flip_vertex_metadata(meshing_voronoi_vertex_id), std::nullopt, flip_voronoi_runtime);
      size_t last_left = last_segments.front().mesh_start_vertex_id;
      size_t last_right = last_segments.back().mesh_end_vertex_id;
      // create one triangle to the event point
      {
        const size_t tris_before = mesh.getTriangleCount();
        segment_builder_.addMeshletTriangle(mesh, last_left, last_right, event_vertex_index, flip_face_metadata);
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

    glm::dvec2 new_boundary_vertex = segment_builder_.kin_del.getPointAt(flip->occurrence_time, opposite_vertex, false, false);

    size_t new_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    // TODO: raw UVs
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { new_boundary_vertex[0], new_boundary_vertex[1], flip->occurrence_time }, centroid, opposite_vertex, flip->occurrence_time, false);

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
  SegmentBuilder::ScopedMetadataCallbackPhase callback_phase(segment_builder_, "after");
  auto& graph = segment_builder_.kin_del.getGraph();

  if (segment_builder_.diagnostics)
  {
    segment_builder_.kin_del.validateFlipAdjacentFaceInsideConsistency(flip->half_edge_id, flip->occurrence_time);
    segment_builder_.logDiagnosticsMonitoredFaceInsideState(flip->occurrence_time, "flip_event");
  }
  logFlipMonitoredEdgeDiagnostics(segment_builder_, graph, *flip, "after_pre_refresh");

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
  const size_t left_voronoi_vertex_id = static_cast<size_t>(graph.halfEdge(even_flip_he).face);
  const size_t left_containing_tri_id = segment_builder_.kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);
  const bool left_inside = segment_builder_.kin_del.getFaceInside(left_containing_tri_id);
  const bool flip_pos_finite = std::isfinite(flip->position[0]) && std::isfinite(flip->position[1]);
  const bool seed_mesh_with_flip_vertex = left_inside && flip_pos_finite;
  // Prefer pre-flip buffered coincidence coords; post-flip strand/face association would reintroduce mismatch.
  const glm::dvec3 event_point = buffered_flip_delaunay_xy_.has_value()
    ? glm::dvec3(buffered_flip_delaunay_xy_->x, buffered_flip_delaunay_xy_->y, flip->occurrence_time)
    : unshiftedFlipEventPoint(segment_builder_.kin_del, graph, flip->half_edge_id, flip->occurrence_time);
  const auto flip_vertex_metadata = [](size_t voronoi_vertex_id)
  {
    return SegmentBuilder::MetadataBuilder()
      .addString("event_type", "flip_event")
      .addString("source", "Voronoi vertex")
      .addSize("voronoi_vertex_id", voronoi_vertex_id)
      .build();
  };
  const std::string flip_face_metadata = segment_builder_.store_mesh_metadata
    ? SegmentBuilder::MetadataBuilder()
        .addString("event_type", "flip_event")
        .addDouble("t", flip->occurrence_time)
        .build()
    : std::string {};
  SegmentBuilder::MeshletVertexRuntimeInfo flip_voronoi_runtime;
  flip_voronoi_runtime.explicit_mesh_position = buffered_flip_mesh_position_;
  flip_voronoi_runtime.explicit_delaunay_xy = buffered_flip_delaunay_xy_;

  // For now also create a mesh, but this might be changed later
  VoronoiMesh mesh;

  // add last vertex indices
  segment_builder_.segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  if (seed_mesh_with_flip_vertex)
  {
    const size_t meshing_voronoi_vertex_id = buffered_flip_voronoi_vertex_id_.value_or(
      canonicalFlipEdgeVoronoiVertexIdForMeshing(graph, flip->half_edge_id, left_voronoi_vertex_id));
    size_t index = segment_builder_.addMeshletVertex(mesh, boundary_polygon, centroid, event_point, vertex,
      flip->occurrence_time, false, std::optional<size_t>(meshing_voronoi_vertex_id),
      flip_vertex_metadata(meshing_voronoi_vertex_id), std::nullopt, flip_voronoi_runtime);
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

    // Open endpoint of this directed half-edge: its face is a flip-edge Voronoi vertex. Reuse the
    // pre-flip buffered coincidence coords so before/after meshlets share identical positions.
    const int neighbor_voronoi_vertex_i = graph.halfEdge(he_id).face;
    if (neighbor_voronoi_vertex_i < 0)
    {
      KINDS_WARNING("Flip neighbor extension: half-edge " << he_id
                                                           << " has no incident Voronoi vertex at t="
                                                           << flip->occurrence_time << "; skipping.");
      continue;
    }
    const size_t neighbor_voronoi_vertex_id = static_cast<size_t>(neighbor_voronoi_vertex_i);
    const size_t meshing_voronoi_vertex_id = buffered_flip_voronoi_vertex_id_.value_or(
      canonicalFlipEdgeVoronoiVertexIdForMeshing(graph, flip->half_edge_id, neighbor_voronoi_vertex_id));
    size_t new_vertex_index
      = segment_builder_.addMeshletVertex(mesh_ref, boundary_polygon, centroid, event_point, vertex,
        flip->occurrence_time, false, std::make_optional(meshing_voronoi_vertex_id),
        flip_vertex_metadata(meshing_voronoi_vertex_id), std::nullopt, flip_voronoi_runtime);

    int he_id_left = segments.front().start_half_edge_id;
    int he_id_right = segments.back().end_half_edge_id;

    if (he_id % 2 == 0 && he_id_left == -1)
    {
      size_t last_left = segments.front().mesh_start_vertex_id;
      size_t last_right = segments.front().mesh_end_vertex_id;
      {
        const size_t tris_before = mesh_ref.getTriangleCount();
        segment_builder_.addMeshletTriangle(mesh_ref, last_left, last_right, new_vertex_index, flip_face_metadata);
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
        segment_builder_.addMeshletTriangle(mesh_ref, last_left, last_right, new_vertex_index, flip_face_metadata);
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
    runtimeBranchIdForFlipEdge(segment_builder_.kin_del, graph, flip->half_edge_id),
    /*separation_offset_segments=*/nullptr, /*seam_outlines=*/nullptr, /*explicit_runtime_branch_ids=*/nullptr,
    flip->creation_time);

  if (segment_builder_.kin_del.isOnComponentBoundary(flip->half_edge_id))
  {
    size_t outer_he_id = segment_builder_.kin_del.isOnComponentBoundaryOutside(flip->half_edge_id) ? flip->half_edge_id
                                                                                                     : graph.twin(flip->half_edge_id);
    size_t inner_he_id = outer_he_id ^ 1;

    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);
    const auto& boundary_last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id];

    glm::dvec2 old_boundary_vertex = segment_builder_.kin_del.getPointAt(flip->occurrence_time, opposite_vertex, false, false);

    size_t old_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { old_boundary_vertex[0], old_boundary_vertex[1], flip->occurrence_time }, centroid, opposite_vertex, flip->occurrence_time, false);

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
  logFlipMonitoredEdgeDiagnostics(segment_builder_, graph, *flip, "after_post_refresh");

  // Force-recompute crossing params on every Delaunay edge of the flipped quad and restore list order.
  // Params stamped at this same timestamp before the flip can still be geometrically wrong; mesh-pair
  // edge links are all unset (-1) here, so sorting cannot break them.
  {
    std::unordered_set<size_t> affected_delaunay_edges;
    affected_delaunay_edges.insert(flip->half_edge_id / 2);
    for (size_t he_id : graph.getQuadBoundaryHalfEdgeIndices(flip->half_edge_id / 2))
    {
      affected_delaunay_edges.insert(he_id / 2);
    }
    for (size_t delaunay_edge_id : affected_delaunay_edges)
    {
      if (graph.isInfinite(2 * delaunay_edge_id) && segment_builder_.kin_del.computeBoundaryOnTheFly())
      {
        continue;
      }
      segment_builder_.kin_del.refreshAndSortDelaunayEdgeIntersectionParams(
        delaunay_edge_id, flip->occurrence_time);
    }
  }

  buffered_flip_voronoi_vertex_id_.reset();
  buffered_flip_mesh_position_.reset();
  buffered_flip_delaunay_xy_.reset();
}
} // namespace kinDS

