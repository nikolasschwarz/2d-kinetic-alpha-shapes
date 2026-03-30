#include "SegmentBuilderCrossingCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

namespace kinDS
{
void SegmentBuilderCrossingCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&e);
  if (!crossing)
  {
    return;
  }
  while (segment_builder_.subdivision_index < segment_builder_.subdivisions.size()
    && segment_builder_.subdivisions[segment_builder_.subdivision_index].second <= crossing->occurrence_time)
  {
    segment_builder_.insertSubdivision(segment_builder_.subdivisions[segment_builder_.subdivision_index].first,
      segment_builder_.subdivisions[segment_builder_.subdivision_index].second);
    segment_builder_.subdivision_index++;
  }
}

void SegmentBuilderCrossingCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&e);
  if (!crossing)
  {
    return;
  }
  auto& graph = segment_builder_.kin_del.getGraph();

  if (segment_builder_.visual_debug)
  {
    std::vector<glm::dvec2> points = segment_builder_.kin_del.getPointsAt(crossing->occurrence_time);
    size_t old_tri = graph.getHalfEdges()[crossing->half_edge_id].face;
    size_t new_tri = graph.getHalfEdges()[crossing->half_edge_id ^ 1].face;
    std::string filename = "t" + std::to_string(crossing->occurrence_time) + "_segmentbuilder_after_crossing_v"
      + std::to_string(crossing->voronoi_vertex_id) + "_" + std::to_string(old_tri) + "_to_" + std::to_string(new_tri) + ".svg";
    const auto& containing_tri_ids = segment_builder_.kin_del.getCrossingData().getContainingTriIds();
    auto intersection_debug_data = segment_builder_.kin_del.getCrossingIntersectionDebugData();
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &segment_builder_.kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data);
    KINDS_INFO("SegmentBuilder wrote SVG: " << filename);
  }

  size_t voronoi_vertex_id = crossing->voronoi_vertex_id;
  glm::dvec3 voronoi_vertex_position = glm::dvec3(crossing->position, crossing->occurrence_time);
  auto half_edges = graph.getFaces()[voronoi_vertex_id].half_edges;
  if (!segment_builder_.kin_del.isOnComponentBoundary(crossing->half_edge_id))
  {
    return;
  }

  size_t inside_he_id = crossing->half_edge_id;
  bool entering_boundary = false;
  if (segment_builder_.kin_del.isOnComponentBoundaryOutside(crossing->half_edge_id))
  {
    inside_he_id = crossing->half_edge_id ^ 1;
    entering_boundary = true;
  }

  size_t component_id = segment_builder_.kin_del.component_data.component_map[graph.destination(inside_he_id)];
  auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
  auto centroid = segment_builder_.kin_del.component_data.component_centroids[component_id];

  for (size_t voronoi_he_id : half_edges)
  {
    auto& segment_mesh_pair_index = segment_builder_.half_edge_index_to_segment_mesh_pair_index[voronoi_he_id];
    VoronoiMesh& mesh = segment_builder_.meshes[segment_mesh_pair_index];

    auto& segments = segment_builder_.segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    auto it = std::find_if(segments.begin(), segments.end(),
      [inside_he_id](const SegmentBuilder::MeshingData& data)
      {
        return data.end_half_edge_id == static_cast<int>(inside_he_id)
          || data.start_half_edge_id == static_cast<int>(inside_he_id);
      });
    if (it != segments.end())
    {
      if (!entering_boundary)
      {
        size_t new_vertex_index = segment_builder_.addMeshletVertex(
          mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
        int last_left = it->mesh_start_vertex_id;
        int last_right = it->mesh_end_vertex_id;
        segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
        segments.erase(it);
      }
      else
      {
        if (it->start_half_edge_id == static_cast<int>(inside_he_id))
        {
          it->start_half_edge_id = -1;
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
          it->mesh_start_vertex_id = static_cast<int>(new_vertex_index);
        }
        else
        {
          assert(it->end_half_edge_id == static_cast<int>(inside_he_id));
          it->end_half_edge_id = -1;
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          int last_left = it->mesh_start_vertex_id;
          int last_right = it->mesh_end_vertex_id;
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
          it->mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }
    else
    {
      size_t even_id = voronoi_he_id & ~1;
      size_t odd_id = even_id + 1;
      size_t end_voronoi_vertex_id = graph.getHalfEdges()[odd_id].face;

      if (end_voronoi_vertex_id == voronoi_vertex_id)
      {
        if (entering_boundary)
        {
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          segments.emplace_back(SegmentBuilder::MeshingData { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index),
            static_cast<int>(inside_he_id), -1 });
        }
        else
        {
          assert(segments.back().end_half_edge_id == -1);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          int last_left = segments.back().mesh_start_vertex_id;
          int last_right = segments.back().mesh_end_vertex_id;
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
          segments.back().end_half_edge_id = static_cast<int>(inside_he_id);
          segments.back().mesh_end_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
      else
      {
        if (entering_boundary)
        {
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          segments.emplace_front(SegmentBuilder::MeshingData { static_cast<int>(new_vertex_index), static_cast<int>(new_vertex_index),
            -1, static_cast<int>(inside_he_id) });
        }
        else
        {
          assert(segments.front().start_half_edge_id == -1);
          size_t new_vertex_index = segment_builder_.addMeshletVertex(
            mesh, boundary_polygon, centroid, voronoi_vertex_position, voronoi_vertex_id, crossing->occurrence_time);
          int last_left = segments.front().mesh_start_vertex_id;
          int last_right = segments.front().mesh_end_vertex_id;
          segment_builder_.addMeshletTriangle(mesh, last_left, last_right, new_vertex_index);
          segments.front().start_half_edge_id = static_cast<int>(inside_he_id);
          segments.front().mesh_start_vertex_id = static_cast<int>(new_vertex_index);
        }
      }
    }
  }
}
} // namespace kinDS

