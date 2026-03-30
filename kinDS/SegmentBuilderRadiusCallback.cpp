#include "SegmentBuilderRadiusCallback.hpp"

#include "SegmentBuilder.hpp"
#include "KineticDelaunayCrossingEvent.hpp"

namespace kinDS
{
void SegmentBuilderRadiusCallback::beforeEvent(KineticDelaunay::Event& e)
{
  auto* radius = dynamic_cast<KineticDelaunay::RadiusEvent*>(&e);
  if (!radius)
  {
    return;
  }
  size_t face_id = segment_builder_.kin_del.getGraph().getHalfEdges()[radius->half_edge_id].face;
  bool is_inside = segment_builder_.kin_del.getFaceInside(face_id);
  auto& graph = segment_builder_.kin_del.getGraph();
  const auto& face_half_edges = graph.getFaces()[face_id].half_edges;

  std::array<bool, 3> is_boundary_edge;
  size_t boundary_edge_count = 0;
  for (size_t i = 0; i < 3; ++i)
  {
    size_t he_id = face_half_edges[i];
    is_boundary_edge[i] = segment_builder_.kin_del.isOnComponentBoundary(he_id);
    if (is_boundary_edge[i])
    {
      boundary_edge_count++;
    }
  }

  switch (boundary_edge_count)
  {
  case 0:
  {
    size_t vertices[3];
    for (size_t i = 0; i < 3; ++i)
    {
      vertices[i] = graph.getHalfEdges()[face_half_edges[i]].origin;

      if (vertices[i] == size_t(-1))
      {
        KINDS_ERROR(
          "Boundary triangle was turned at " << radius->occurrence_time << ", that is impossible and will be ignored!");
        break;
      }
    }
    glm::dvec2 p0 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[0]);
    glm::dvec2 p1 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[1]);
    glm::dvec2 p2 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[2]);
    glm::dvec2 new_point = (p0 + p1 + p2) / 3.0;

    size_t new_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(glm::dvec3 { new_point[0], new_point[1], radius->occurrence_time },
      glm::dvec2 { 0.0, 0.0 }, vertices[0], radius->occurrence_time);

    if (is_inside)
    {
      for (size_t i = 0; i < 3; ++i)
      {
        segment_builder_.boundary_mesh_last_left_and_right_vertex[face_half_edges[i]]
          = std::make_pair(new_vertex_index, new_vertex_index);
      }
    }
    else
    {
      for (size_t i = 0; i < 3; ++i)
      {
        size_t outer_he_id = face_half_edges[i] ^ 1;
        segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id]
          = std::make_pair(new_vertex_index, new_vertex_index);
      }
    }

    break;
  }
  case 1:
  {
    size_t boundary_he_index = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (is_boundary_edge[i])
      {
        boundary_he_index = i;
        break;
      }
    }

    size_t inner_he_id = face_half_edges[boundary_he_index];
    size_t outer_he_id = inner_he_id ^ 1;
    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);

    glm::dvec2 opposite_point = segment_builder_.kin_del.getPointAt(radius->occurrence_time, opposite_vertex);
    size_t u = graph.getHalfEdges()[inner_he_id].origin;
    glm::dvec2 p_u = segment_builder_.kin_del.getPointAt(radius->occurrence_time, u);
    size_t v = graph.getHalfEdges()[outer_he_id].origin;
    glm::dvec2 p_v = segment_builder_.kin_del.getPointAt(radius->occurrence_time, v);

    glm::dvec2 new_boundary_vertex = (opposite_point + p_u + p_v) / 3.0;

    size_t component_id = segment_builder_.kin_del.component_data.component_map[v];
    auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    auto centroid = polygonCentroid(boundary_polygon);

    size_t boundary_he_id = outer_he_id;
    if (!is_inside)
    {
      boundary_he_id = inner_he_id;
    }

    const auto& boundary_last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[boundary_he_id];

    size_t new_boundary_vertex_index = segment_builder_.addBoundaryVertex(
      glm::dvec3 { new_boundary_vertex[0], new_boundary_vertex[1], radius->occurrence_time }, centroid, opposite_vertex, radius->occurrence_time);

    size_t index
      = segment_builder_.addBoundaryTriangle(boundary_last_vertices.first, boundary_last_vertices.second, new_boundary_vertex_index);

    if (index == size_t(-1))
    {
      KINDS_DEBUG("\ninner_he_id: " << inner_he_id << "\nouter_he_id: " << outer_he_id
                                    << "\nlast_vertices: " << boundary_last_vertices.first << ", "
                                    << boundary_last_vertices.second << "\ninner last vertices: "
                                    << segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id].first << ", "
                                    << segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id].second << "\nu: " << u
                                    << ", v: " << v << ", opposite: " << opposite_vertex);
    }

    size_t he1_id = graph.getHalfEdges()[inner_he_id].next;
    size_t he2_id = graph.getHalfEdges()[he1_id].next;

    if (!is_inside)
    {
      he1_id = he1_id ^ 1;
      he2_id = he2_id ^ 1;
    }

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id]
      = std::make_pair(boundary_last_vertices.first, new_boundary_vertex_index);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id]
      = std::make_pair(new_boundary_vertex_index, boundary_last_vertices.second);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[boundary_he_id] = std::make_pair(-1, -1);

    break;
  }
  case 2:
  {
    size_t non_boundary_he_index = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      if (!is_boundary_edge[i])
      {
        non_boundary_he_index = i;
        break;
      }
    }

    size_t inner_he_id = face_half_edges[non_boundary_he_index];
    size_t outer_he_id = inner_he_id ^ 1;

    size_t opposite_vertex = graph.triangleOppositeVertex(inner_he_id);

    glm::dvec2 opposite_point = segment_builder_.kin_del.getPointAt(radius->occurrence_time, opposite_vertex);
    size_t u = graph.getHalfEdges()[inner_he_id].origin;
    glm::dvec2 p_u = segment_builder_.kin_del.getPointAt(radius->occurrence_time, u);
    size_t v = graph.getHalfEdges()[outer_he_id].origin;
    glm::dvec2 p_v = segment_builder_.kin_del.getPointAt(radius->occurrence_time, v);
    glm::dvec2 old_boundary_vertex = (opposite_point + p_u + p_v) / 3.0;

    size_t component_id = segment_builder_.kin_del.component_data.component_map[v];
    auto& boundary_polygon = segment_builder_.kin_del.component_data.component_boundaries[component_id][0];
    auto centroid = polygonCentroid(boundary_polygon);

    size_t old_boundary_vertex_index = segment_builder_.boundary_mesh.getVertices().size();
    segment_builder_.addBoundaryVertex(
      glm::dvec3 { old_boundary_vertex[0], old_boundary_vertex[1], radius->occurrence_time }, centroid, opposite_vertex, radius->occurrence_time);

    size_t he1_id = graph.getHalfEdges()[inner_he_id].next;
    size_t he2_id = graph.getHalfEdges()[he1_id].next;
    if (is_inside)
    {
      he1_id = he1_id ^ 1;
      he2_id = he2_id ^ 1;
    }

    size_t index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second, old_boundary_vertex_index);
    if (index == size_t(-1))
    {
      KINDS_DEBUG("he1_id: " << he1_id << "\ntwin: " << (he1_id ^ 1)
                             << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first << ", "
                             << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].second
                             << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].first
                             << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id ^ 1].second
                             << "\nopposite: " << opposite_vertex);
    }

    index = segment_builder_.addBoundaryTriangle(segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first,
      segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second, old_boundary_vertex_index);
    if (index == size_t(-1))
    {
      KINDS_DEBUG("he2_id: " << he2_id << "\ntwin: " << (he2_id ^ 1)
                             << "\nlast_vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].first << ", "
                             << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second
                             << "\ntwin last vertices: " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].first
                             << ", " << segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id ^ 1].second
                             << "\nopposite: " << opposite_vertex);
    }

    segment_builder_.half_edge_to_boundary_vertex_index[outer_he_id] = old_boundary_vertex_index;

    if (!is_inside)
    {
      segment_builder_.boundary_mesh_last_left_and_right_vertex[outer_he_id]
        = std::make_pair(segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first,
          segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second);
    }
    else
    {
      segment_builder_.boundary_mesh_last_left_and_right_vertex[inner_he_id]
        = std::make_pair(segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id].second,
          segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id].first);
    }

    segment_builder_.boundary_mesh_last_left_and_right_vertex[he1_id] = std::make_pair(-1, -1);
    segment_builder_.boundary_mesh_last_left_and_right_vertex[he2_id] = std::make_pair(-1, -1);
    break;
  }
  case 3:
  {
    size_t vertices[3];
    for (size_t i = 0; i < 3; ++i)
    {
      vertices[i] = graph.getHalfEdges()[face_half_edges[i]].origin;
    }
    glm::dvec2 p0 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[0]);
    glm::dvec2 p1 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[1]);
    glm::dvec2 p2 = segment_builder_.kin_del.getPointAt(radius->occurrence_time, vertices[2]);
    glm::dvec2 new_point = (p0 + p1 + p2) / 3.0;

    size_t new_vertex_index = segment_builder_.addBoundaryVertex(glm::dvec3 { new_point[0], new_point[1], radius->occurrence_time },
      glm::dvec2 { 0.0, 0.0 }, vertices[0], radius->occurrence_time);

    for (size_t i = 0; i < 3; ++i)
    {
      const auto& last_vertices = segment_builder_.boundary_mesh_last_left_and_right_vertex[face_half_edges[i]];
      segment_builder_.addBoundaryTriangle(last_vertices.first, last_vertices.second, new_vertex_index);
    }

    break;
  }
  }
}

void SegmentBuilderRadiusCallback::afterEvent(KineticDelaunay::Event& e)
{
  auto* radius = dynamic_cast<KineticDelaunay::RadiusEvent*>(&e);
  if (!radius)
  {
    return;
  }
  auto& graph = segment_builder_.kin_del.getGraph();
  auto vertices = graph.adjacentTriangleVertices(radius->half_edge_id);
  size_t component_id = segment_builder_.kin_del.component_data.component_map[vertices[0]];

  auto split = segment_builder_.kin_del.checkForSplit(vertices);
  segment_builder_.splitComponent(component_id, split, radius->occurrence_time);

  if (split.empty())
  {
    std::vector<bool> visited(segment_builder_.kin_del.getGraph().getHalfEdges().size(), false);
    segment_builder_.kin_del.component_data.component_boundaries[component_id] = segment_builder_.kin_del.extractComponentBoundaries(
      segment_builder_.kin_del.component_data.components[component_id], radius->occurrence_time, visited);
    segment_builder_.kin_del.component_data.component_centroids[component_id]
      = polygonCentroid(segment_builder_.kin_del.component_data.component_boundaries[component_id][0]);
    segment_builder_.kin_del.component_data.component_last_updated[component_id] = radius->occurrence_time;
  }

  if (segment_builder_.visual_debug)
  {
    std::vector<glm::dvec2> points = segment_builder_.kin_del.getPointsAt(radius->occurrence_time);
    std::string filename = "t" + std::to_string(radius->occurrence_time) + "_segmentbuilder_after_radius_he"
      + std::to_string(radius->half_edge_id) + ".svg";
    const auto& containing_tri_ids = segment_builder_.kin_del.getCrossingData().getContainingTriIds();
    auto intersection_debug_data = segment_builder_.kin_del.getCrossingIntersectionDebugData();
    HalfEdgeDelaunayGraphToSVG::write(points, graph, filename, 0.1, &segment_builder_.kin_del.getFacesInside(), true,
      &containing_tri_ids, &intersection_debug_data);
    KINDS_INFO("SegmentBuilder wrote SVG: " << filename);
  }

  auto triangle_he_ids = graph.getTriangleHalfEdgeIndices(radius->half_edge_id);
  for (size_t tri_he_id : triangle_he_ids)
  {
    auto d_edge_intersections = segment_builder_.kin_del.getCrossingData().delaunay_edge_intersections[tri_he_id / 2];
    for (auto intersection_ref : d_edge_intersections)
    {
      auto v_it = intersection_ref->voronoi_ref;
      (void)v_it;
    }
  }
}
} // namespace kinDS

