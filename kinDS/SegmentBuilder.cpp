#include "SegmentBuilder.hpp"

#include "KineticDelaunayCrossingEvent.hpp"
#include "KineticDelaunayFlipEvent.hpp"
#include "KineticDelaunayRadiusEvent.hpp"
#include "Logger.hpp"
#include "SegmentBuilderCrossingCallback.hpp"
#include "SegmentBuilderFlipCallback.hpp"
#include "SegmentBuilderRadiusCallback.hpp"
#include "SegmentBuilderSectionCallback.hpp"
#include "SegmentBuilderSubdivisionCallback.hpp"

#include <cassert>
#include <cmath>
#include <glm/gtx/exterior_product.hpp>
#include <map>
#include <string>
#include <unordered_map>

using namespace kinDS;

static bool raySegmentIntersection(
  const glm::dvec2& C, const glm::dvec2& D, const glm::dvec2& A, const glm::dvec2& B, double& t_out)
{
  glm::dvec2 E = B - A;
  glm::dvec2 AC = A - C;

  double det = glm::cross(D, E);
  if (std::abs(det) < 1e-12)
    return false;

  double t = glm::cross(AC, E) / det;
  double u = glm::cross(AC, D) / det;

  if (t >= 0.0 && u >= 0.0 && u <= 1.0)
  {
    t_out = t;
    return true;
  }
  return false;
}

static std::vector<double> rayCast(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i].p;
    const glm::dvec2& B = polygon[(i + 1) % n].p;

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static std::vector<double> rayCast(
  const std::vector<glm::dvec2>& polygon, const glm::dvec2& origin, const glm::dvec2& dir)
{
  double lenCP = glm::length(dir);

  if (lenCP < 1e-12)
    return {};

  std::vector<double> hits;
  const size_t n = polygon.size();
  for (size_t i = 0; i < n; ++i)
  {
    const glm::dvec2& A = polygon[i];
    const glm::dvec2& B = polygon[(i + 1) % n];

    double t;
    if (raySegmentIntersection(origin, dir, A, B, t))
    {
      hits.emplace_back(t);
    }
  }

  return hits;
}

static double relativeDistanceFromCenter(
  const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  glm::dvec2 dir = point - center;
  auto hits = rayCast(polygon, center, dir);

  if (hits.empty())
    return std::numeric_limits<double>::quiet_NaN();

  double t_max = -std::numeric_limits<double>::infinity();

  for (double t : hits)
  {
    t_max = std::max(t_max, t);
  }

  if (t_max < 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // |B - C| = t_max * |D|
  return 1.0 / t_max;
}

static bool isInside(const std::vector<BoundaryPoint>& polygon, const glm::dvec2& center, const glm::dvec2& point)
{
  double rel_dist = relativeDistanceFromCenter(polygon, center, point);
  if (std::isnan(rel_dist))
  {
    throw std::runtime_error("Center lies outside of polygon");
  }
  return rel_dist <= 1.0;
}

[[nodiscard]] glm::dvec3 kinDS::SegmentBuilder::computeVoronoiVertex(size_t half_edge_id, double t) const
{
  return kin_del.computeVoronoiVertexClampedInfinity(half_edge_id, t);
}

bool clampVoronoiVertices(glm::dvec3& left_vertex, glm::dvec3& right_vertex,
  const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  // don't do this for now
  return true;

  bool left_inside = isInside(boundary_points, centroid, glm::dvec2 { left_vertex[0], left_vertex[1] });
  bool right_inside = isInside(boundary_points, centroid, glm::dvec2 { right_vertex[0], right_vertex[1] });

  if (left_inside && !right_inside)
  {
    KINDS_DEBUG(
      "Clamping right vertex: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
    // clamp right to the boundary
    glm::dvec2 origin { left_vertex[0], left_vertex[1] };
    glm::dvec2 dir { right_vertex[0] - left_vertex[0], right_vertex[1] - left_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto right_vertex_2d = origin + s_min * dir;
    right_vertex = glm::dvec3 { right_vertex_2d[0], right_vertex_2d[1], right_vertex[2] };
    KINDS_DEBUG(
      "Clamped right vertex to: (" << right_vertex[0] << ", " << right_vertex[1] << ", " << right_vertex[2] << ")");
  }
  else if (!left_inside && right_inside)
  {
    KINDS_DEBUG("Clamping left vertex: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
    // clamp left to the boundary
    glm::dvec2 origin { right_vertex[0], right_vertex[1] };
    glm::dvec2 dir { left_vertex[0] - right_vertex[0], left_vertex[1] - right_vertex[1] };
    auto hits = rayCast(boundary_points, origin, dir);

    double s_min = std::numeric_limits<double>::infinity();

    for (double s : hits)
    {
      if (s >= 0)
      {
        s_min = std::min(s_min, s);
      }
    }

    auto left_vertex_2d = origin + s_min * dir;
    left_vertex = glm::dvec3 { left_vertex_2d[0], left_vertex_2d[1], left_vertex[2] };
    KINDS_DEBUG(
      "Clamped left vertex to: (" << left_vertex[0] << ", " << left_vertex[1] << ", " << left_vertex[2] << ")");
  }
  else if (!left_inside && !right_inside)
  {
    // I'm not sure yet if this will work, especially for connections, but for now we will just discard it
    return false;
  }
  return true;
}

void kinDS::SegmentBuilder::finishMesh(size_t he_id, double t, const std::vector<BoundaryPoint>& boundary_points)
{
  size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
  auto& last_segments = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  if (last_segments.empty())
  {
    // this can happen if the first Voronoi vertex was discarded because it was outside of the boundary, in this case we
    // just don't create any triangles
    return;
  }

  // Get corresponding mesh
  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
  // Insert Voronoi vertex

  glm::dvec3 left_vertex = computeVoronoiVertex(he_id & ~1, t);
  glm::dvec3 right_vertex = computeVoronoiVertex((he_id & ~1) + 1, t);
  auto& he = kin_del.getGraph().getHalfEdges()[he_id & ~1];
  glm::dvec2 centroid = polygonCentroid(boundary_points);

  size_t v = he.origin;
  if (v == -1)
  {
    v = kin_del.getGraph().destination(he_id & ~1);
  }

  // TODO: Compute UVs here
  size_t new_left_vertex_index = addMeshletVertex(mesh, boundary_points, centroid, left_vertex, v, t);
  size_t new_right_vertex_index = addMeshletVertex(mesh, boundary_points, centroid, right_vertex, v, t);
  // build triangles

  for (auto& segment : last_segments)
  {

    size_t last_left = segment.mesh_start_vertex_id;
    size_t last_right = segment.mesh_end_vertex_id;
    // create two triangles
    // split quad differently depending on which side is closer
    if (last_left == last_right)
    {
      addMeshletTriangle(mesh, new_left_vertex_index, last_right, new_right_vertex_index);
    }
    else if (mesh.getVertices()[last_left][2] < mesh.getVertices()[last_right][2])
    {
      addMeshletTriangle(mesh, last_left, last_right, new_left_vertex_index);
      addMeshletTriangle(mesh, new_left_vertex_index, last_right, new_right_vertex_index);
    }
    else
    {
      addMeshletTriangle(mesh, last_left, last_right, new_right_vertex_index);
      addMeshletTriangle(mesh, last_left, new_right_vertex_index, new_left_vertex_index);
    }

    // update last vertex indices
    segment.mesh_start_vertex_id = static_cast<int>(new_left_vertex_index);
    segment.mesh_end_vertex_id = static_cast<int>(new_right_vertex_index);
  }
}

SegmentBuilder::SegmentBuilder(
  KineticDelaunay& kin_del, std::vector<std::pair<size_t, double>> subdivisions, bool create_transformed_mesh)
  : kin_del(kin_del)
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  assert(std::is_sorted(subdivisions.begin(), subdivisions.end(),
    [](const auto& a, const auto& b) { return a.second < b.second; }));
  kin_del.setSubdivisionSchedule(std::move(subdivisions));
}

SegmentBuilder::SegmentBuilder(KineticDelaunay& kin_del, bool create_transformed_mesh)
  : kin_del(kin_del)
  , section_callback_(std::make_unique<SegmentBuilderSectionCallback>(*this))
  , flip_callback_(std::make_unique<SegmentBuilderFlipCallback>(*this))
  , radius_callback_(std::make_unique<SegmentBuilderRadiusCallback>(*this))
  , crossing_callback_(std::make_unique<SegmentBuilderCrossingCallback>(*this))
  , subdivision_callback_(std::make_unique<SegmentBuilderSubdivisionCallback>(*this))
  , create_transformed_mesh(create_transformed_mesh)
  , boundary_mesh({ "bark", "interior" })
{
  kin_del.setSubdivisionSchedule({});
}

SegmentBuilder::~SegmentBuilder() = default;

static glm::dvec2 intersectSegments(
  const glm::dvec2& p, const glm::dvec2& p2, const glm::dvec2& q, const glm::dvec2& q2)
{
  glm::dvec2 r = p2 - p;
  glm::dvec2 s = q2 - q;

  double rxs = r.x * s.y - r.y * s.x;
  double qpxr = (q - p).x * r.y - (q - p).y * r.x;

  if (rxs == 0.0 && qpxr == 0.0)
    return glm::dvec2(std::numeric_limits<double>::quiet_NaN()); // collinear

  if (rxs == 0.0 && qpxr != 0.0)
    return glm::dvec2(std::numeric_limits<double>::infinity()); // parallel

  double t = ((q - p).x * s.y - (q - p).y * s.x) / rxs;
  double u = ((q - p).x * r.y - (q - p).y * r.x) / rxs;

  if (t < 0 || t > 1 || u < 0 || u > 1)
  {
    KINDS_WARNING("Segments do not intersect within the segment bounds");
  }

  return p + t * r;
}

void SegmentBuilder::startNewMesh(size_t half_edge_id, double t)
{
  size_t even_id = half_edge_id & ~1;
  size_t odd_id = even_id + 1;

  const auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[even_id];
  const auto& twin_he = graph.getHalfEdges()[odd_id];

  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
  segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();

  half_edge_index_to_segment_mesh_pair_index[even_id] = segment_mesh_pairs.size();
  half_edge_index_to_segment_mesh_pair_index[odd_id] = segment_mesh_pairs.size();

  segment_mesh_pairs.push_back(segment_mesh_pair);

  size_t vertex = std::max(he.origin, twin_he.origin);
  size_t component_id = kin_del.component_data.component_map[vertex];

  std::vector<bool> he_visited(graph.getHalfEdges().size(), false);
  updateBoundary(t, he_visited, component_id);

  auto& boundary_polygon = kin_del.component_data.component_boundaries[component_id][0];
  auto& centroid = kin_del.component_data.component_centroids[component_id];

  // For now also create a mesh, but this might be changed later
  VoronoiMesh mesh;

  glm::dvec3 left_vertex = computeVoronoiVertex(even_id, t);
  glm::dvec3 right_vertex = computeVoronoiVertex(odd_id, t);

  // Track how the Voronoi edge between these two vertices intersects with the boundary, we will need this information
  // to correctly build the boundary mesh
  size_t left_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[even_id].face;
  size_t left_containing_tri_id = kin_del.getCrossingDataContainingTriId(left_voronoi_vertex_id);

  size_t right_voronoi_vertex_id = kin_del.getGraph().getHalfEdges()[odd_id].face;

  // Now go through all faces
  bool inside = kin_del.getFaceInside(left_containing_tri_id);

  segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  // If already inside, the left Voronoi vertex is the first one to be added to the mesh
  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    size_t vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, left_vertex, he.origin, t);
    segment_mesh_pair_last_left_and_right_vertex.back().emplace_back(
      MeshingData { static_cast<int>(vertex_index), -1, -1, -1 });
  }

  if (kin_del.computeBoundaryOnTheFly())
  {
    // KINDS_DEBUG("Determining edge crossings of voronoi edge from " << left_voronoi_vertex_id << " to " <<
    // right_voronoi_vertex_id);
    // TODO: I don't think we need this anymore, we can use the edge intersection list instead
    std::vector<size_t> crossed_half_edges
      = kin_del.computeCrossedHalfEdges(left_containing_tri_id, right_vertex, left_vertex, t).first;
    // KINDS_DEBUG("Number of crossed half edges: " << crossed_half_edges.size() << " with voronoi vertex ids: " <<
    // left_voronoi_vertex_id << " and " << right_voronoi_vertex_id);
    for (size_t crossed_he_id : crossed_half_edges)
    {
      // KINDS_DEBUG("Crossed half edge: " << crossed_he_id);
      size_t prev_face_id = graph.getHalfEdges()[crossed_he_id].face;
      size_t next_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;

      bool next_inside = kin_del.getFaceInside(next_face_id);

      if (inside != next_inside)
      {
        // Compute the intersection point of the Voronoi edge with the boundary corresponding to this half-edge
        glm::dvec2 intersection_point
          = intersectSegments(glm::dvec2(left_vertex[0], left_vertex[1]), glm::dvec2(right_vertex[0], right_vertex[1]),
            kin_del.getPointAt(t, graph.getHalfEdges()[crossed_he_id].origin),
            kin_del.getPointAt(t, graph.getHalfEdges()[crossed_he_id ^ 1].origin));

        // Add the intersection point as a vertex to the mesh
        size_t vertex_index
          = addMeshletVertex(mesh, boundary_polygon, centroid, glm::dvec3(intersection_point, t), -1, t);

        if (inside)
        {
          // leaving the inside, add the intersection point as second of the last added pair
          segment_mesh_pair_last_left_and_right_vertex.back().back().mesh_end_vertex_id
            = static_cast<int>(vertex_index);
          segment_mesh_pair_last_left_and_right_vertex.back().back().end_half_edge_id = static_cast<int>(crossed_he_id);
        }
        else
        {
          // Use the twin half-edge as it is located inside the component.
          segment_mesh_pair_last_left_and_right_vertex.back().emplace_back(
            MeshingData { static_cast<int>(vertex_index), -1, static_cast<int>(crossed_he_id ^ 1), -1 });
        }
      }

      // current_face_id = next_face_id;
      inside = next_inside;
    }
  }

  if (inside || !kin_del.computeBoundaryOnTheFly())
  {
    // If we end inside, the right Voronoi vertex is the last one to be added to the mesh
    size_t vertex_index = addMeshletVertex(mesh, boundary_polygon, centroid, right_vertex, twin_he.origin, t);
    segment_mesh_pair_last_left_and_right_vertex.back().back().mesh_end_vertex_id = static_cast<int>(vertex_index);
  }

  meshes.push_back(mesh);

  assert(segment_mesh_pairs.size() == segment_mesh_pair_last_left_and_right_vertex.size());
}

void kinDS::SegmentBuilder::completeBoundaryMeshSection(size_t he_id, size_t new_left, size_t new_right)
{
  const auto& last_left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
  if (last_left_and_right.first != -1)
  {
    // distinguish the case that we have previously flipped an infinite edge that became a boundary edge
    if (half_edge_to_boundary_vertex_index[he_id] == -1)
    {
      // no edge flip
      addBoundaryTriangle(last_left_and_right.first, new_right, new_left);
      if (last_left_and_right.second != -1)
      {
        addBoundaryTriangle(last_left_and_right.second, new_right, last_left_and_right.first);
      }
    }
    else
    {
      assert(last_left_and_right.second != -1);
      // he_id was previously flipped, it's corresponding vertex is no longer part of the boundary
      addBoundaryTriangle(last_left_and_right.second, new_right, half_edge_to_boundary_vertex_index[he_id]);
      addBoundaryTriangle(new_left, last_left_and_right.first, half_edge_to_boundary_vertex_index[he_id]);
      addBoundaryTriangle(new_left, half_edge_to_boundary_vertex_index[he_id], new_right);

      // reset the half-edge to boundary vertex index
      half_edge_to_boundary_vertex_index[he_id] = -1;
    }
  }
  else
  {
    assert(last_left_and_right.second == -1);
  }
}

size_t kinDS::SegmentBuilder::addBoundaryTriangle(size_t u, size_t v, size_t w)
{
  // check bounds
  if (u >= boundary_mesh.getVertexCount() || v >= boundary_mesh.getVertexCount() || w >= boundary_mesh.getVertexCount())
  {
    KINDS_ERROR("Vertex index out of boundary mesh range.");
    return -1;
  }

  if (u >= boundary_mesh_raw_uvs.size() || v >= boundary_mesh_raw_uvs.size() || w >= boundary_mesh_raw_uvs.size())
  {
    KINDS_ERROR("Vertex index out of raw uv range.");
    return -1;
  }

  // get raw UVs
  glm::dvec3 uv_u = { boundary_mesh_raw_uvs[u][0], boundary_mesh_raw_uvs[u][1], 0.0 };
  glm::dvec3 uv_v = { boundary_mesh_raw_uvs[v][0], boundary_mesh_raw_uvs[v][1], 0.0 };
  glm::dvec3 uv_w = { boundary_mesh_raw_uvs[w][0], boundary_mesh_raw_uvs[w][1], 0.0 };

  // output UVs
  /*KINDS_DEBUG("Adding boundary triangle with raw UVs: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1])
     +
                "), v(" + std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) +
                ", " + std::to_string(uv_w[1]) + ")");*/

  // adjust UVs to avoid seams, first coordinate is the angle normalized to [-0.5, 0.5]
  // As a heuristic, we take the first angle and adjust the others such that they have less than 0.5 difference
  double base_angle = uv_u[0];
  double& angle_v = uv_v[0];
  double diff_v = angle_v - base_angle;
  double adjustment = std::round(diff_v);
  angle_v -= adjustment;

  double& angle_w = uv_w[0];
  double diff_w = angle_w - base_angle;
  adjustment = std::round(diff_w);
  angle_w -= adjustment;

  uv_u[0] *= uv_circum_factor;
  uv_v[0] *= uv_circum_factor;
  uv_w[0] *= uv_circum_factor;
  uv_u[1] *= uv_height_factor;
  uv_v[1] *= uv_height_factor;
  uv_w[1] *= uv_height_factor;

  // add adjusted UVs
  size_t uv_index_u = boundary_mesh.addUV(uv_u);
  size_t uv_index_v = boundary_mesh.addUV(uv_v);
  size_t uv_index_w = boundary_mesh.addUV(uv_w);

  /*KINDS_DEBUG("UVs after adjustment: u(" + std::to_string(uv_u[0]) + ", " + std::to_string(uv_u[1]) + "), v(" +
                std::to_string(uv_v[0]) + ", " + std::to_string(uv_v[1]) + "), w(" + std::to_string(uv_w[0]) + ", " +
                std::to_string(uv_w[1]) + ")");*/
  return boundary_mesh.addTriangle(u, v, w, uv_index_u, uv_index_v, uv_index_w, 0);
}

size_t kinDS::SegmentBuilder::addBoundaryVertex(glm::dvec3 vertex, glm::dvec2 centroid, size_t strand_id, double t)
{
  double angle = std::atan2(centroid[1] - vertex[1], centroid[0] - vertex[0]);

  glm::dvec2 raw_uv { angle / (2.0 * glm::pi<double>()), vertex[2] };
  if (create_transformed_mesh)
  {
    vertex = kin_del.getStrandTree().transformToObjectSpace(vertex, strand_id, t);
  }

  size_t index = boundary_mesh.addVertex(vertex);
  boundary_vertex_to_strand_id.push_back(strand_id);
  boundary_mesh_raw_uvs.resize(index + 1, glm::dvec2 {});
  boundary_mesh_raw_uvs[index] = raw_uv;

  return index;
}

size_t kinDS::SegmentBuilder::addMeshletTriangle(VoronoiMesh& mesh, size_t u, size_t v, size_t w)
{
  return mesh.addTriangle(u, v, w, u, v, w); // For meshlets, the UVs are assigned per vertex so the indices match
}

size_t kinDS::SegmentBuilder::addMeshletVertex(VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, glm::dvec3 vertex, size_t strand_id, double t)
{
  if (create_transformed_mesh)
  {
    vertex = kin_del.getStrandTree().transformToObjectSpace(vertex, strand_id, t);
  }
  size_t index = mesh.addVertex(vertex);
  double rel_dist = relativeDistanceFromCenter(boundary_polygon, centroid, glm::dvec2 { vertex[0], vertex[1] });

  /*if (rel_dist > 1.0 + std::numeric_limits<double>::epsilon()) {
    KINDS_WARNING("Adding vertex that is too far outside, relative distance: " << rel_dist);
  }*/

  // TODO: this can be simplified to not use trigonometric functions
  double angle = std::atan2(centroid[1] - vertex[1], centroid[0] - vertex[0]);
  double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
  double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
  mesh.addUV(u, v, vertex[2] * uv_height_factor);
  return index;
}

void kinDS::SegmentBuilder::addVoronoiTriangulationToBoundaryMesh(double t, bool invert_orientation, double offset)
{
  auto& graph = kin_del.getGraph();

  size_t index_offset = boundary_mesh.getVertices().size();
  std::vector<double> relative_center_distances;
  // add all vertices
  for (size_t i = 0; i < graph.getVertexCount(); i++)
  {
    glm::dvec2 vertex = kin_del.getPointAt(t, i);

    auto component_index = kin_del.component_data.component_map[i];
    auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    addBoundaryVertex(glm::dvec3 { vertex[0], vertex[1], t + offset }, centroid, i, t);
    // KINDS_DEBUG("New raw uv: " << raw_uv[0] << ", " << raw_uv[1] << " for vertex: " << vertex_index);

    double rel_dist = relativeDistanceFromCenter(boundary_polygon, centroid, vertex);
    relative_center_distances.push_back(rel_dist);
  }
  // add all triangles
  // for (const auto& triangle : graph.getFaces()) {
  for (size_t face_index = 0; face_index < graph.getFaces().size(); face_index++)
  {
    const auto& triangle = graph.getFaces()[face_index];
    const auto& he_ids = triangle.half_edges;
    auto vertices = graph.adjacentTriangleVertices(triangle.half_edges[0]);

    // check if on component boundary and update last left and right vertices accordingly
    // store last left and right
    for (size_t i = 0; i < 3; i++)
    {
      if (kin_del.isOnComponentBoundaryOutside(he_ids[i]))
      {
        completeBoundaryMeshSection(he_ids[i], index_offset + vertices[i], index_offset + vertices[(i + 1) % 3]);
        // KINDS_DEBUG("Assigning boundary last vertices for he_id " << he_ids[i]);
        boundary_mesh_last_left_and_right_vertex[he_ids[i]]
          = std::make_pair(index_offset + vertices[i], index_offset + vertices[(i + 1) % 3]);
      }
    }
    // skip faces that are outside
    if (!kin_del.getFaceInside(face_index))
    {
      continue;
    }

    // check for infinite vertices
    if (vertices[0] == -1 || vertices[1] == -1 || vertices[2] == -1)
    {
      continue; // skip triangles with infinite vertices
    }

    if (invert_orientation)
    {
      std::swap(vertices[1], vertices[2]);
    }

    // next, get the angles from the raw UVs and convert back to cartesian coordinates centered at (0.5, 0.5)
    size_t uv_indices[3];

    for (size_t i = 0; i < 3; i++)
    {
      double rel_dist = relative_center_distances[vertices[i]];
      double angle = boundary_mesh_raw_uvs[index_offset + vertices[i]][0] * 2.0 * glm::pi<double>();
      double u = 0.5 + texture_diameter * rel_dist * 0.5 * std::cos(angle);
      double v = 0.5 + texture_diameter * rel_dist * 0.5 * std::sin(angle);
      uv_indices[i] = boundary_mesh.addUV(u, v, 0.0);
    }

    // as an exception, we directly add the triangle here to have access to the UV indices
    boundary_mesh.addTriangle(index_offset + vertices[0], index_offset + vertices[1], index_offset + vertices[2],
      uv_indices[0], uv_indices[1], uv_indices[2], 1);
  }
}

std::vector<BoundaryPoint> kinDS::SegmentBuilder::traceConvexHull(double t) const
{
  const auto& graph = kin_del.getGraph();
  std::vector<BoundaryPoint> convex_hull_points;
  for (HalfEdgeDelaunayGraph::ConvexHullEdgeIterator it = graph.boundaryEdgesBegin(), end = graph.boundaryEdgesEnd();
    it != end; ++it)
  {
    size_t he_id = *it;
    size_t strand_index = graph.getHalfEdges()[he_id].origin;

    glm::dvec2 convex_hull_point = kin_del.getPointAt(t, strand_index);
    convex_hull_points.push_back({ strand_index, he_id, convex_hull_point });
  }

  return convex_hull_points;
}

void kinDS::SegmentBuilder::advanceBoundaryMesh(
  double t, const std::vector<BoundaryPoint>& boundary_points, const glm::dvec2& centroid)
{
  std::vector<size_t> new_vertex_indices;

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    glm::dvec2 boundary_point = boundary_points[i].p;

    new_vertex_indices.push_back(addBoundaryVertex(
      glm::dvec3 { boundary_point[0], boundary_point[1], t }, centroid, boundary_points[i].vertex_id, t));
  }

  for (size_t i = 0; i < boundary_points.size(); i++)
  {
    size_t he_id = boundary_points[i].he_id;
    size_t left_vertex_index = new_vertex_indices[i];
    size_t right_vertex_index = new_vertex_indices[(i + 1) % boundary_points.size()];
    auto& left_and_right = boundary_mesh_last_left_and_right_vertex[he_id];
    completeBoundaryMeshSection(he_id, left_vertex_index, right_vertex_index);

    // KINDS_DEBUG("Assigning boundary last vertices for he_id " << he_id);
    left_and_right.first = left_vertex_index;
    left_and_right.second = right_vertex_index;
  }
}

void kinDS::SegmentBuilder::updateBoundary(double t, std::vector<bool>& visited, size_t component_index)
{
  if (kin_del.component_data.component_last_updated[component_index] != t)
  {
    kin_del.component_data.component_boundaries[component_index]
      = kin_del.extractComponentBoundaries(kin_del.component_data.components[component_index], t, visited);
    kin_del.component_data.component_centroids[component_index]
      = polygonCentroid(kin_del.component_data.component_boundaries[component_index][0]);
    kin_del.component_data.component_last_updated[component_index] = t;
  }
}

void kinDS::SegmentBuilder::updateBoundaries(double t)
{
  std::vector<bool> visited(kin_del.getGraph().getHalfEdges().size(), false);

  for (size_t component_index = 0; component_index < kin_del.component_data.components.size(); component_index++)
  {
    updateBoundary(t, visited, component_index);
  }
}

void kinDS::SegmentBuilder::advanceBoundaryMeshes(double t)
{
  updateBoundaries(t);

  for (size_t component_index = 0; component_index < kin_del.component_data.components.size(); component_index++)
  {
    auto& boundaries = kin_del.component_data.component_boundaries[component_index];
    auto& centroid = kin_del.component_data.component_centroids[component_index];

    for (auto& boundary_points : boundaries)
    {
      advanceBoundaryMesh(t, boundary_points, centroid);
    }
  }
}

std::vector<glm::dvec3> SegmentBuilder::computeClampedVoronoiVertices(
  size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid)
{
  auto& graph = kin_del.getGraph();
  std::vector<glm::dvec3> voronoi_vertices;

  // we just create a triangle fan because the Voronoi cell is convex
  // iterate over all segment indices of the strand
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    glm::dvec3 voronoi_vertex = computeVoronoiVertex(*it, t);
    voronoi_vertices.emplace_back(voronoi_vertex);
    // addMeshletVertex(mesh, boundary_polygon, centroid, voronoi_vertex);
  }

  std::vector<glm::dvec3> clamped_voronoi_vertices;
  // go through the vertices and clamp them
  for (size_t i = 0; i < voronoi_vertices.size(); i++)
  {
    auto left = voronoi_vertices[i];
    auto right = voronoi_vertices[(i + 1) % voronoi_vertices.size()];

    bool inside = clampVoronoiVertices(left, right, boundary_polygon, centroid);

    if (inside)
    {
      clamped_voronoi_vertices.push_back(left);
      clamped_voronoi_vertices.push_back(right);
    }
  }

  voronoi_vertices.clear();
  // Put back into original vector and remove duplicates
  voronoi_vertices.emplace_back(clamped_voronoi_vertices.front());
  for (size_t i = 1; i < clamped_voronoi_vertices.size() - 1; i++)
  {
    if (clamped_voronoi_vertices[i] != voronoi_vertices.back())
    {
      voronoi_vertices.push_back(clamped_voronoi_vertices[i]);
    }
  }
  if (clamped_voronoi_vertices.front() != clamped_voronoi_vertices.back())
  {
    voronoi_vertices.push_back(clamped_voronoi_vertices.back());
  }

  return voronoi_vertices;
}

size_t kinDS::SegmentBuilder::closingMeshCountStrandIncidentEdges(size_t strand_id) const
{
  auto& graph = kin_del.getGraph();
  size_t num = 0;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++num;
  }
  return num;
}

int kinDS::SegmentBuilder::closingMeshAppendVertex(VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, size_t strand_id, double t,
  const glm::dvec3& position)
{
  const size_t vertex_id = addMeshletVertex(mesh, boundary_polygon, centroid, position, strand_id, t);
  return static_cast<int>(vertex_id);
}

const void* kinDS::SegmentBuilder::closingMeshFindVoronoiEdgeIntersectionPtr(
  size_t voronoi_edge_id, size_t crossed_delaunay_he_id) const
{
  const size_t d = crossed_delaunay_he_id / 2;
  const auto& crossing_data = kin_del.getCrossingData();
  if (voronoi_edge_id >= crossing_data.voronoi_edge_intersections.size())
  {
    return nullptr;
  }
  for (const auto& ref : crossing_data.voronoi_edge_intersections[voronoi_edge_id])
  {
    if (ref->delaunay_edge_id == d)
    {
      return &(*ref);
    }
  }
  return nullptr;
}

glm::dvec3 kinDS::SegmentBuilder::closingMeshVoronoiDelaunayCrossingPosition(
  double t, size_t voronoi_edge_id, size_t delaunay_edge_id) const
{
  auto& graph = kin_del.getGraph();
  const size_t voronoi_he0 = 2 * voronoi_edge_id;
  const size_t voronoi_he1 = voronoi_he0 + 1;

  const glm::dvec3 left_vertex = computeVoronoiVertex(voronoi_he0, t);
  const glm::dvec3 right_vertex = computeVoronoiVertex(voronoi_he1, t);
  const glm::dvec2 left2(left_vertex.x, left_vertex.y);
  const glm::dvec2 right2(right_vertex.x, right_vertex.y);

  const size_t d_he0 = 2 * delaunay_edge_id;
  const size_t d_he1 = d_he0 + 1;
  const size_t d_a = graph.getHalfEdges()[d_he0].origin;
  const size_t d_b = graph.getHalfEdges()[d_he1].origin;
  if (d_a == size_t(-1) || d_b == size_t(-1))
  {
    return glm::dvec3(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), t);
  }

  const glm::dvec2 d0 = kin_del.getPointAt(t, d_a);
  const glm::dvec2 d1 = kin_del.getPointAt(t, d_b);
  const glm::dvec2 p = intersectSegments(left2, right2, d0, d1);
  return glm::dvec3(p, t);
}

auto kinDS::SegmentBuilder::extractSegmentsForVoronoiEdge(double t, int incident_edge_index, size_t voronoi_edge_id,
  const std::function<int(const glm::dvec3&)>& track_vertex, bool reverse) -> std::vector<MeshingData>
{
  std::vector<MeshingData> out_segments;
  auto& graph = kin_del.getGraph();

  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return out_segments;
  }

  // `reverse`: when false, walk even circumcenter → odd; when true, odd → even (strand-side-first in closingMesh).
  const size_t voronoi_he_strand_side = reverse ? voronoi_he_odd : voronoi_he_even;
  const size_t voronoi_he_other_side = voronoi_he_strand_side ^ 1;
  const bool strand_at_voronoi_even_he = !reverse;

  const glm::dvec3 strand_cm_vertex = computeVoronoiVertex(voronoi_he_strand_side, t);
  const glm::dvec3 other_cm_vertex = computeVoronoiVertex(voronoi_he_other_side, t);
  const glm::dvec2 strand2 = glm::dvec2(strand_cm_vertex);
  const glm::dvec2 other2 = glm::dvec2(other_cm_vertex);
  const glm::dvec2 vor_dir = other2 - strand2;
  const double vor_len2 = glm::dot(vor_dir, vor_dir);
  if (vor_len2 < 1e-14)
  {
    return out_segments;
  }

  const size_t strand_face_id = graph.getHalfEdges()[voronoi_he_strand_side].face;
  const size_t strand_containing_tri_id = kin_del.getCrossingDataContainingTriId(strand_face_id);
  bool inside = kin_del.getFaceInside(strand_containing_tri_id);
  if (inside)
  {
    out_segments.push_back(MeshingData { track_vertex(strand_cm_vertex), -1, -1, -1 });
    out_segments.back().closing_incident_edge_index = incident_edge_index;
    out_segments.back().closing_voronoi_edge_id = voronoi_edge_id;
    out_segments.back().closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
  }

  if (!kin_del.computeBoundaryOnTheFly())
  {
    if (inside && !out_segments.empty())
    {
      out_segments.back().mesh_end_vertex_id = track_vertex(other_cm_vertex);
    }
    return out_segments;
  }

  const auto crossed_edges_params = kin_del.computeCrossedHalfEdges(strand_containing_tri_id, other2, strand2, t);
  for (size_t crossed_he_id : crossed_edges_params.first)
  {
    const size_t next_face_id = graph.getHalfEdges()[crossed_he_id ^ 1].face;
    const bool next_inside = kin_del.getFaceInside(next_face_id);

    if (inside != next_inside)
    {
      const size_t da = graph.getHalfEdges()[crossed_he_id].origin;
      const size_t db = graph.getHalfEdges()[crossed_he_id ^ 1].origin;
      if (da != size_t(-1) && db != size_t(-1))
      {
        const glm::dvec2 d0 = kin_del.getPointAt(t, da);
        const glm::dvec2 d1 = kin_del.getPointAt(t, db);
        const glm::dvec2 p = intersectSegments(strand2, other2, d0, d1);

        if (std::isfinite(p.x) && std::isfinite(p.y))
        {
          const int mesh_vertex_id = track_vertex(glm::dvec3(p, t));
          const void* intersection_ptr = closingMeshFindVoronoiEdgeIntersectionPtr(voronoi_edge_id, crossed_he_id);

          if (inside)
          {
            if (!out_segments.empty())
            {
              auto& segment = out_segments.back();
              segment.mesh_end_vertex_id = mesh_vertex_id;
              segment.end_half_edge_id = static_cast<int>(crossed_he_id);
              segment.end_intersection_ref = intersection_ptr;
            }
          }
          else
          {
            MeshingData segment { mesh_vertex_id, -1, static_cast<int>(crossed_he_id ^ 1), -1 };
            segment.start_intersection_ref = intersection_ptr;
            segment.closing_incident_edge_index = incident_edge_index;
            segment.closing_voronoi_edge_id = voronoi_edge_id;
            segment.closing_strand_at_voronoi_even_he = strand_at_voronoi_even_he;
            out_segments.push_back(segment);
          }
        }
      }
    }

    inside = next_inside;
  }

  if (inside && !out_segments.empty())
  {
    out_segments.back().mesh_end_vertex_id = track_vertex(other_cm_vertex);
  }

  return out_segments;
}

auto kinDS::SegmentBuilder::closingMeshExtractRawSegmentsForVoronoiEdge(size_t strand_id, double t, VoronoiMesh& mesh,
  const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid, int incident_edge_index,
  size_t incident_he, const std::function<int(const glm::dvec3&)>& track_vertex) -> std::vector<MeshingData>
{
  auto& graph = kin_del.getGraph();

  // Undirected Delaunay edge id (matches `CrossingData` / `computeEdgeIntersections`: voronoi_edge_id == he/2).
  const size_t voronoi_edge_id = incident_he / 2;
  const size_t voronoi_he_even = 2 * voronoi_edge_id;
  const size_t voronoi_he_odd = voronoi_he_even + 1;

  if (graph.isInfinite(voronoi_he_even))
  {
    return {};
  }

  const int origin_even = graph.getHalfEdges()[voronoi_he_even].origin;
  const int origin_odd = graph.getHalfEdges()[voronoi_he_odd].origin;
  const auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
  if (!finite_origin_is_strand(origin_even) && !finite_origin_is_strand(origin_odd))
  {
    KINDS_ERROR("createClosingMesh: Voronoi edge for dual Delaunay edge_id "
      << voronoi_edge_id << " (half-edges " << voronoi_he_even << ", " << voronoi_he_odd << ", incident_he "
      << incident_he << ") has no finite half-edge with origin strand_id " << strand_id << " (origins " << origin_even
      << ", " << origin_odd << ")");
  }

  // Orient the inside trace from the Voronoi half-edge whose origin is the strand (matches how we leave the strand
  // along the dual). `CrossingData::computeEdgeIntersections` still stores the Voronoi list in even->odd order; when
  // the strand is on the odd half-edge we reverse traversal of `delaunay_edge_intersections` during the boundary
  // walk.
  const bool strand_at_voronoi_even_he = finite_origin_is_strand(origin_even);
  const bool reverse = !strand_at_voronoi_even_he;

  return extractSegmentsForVoronoiEdge(t, incident_edge_index, voronoi_edge_id, track_vertex, reverse);
}

SegmentBuilder::ClosingMeshOrderedIndex kinDS::SegmentBuilder::closingMeshBuildOrderedIndex(
  std::list<MeshingData>& closing_segments)
{
  ClosingMeshOrderedIndex out;
  out.ordered_segments.reserve(closing_segments.size());
  for (auto& seg : closing_segments)
  {
    if (seg.mesh_start_vertex_id != -1 && seg.mesh_end_vertex_id != -1)
    {
      out.ordered_segments.push_back(&seg);
    }
  }

  out.start_ref_to_segment.reserve(out.ordered_segments.size());
  for (size_t i = 0; i < out.ordered_segments.size(); ++i)
  {
    if (out.ordered_segments[i]->start_intersection_ref != nullptr)
    {
      out.start_ref_to_segment[out.ordered_segments[i]->start_intersection_ref] = i;
    }
  }
  return out;
}

void kinDS::SegmentBuilder::closingMeshLogDiscoveredSegments(size_t strand_id, double t,
  const std::list<MeshingData>& closing_segments, const std::vector<MeshingData*>& ordered_segments)
{
  KINDS_INFO("createClosingMesh strand " << strand_id << " t=" << t << ": " << closing_segments.size()
                                         << " raw closing segment(s), " << ordered_segments.size()
                                         << " ordered (complete mesh vertex ids).");
  size_t skipped_incomplete = 0;
  for (const auto& seg : closing_segments)
  {
    if (seg.mesh_start_vertex_id == -1 || seg.mesh_end_vertex_id == -1)
    {
      ++skipped_incomplete;
      KINDS_INFO("  skipped incomplete closing segment: incident_edge="
        << seg.closing_incident_edge_index << " voronoi_edge_id="
        << (seg.closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(seg.closing_voronoi_edge_id))
        << " mesh_v(" << seg.mesh_start_vertex_id << "->" << seg.mesh_end_vertex_id << ") he(" << seg.start_half_edge_id
        << "->" << seg.end_half_edge_id
        << ") start_ref=" << kin_del.formatCrossingIntersectionForLog(seg.start_intersection_ref)
        << " end_ref=" << kin_del.formatCrossingIntersectionForLog(seg.end_intersection_ref));
    }
  }
  if (skipped_incomplete > 0)
  {
    KINDS_INFO("  (" << skipped_incomplete << " incomplete segment(s) omitted from ordered_segments / matching.)");
  }
  for (size_t i = 0; i < ordered_segments.size(); ++i)
  {
    const MeshingData* s = ordered_segments[i];
    KINDS_INFO("  ordered Voronoi segment ["
      << i << "] incident_edge=" << s->closing_incident_edge_index << " voronoi_edge_id="
      << (s->closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(s->closing_voronoi_edge_id))
      << " dual_delaunay_edge_id="
      << (s->closing_voronoi_edge_id == static_cast<size_t>(-1) ? -1 : static_cast<int>(s->closing_voronoi_edge_id))
      << " mesh_v(" << s->mesh_start_vertex_id << "->" << s->mesh_end_vertex_id << ") he(" << s->start_half_edge_id
      << "->" << s->end_half_edge_id << ") strand_at_voronoi_even_he=" << (s->closing_strand_at_voronoi_even_he ? 1 : 0)
      << " start_ref=" << kin_del.formatCrossingIntersectionForLog(s->start_intersection_ref)
      << " end_ref=" << kin_del.formatCrossingIntersectionForLog(s->end_intersection_ref));
  }
}

void kinDS::SegmentBuilder::closingMeshValidateOrderedSegmentGeometry(
  double t, const VoronoiMesh& mesh, const std::vector<MeshingData*>& ordered_segments)
{
  constexpr double k_closing_cap_geom_eps = 1e-5;
  auto mesh_vertex_xy = [&](int mesh_vid) -> glm::dvec2
  {
    const glm::dvec3& v = mesh.getVertices()[static_cast<size_t>(mesh_vid)];
    return glm::dvec2(v.x, v.y);
  };

  for (size_t si = 0; si < ordered_segments.size(); ++si)
  {
    const MeshingData* s = ordered_segments[si];
    if (s->closing_voronoi_edge_id == static_cast<size_t>(-1))
    {
      continue;
    }
    const size_t ve = s->closing_voronoi_edge_id;
    {
      const std::string ctxs = std::string("ordered_seg[") + std::to_string(si)
        + "] incident=" + std::to_string(s->closing_incident_edge_index);
      kin_del.validateClosingCapCrossingRef(
        (ctxs + " start_ref").c_str(), s->start_intersection_ref, ve, s->start_half_edge_id);
      kin_del.validateClosingCapCrossingRef(
        (ctxs + " end_ref").c_str(), s->end_intersection_ref, ve, s->end_half_edge_id);
    }

    const glm::dvec3 L3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve, t);
    const glm::dvec3 R3 = kin_del.computeVoronoiVertexClampedInfinity(2 * ve + 1, t);
    const glm::dvec2 L(L3.x, L3.y);
    const glm::dvec2 R(R3.x, R3.y);
    const glm::dvec2 axis = R - L;
    const double ax2 = glm::dot(axis, axis);
    if (ax2 < 1e-24)
    {
      KINDS_ERROR("createClosingMesh ordered_seg[" << si << "]: degenerate Voronoi edge " << ve << " (zero length).");
      continue;
    }
    const glm::dvec2 ps = mesh_vertex_xy(s->mesh_start_vertex_id);
    const glm::dvec2 pe = mesh_vertex_xy(s->mesh_end_vertex_id);
    const double axis_len = glm::length(axis);
    const auto perp_dist
      = [&](const glm::dvec2& p) -> double { return std::abs((p.x - L.x) * axis.y - (p.y - L.y) * axis.x) / axis_len; };
    if (perp_dist(ps) > k_closing_cap_geom_eps || perp_dist(pe) > k_closing_cap_geom_eps)
    {
      KINDS_ERROR("createClosingMesh ordered_seg["
        << si << "]: mesh_start/end not collinear with canonical Voronoi "
        << "edge " << ve << " (even->odd circumcenters); check segment direction vs inside walk.");
    }
    const double ts = glm::dot(ps - L, axis) / ax2;
    const double te = glm::dot(pe - L, axis) / ax2;
    // mesh_start is always the strand-side circumcenter; along the canonical even->odd axis that is the smaller t iff
    // the strand Voronoi half-edge is even.
    if (s->closing_strand_at_voronoi_even_he)
    {
      if (ts > te + 1e-8)
      {
        KINDS_ERROR("createClosingMesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                                     << " (t_start=" << ts << " > t_end=" << te
                                                     << "); expected mesh_start nearer voronoi_he_even.");
      }
    }
    else if (te > ts + 1e-8)
    {
      KINDS_ERROR("createClosingMesh ordered_seg[" << si << "]: mesh vertices reversed along Voronoi edge " << ve
                                                   << " (t_start=" << ts << " > t_end=" << te
                                                   << " along even->odd); expected mesh_start nearer "
                                                   << "voronoi_he_odd (strand-side).");
    }

    glm::dvec2 ip;
    if (s->start_intersection_ref != nullptr
      && kin_del.tryComputeCrossingIntersectionPosition2D(s->start_intersection_ref, t, ip))
    {
      if (glm::distance(ip, ps) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("createClosingMesh ordered_seg[" << si << "]: start_ref 2D position does not match mesh_start.");
      }
    }
    if (s->end_intersection_ref != nullptr
      && kin_del.tryComputeCrossingIntersectionPosition2D(s->end_intersection_ref, t, ip))
    {
      if (glm::distance(ip, pe) > k_closing_cap_geom_eps)
      {
        KINDS_ERROR("createClosingMesh ordered_seg[" << si << "]: end_ref 2D position does not match mesh_end.");
      }
    }
  }
}

SegmentBuilder::ClosingMeshPolygonsTraceResult kinDS::SegmentBuilder::closingMeshTraceCapPolygons(size_t strand_id,
  double t, size_t num_incident_edges, VoronoiMesh& mesh, const std::vector<BoundaryPoint>& boundary_polygon,
  const glm::dvec2& centroid, std::vector<size_t> mesh_vertex_ids, const std::vector<MeshingData*>& ordered_segments,
  const std::unordered_map<const void*, size_t>& start_ref_to_segment)
{
  ClosingMeshPolygonsTraceResult result;
  result.mesh_vertex_ids = std::move(mesh_vertex_ids);
  result.segment_used.assign(ordered_segments.size(), false);

  auto& graph = kin_del.getGraph();
  const auto& crossing_data = kin_del.getCrossingData();

  auto& cap_vertex_ids = result.mesh_vertex_ids;
  auto add_mesh_vertex = [&](const glm::dvec3& v) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, v);
    cap_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  // Only connect along a Delaunay boundary edge to intersections whose Voronoi edge is dual to an edge incident to
  // this strand; otherwise the first CrossingData hit on that Delaunay edge can belong to another cell and we
  // incorrectly skip the strand's Voronoi continuation.
  auto voronoi_edge_incident_to_strand = [&](size_t voronoi_edge_id) -> bool
  {
    const size_t he0 = 2 * voronoi_edge_id;
    const size_t he1 = he0 + 1;
    if (he1 >= graph.getHalfEdges().size())
    {
      return false;
    }
    const int o0 = graph.getHalfEdges()[he0].origin;
    const int o1 = graph.getHalfEdges()[he1].origin;
    auto finite_origin_is_strand = [&](int o) { return o >= 0 && static_cast<size_t>(o) == strand_id; };
    return finite_origin_is_strand(o0) || finite_origin_is_strand(o1);
  };

  auto closing_cap_edge_ids = [](const MeshingData* s) -> std::string
  {
    if (s->closing_voronoi_edge_id != static_cast<size_t>(-1))
    {
      const size_t ve = s->closing_voronoi_edge_id;
      return std::string("voronoi_edge_id=") + std::to_string(ve) + " dual_delaunay_edge_id=" + std::to_string(ve);
    }
    return std::string("voronoi_edge_id=? dual_delaunay_edge_id=?");
  };

  // For a single vertex append: only the ref that corresponds to that endpoint of the segment.
  auto closing_cap_vertex_detail = [&](const MeshingData* s, bool is_segment_mesh_start) -> std::string
  {
    if (is_segment_mesh_start)
    {
      return closing_cap_edge_ids(s)
        + " start_ref=" + kin_del.formatCrossingIntersectionForLog(s->start_intersection_ref);
    }
    return closing_cap_edge_ids(s) + " end_ref=" + kin_del.formatCrossingIntersectionForLog(s->end_intersection_ref);
  };

  size_t closing_polygon_walk_seq = 0;
  for (;;)
  {
    size_t seed_index = static_cast<size_t>(-1);
    for (size_t i = 0; i < ordered_segments.size(); ++i)
    {
      if (!result.segment_used[i])
      {
        seed_index = i;
        break;
      }
    }
    if (seed_index == static_cast<size_t>(-1))
    {
      break;
    }

    const size_t walk_id = closing_polygon_walk_seq++;
    std::vector<size_t> polygon;
    size_t walk_vertex_ordinal = 0;
    auto walk_append = [&](const char* category, const std::string& detail, size_t vertex_id)
    {
      if (polygon.empty() || polygon.back() != vertex_id)
      {
        polygon.push_back(vertex_id);
        const glm::dvec3& v = mesh.getVertices()[vertex_id];
        KINDS_INFO("createClosingMesh strand " << strand_id << " t=" << t << " walk " << walk_id << " add_vertex["
                                               << walk_vertex_ordinal++ << "] " << category << " mesh_id=" << vertex_id
                                               << " pos=(" << v.x << "," << v.y << "," << v.z << ") " << detail);
      }
    };

    MeshingData* seed_seg = ordered_segments[seed_index];
    KINDS_INFO("createClosingMesh strand "
      << strand_id << " t=" << t << " walk " << walk_id << " START from ordered_seg[" << seed_index
      << "] incident_edge=" << seed_seg->closing_incident_edge_index << " mesh_v(" << seed_seg->mesh_start_vertex_id
      << "->" << seed_seg->mesh_end_vertex_id << ")"
      << " he(start=" << seed_seg->start_half_edge_id << " end=" << seed_seg->end_half_edge_id << ")");

    size_t current_segment_index = seed_index;
    size_t guard = 0;

    while (guard++ < ordered_segments.size() * 4)
    {
      bool closed_walk_at_seed = false;
      MeshingData* current_segment = ordered_segments[current_segment_index];
      result.segment_used[current_segment_index] = true;

      walk_append("voronoi_inside",
        std::string("ordered_seg[") + std::to_string(current_segment_index) + "] along_inside_voronoi START; "
          + closing_cap_vertex_detail(current_segment, true),
        static_cast<size_t>(current_segment->mesh_start_vertex_id));
      walk_append("voronoi_inside",
        std::string("ordered_seg[") + std::to_string(current_segment_index)
          + "] along_inside_voronoi END (to boundary stitch); " + closing_cap_vertex_detail(current_segment, false),
        static_cast<size_t>(current_segment->mesh_end_vertex_id));

      const bool incomplete_at_boundary_exit
        = (current_segment->end_intersection_ref == nullptr || current_segment->end_half_edge_id == -1);
      if (incomplete_at_boundary_exit)
      {
        // End of inside polyline on this Voronoi edge (free endpoint at circumcenter): continue along the strand cap
        // to the next incident Voronoi edge without requiring a boundary intersection on this edge.
        if (current_segment->closing_incident_edge_index < 0 || num_incident_edges == 0)
        {
          break;
        }
        const int next_edge = (current_segment->closing_incident_edge_index + 1) % static_cast<int>(num_incident_edges);
        size_t next_segment_index = static_cast<size_t>(-1);
        for (size_t j = 0; j < ordered_segments.size(); ++j)
        {
          if (result.segment_used[j])
          {
            continue;
          }
          if (ordered_segments[j]->closing_incident_edge_index == next_edge)
          {
            next_segment_index = j;
            break;
          }
        }
        if (next_segment_index == static_cast<size_t>(-1))
        {
          break;
        }
        if (next_segment_index == seed_index)
        {
          break;
        }
        current_segment_index = next_segment_index;
        continue;
      }

      const int exit_he_id = current_segment->end_half_edge_id;
      const void* exit_intersection_ref = current_segment->end_intersection_ref;

      const size_t start_boundary_he = kin_del.isOnComponentBoundaryOutside(static_cast<size_t>(exit_he_id))
        ? static_cast<size_t>(exit_he_id)
        : static_cast<size_t>(exit_he_id) ^ 1;

      size_t boundary_he = start_boundary_he;
      const void* current_ref_ptr = exit_intersection_ref;
      bool found_next_segment = false;
      size_t boundary_guard = 0;

      constexpr double k_strand_corner_dist_eps = 1e-9;
      auto append_strand_corner_if_needed = [&]()
      {
        const glm::dvec2 xy = kin_del.getPointAt(t, strand_id);
        const glm::dvec3 strand_pos(xy, t);
        if (!polygon.empty())
        {
          const glm::dvec3& last = mesh.getVertices()[polygon.back()];
          if (glm::distance(glm::dvec2(last.x, last.y), xy) <= k_strand_corner_dist_eps)
          {
            return;
          }
        }
        const int vid = add_mesh_vertex(strand_pos);
        walk_append("corner_strand",
          std::string("strand_id=") + std::to_string(strand_id) + " boundary_he=" + std::to_string(boundary_he)
            + " delaunay_edge_id=" + std::to_string(boundary_he / 2) + " (Delaunay boundary corner at strand)",
          static_cast<size_t>(vid));
      };

      while (boundary_guard++ < graph.getHalfEdges().size() * 2)
      {
        const int he_origin = graph.getHalfEdges()[boundary_he].origin;
        if (he_origin >= 0 && static_cast<size_t>(he_origin) == strand_id)
        {
          append_strand_corner_if_needed();
        }

        const auto& d_intersections = crossing_data.delaunay_edge_intersections[boundary_he / 2];
        const bool forward_on_delaunay = (boundary_he % 2 == 0);
        // Crossing lists are ordered for the even Delaunay half-edge; flip when this cap segment's strand attaches on
        // the odd Voronoi half-edge (dual walk is strand-side -> other, opposite of even->odd storage).
        const bool effective_list_forward = forward_on_delaunay ^ !current_segment->closing_strand_at_voronoi_even_he;

        // Walk every crossing along this directed Delaunay half-edge after `current_ref_ptr`: each crossing gets a
        // mesh vertex. Strand-incident crossings that match an ordered segment start_ref hand off to that segment;
        // others only advance along the same Delaunay edge (no skipping to the corner).
        bool exit_boundary_chain_to_voronoi = false;
        if (!d_intersections.empty())
        {
          auto find_next_intersection = [&](const void* ref_ptr) -> decltype(d_intersections.begin())
          {
            if (d_intersections.empty())
            {
              return d_intersections.end();
            }
            if (ref_ptr == nullptr)
            {
              if (effective_list_forward)
              {
                return d_intersections.begin();
              }
              return std::prev(d_intersections.end());
            }
            for (auto it = d_intersections.begin(); it != d_intersections.end(); ++it)
            {
              if (&(*(*it)) == ref_ptr)
              {
                if (effective_list_forward)
                {
                  return std::next(it);
                }
                if (it == d_intersections.begin())
                {
                  return d_intersections.end();
                }
                return std::prev(it);
              }
            }
            return d_intersections.end();
          };

          for (size_t hop = 0; hop <= d_intersections.size(); ++hop)
          {
            const auto next_it = find_next_intersection(current_ref_ptr);
            if (next_it == d_intersections.end())
            {
              break;
            }

            const auto& candidate_ref = *next_it;
            const void* cand_ptr = &(*candidate_ref);
            const std::string cross_base = std::string("crossing=") + kin_del.formatCrossingIntersectionForLog(cand_ptr)
              + " boundary_he=" + std::to_string(boundary_he) + " delaunay_edge(boundary_he/2)="
              + std::to_string(boundary_he / 2) + " forward_on_delaunay=" + (forward_on_delaunay ? "1" : "0")
              + " effective_list_forward=" + (effective_list_forward ? "1" : "0");

            MeshingData* seed_ordered_seg = ordered_segments[seed_index];
            const bool strand_inc = voronoi_edge_incident_to_strand(candidate_ref->voronoi_edge_id);
            const auto seg_it_start = strand_inc ? start_ref_to_segment.find(cand_ptr) : start_ref_to_segment.end();
            size_t next_segment_index = static_cast<size_t>(-1);
            if (seg_it_start != start_ref_to_segment.end())
            {
              next_segment_index = seg_it_start->second;
            }
            else if (strand_inc && cand_ptr == seed_ordered_seg->end_intersection_ref)
            {
              // Completing the loop can return to the seed at its boundary end without a start_ref map entry.
              next_segment_index = seed_index;
            }
            const bool closes_polygon_at_seed_crossing = (next_segment_index == seed_index)
              && (cand_ptr == seed_ordered_seg->start_intersection_ref
                || cand_ptr == seed_ordered_seg->end_intersection_ref);

            // Walking a new boundary half-edge resets current_ref_ptr to null; the first listed crossing on that edge
            // can be the seed segment's start crossing — already emitted as voronoi_inside mesh_start. Closing the loop
            // must not add another mesh vertex there.
            if (!closes_polygon_at_seed_crossing)
            {
              const glm::dvec3 inter_pos = closingMeshVoronoiDelaunayCrossingPosition(
                t, candidate_ref->voronoi_edge_id, candidate_ref->delaunay_edge_id);
              if (std::isfinite(inter_pos.x) && std::isfinite(inter_pos.y))
              {
                walk_append("intersection", cross_base + " (vertex at every Delaunay crossing)",
                  static_cast<size_t>(add_mesh_vertex(inter_pos)));
              }
              else
              {
                KINDS_WARNING("createClosingMesh strand "
                  << strand_id << " t=" << t << " walk " << walk_id
                  << " crossing on Delaunay edge but non-finite position: " << cross_base);
              }
            }

            if (strand_inc)
            {
              if (next_segment_index != static_cast<size_t>(-1))
              {
                if (closes_polygon_at_seed_crossing)
                {
                  found_next_segment = true;
                  closed_walk_at_seed = true;
                  exit_boundary_chain_to_voronoi = true;
                  break;
                }
                // Cyclic incident index: only when handing off to a *different* ordered segment (not loop closure).
                if (next_segment_index != seed_index && num_incident_edges > 0
                  && current_segment->closing_incident_edge_index >= 0)
                {
                  const int expected_next_incident
                    = (current_segment->closing_incident_edge_index + 1) % static_cast<int>(num_incident_edges);
                  const MeshingData* handed = ordered_segments[next_segment_index];
                  if (handed->closing_incident_edge_index < 0
                    || handed->closing_incident_edge_index != expected_next_incident)
                  {
                    KINDS_ERROR("createClosingMesh strand "
                      << strand_id
                      << " closing cap: Delaunay boundary handoff "
                         "crossing matches ordered_seg["
                      << next_segment_index
                      << "] (start_ref) but that segment's "
                         "closing_incident_edge_index="
                      << handed->closing_incident_edge_index
                      << " is not the next strand Voronoi edge after current strip "
                         "(current closing_incident_edge_index="
                      << current_segment->closing_incident_edge_index
                      << ", expected next index=" << expected_next_incident << " in [0," << (num_incident_edges - 1)
                      << "] cyclic around strand). Crossing: " << cross_base);
                  }
                }
                found_next_segment = true;
                if (next_segment_index == seed_index)
                {
                  closed_walk_at_seed = true;
                }
                else
                {
                  current_segment_index = next_segment_index;
                }
                exit_boundary_chain_to_voronoi = true;
                break;
              }
              KINDS_ERROR("createClosingMesh strand "
                << strand_id
                << " closing cap: strand-incident Voronoi crossing on "
                   "Delaunay boundary matches no ordered segment start_ref (and is not seed end_ref loop close): "
                << cross_base);
            }

            current_ref_ptr = cand_ptr;
          }
        }

        if (exit_boundary_chain_to_voronoi)
        {
          break;
        }

        const int corner_vertex_id = graph.destination(boundary_he);
        if (corner_vertex_id >= 0)
        {
          if (static_cast<size_t>(corner_vertex_id) != strand_id)
          {
            KINDS_ERROR("createClosingMesh: strand " << strand_id << " closing cap: Delaunay boundary walk (half-edge "
                                                     << boundary_he << ") reached corner vertex " << corner_vertex_id
                                                     << " which is not the strand vertex (expected " << strand_id
                                                     << ").");
            const glm::dvec2 corner = kin_del.getPointAt(t, static_cast<size_t>(corner_vertex_id));
            const int cv = add_mesh_vertex(glm::dvec3(corner, t));
            walk_append("corner_delaunay",
              std::string("unexpected Delaunay vertex_id=") + std::to_string(corner_vertex_id) + " boundary_he="
                + std::to_string(boundary_he) + " delaunay_edge_id=" + std::to_string(boundary_he / 2),
              static_cast<size_t>(cv));
          }
          else
          {
            append_strand_corner_if_needed();
          }
        }

        boundary_he = kin_del.nextOnComponentBoundaryId(boundary_he);
        current_ref_ptr = nullptr;
      }

      if (!found_next_segment || closed_walk_at_seed)
      {
        break;
      }
    }

    if (polygon.size() >= 3)
    {
      double area2 = 0.0;
      for (size_t i = 0; i < polygon.size(); ++i)
      {
        const glm::dvec3& p0 = mesh.getVertices()[polygon[i]];
        const glm::dvec3& p1 = mesh.getVertices()[polygon[(i + 1) % polygon.size()]];
        area2 += p0.x * p1.y - p1.x * p0.y;
      }
      if (area2 < 0.0)
      {
        std::reverse(polygon.begin() + 1, polygon.end());
      }
      result.polygons.push_back(std::move(polygon));
    }
  }

  return result;
}

void kinDS::SegmentBuilder::closingMeshLogUnmatchedOrderedSegments(
  size_t strand_id, double t, const std::vector<MeshingData*>& ordered_segments, const std::vector<bool>& segment_used)
{
  auto closing_cap_segment_summary = [&](const MeshingData* s) -> std::string
  {
    auto edge_ids = [](const MeshingData* seg) -> std::string
    {
      if (seg->closing_voronoi_edge_id != static_cast<size_t>(-1))
      {
        const size_t ve = seg->closing_voronoi_edge_id;
        return std::string("voronoi_edge_id=") + std::to_string(ve) + " dual_delaunay_edge_id=" + std::to_string(ve);
      }
      return std::string("voronoi_edge_id=? dual_delaunay_edge_id=?");
    };
    return edge_ids(s) + " start_ref=" + kin_del.formatCrossingIntersectionForLog(s->start_intersection_ref)
      + " end_ref=" + kin_del.formatCrossingIntersectionForLog(s->end_intersection_ref);
  };

  size_t unmatched_count = 0;
  for (size_t i = 0; i < segment_used.size(); ++i)
  {
    if (!segment_used[i])
    {
      ++unmatched_count;
      const MeshingData* s = ordered_segments[i];
      KINDS_WARNING("createClosingMesh strand "
        << strand_id << " t=" << t << ": unmatched ordered segment [" << i << "] incident_edge="
        << s->closing_incident_edge_index << " mesh_v(" << s->mesh_start_vertex_id << "->" << s->mesh_end_vertex_id
        << ") he(" << s->start_half_edge_id << "->" << s->end_half_edge_id << ") " << closing_cap_segment_summary(s));
    }
  }
  if (unmatched_count > 0)
  {
    KINDS_WARNING("createClosingMesh strand " << strand_id << " t=" << t << ": " << unmatched_count << " of "
                                              << segment_used.size()
                                              << " ordered Voronoi segment(s) were not reached by polygon tracing.");
  }
}

void kinDS::SegmentBuilder::closingMeshTriangulatePolygonsFan(
  VoronoiMesh& mesh, const std::vector<std::vector<size_t>>& polygons)
{
  for (const auto& polygon : polygons)
  {
    for (size_t i = 2; i < polygon.size(); ++i)
    {
      addMeshletTriangle(mesh, polygon[0], polygon[i - 1], polygon[i]);
    }
  }
}

size_t kinDS::SegmentBuilder::createClosingMesh(
  size_t strand_id, double t, const std::vector<BoundaryPoint>& boundary_polygon, const glm::dvec2& centroid)
{
  const size_t num_incident_edges = closingMeshCountStrandIncidentEdges(strand_id);

  segment_mesh_pairs.push_back(MeshStructure::SegmentMeshPair {});

  VoronoiMesh mesh;
  std::vector<size_t> mesh_vertex_ids;
  mesh_vertex_ids.reserve(32);
  auto track_vertex = [&](const glm::dvec3& pos) -> int
  {
    const int id = closingMeshAppendVertex(mesh, boundary_polygon, centroid, strand_id, t, pos);
    mesh_vertex_ids.push_back(static_cast<size_t>(id));
    return id;
  };

  auto& graph = kin_del.getGraph();
  // Raw segments keyed by undirected Voronoi (dual) edge id; one bucket per incident finite edge.
  std::map<size_t, std::vector<MeshingData>> segments_by_voronoi_edge;
  std::list<MeshingData> closing_segments;

  // Trace every Voronoi edge incident to this strand and keep only portions that lie inside.
  // Use the same half-edge walk as `startNewMesh` (`computeCrossedHalfEdges`) so face inside/outside stays in sync;
  // then attach optional `CrossingData` intersection pointers for Delaunay boundary tracing.
  int incident_edge_index = -1;
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id),
                                                   end = graph.incidentEdgesEnd(strand_id);
    it != end; ++it)
  {
    ++incident_edge_index;
    const size_t incident_he = *it;
    const size_t voronoi_edge_id = incident_he / 2;
    std::vector<MeshingData> bucket = closingMeshExtractRawSegmentsForVoronoiEdge(strand_id, t, mesh, boundary_polygon,
      centroid, incident_edge_index, incident_he, track_vertex);
    segments_by_voronoi_edge[voronoi_edge_id] = bucket;
    closing_segments.insert(closing_segments.end(), std::make_move_iterator(bucket.begin()),
      std::make_move_iterator(bucket.end()));
  }

  ClosingMeshOrderedIndex index_data = closingMeshBuildOrderedIndex(closing_segments);

  closingMeshLogDiscoveredSegments(strand_id, t, closing_segments, index_data.ordered_segments);
  closingMeshValidateOrderedSegmentGeometry(t, mesh, index_data.ordered_segments);

  ClosingMeshPolygonsTraceResult trace
    = closingMeshTraceCapPolygons(strand_id, t, num_incident_edges, mesh, boundary_polygon, centroid,
      std::move(mesh_vertex_ids), index_data.ordered_segments, index_data.start_ref_to_segment);

  closingMeshLogUnmatchedOrderedSegments(strand_id, t, index_data.ordered_segments, trace.segment_used);
  closingMeshTriangulatePolygonsFan(mesh, trace.polygons);

  const size_t index = meshes.size();
  meshes.push_back(std::move(mesh));
  segment_mesh_pair_last_left_and_right_vertex.emplace_back();
  segment_mesh_pair_last_left_and_right_vertex.back() = std::move(closing_segments);

  return index;
}

void kinDS::SegmentBuilder::accumulateSegmentProperties()
{
  // Iterate through all pairs and accumulate properties
  for (size_t pair_id = 0; pair_id < segment_mesh_pairs.size(); ++pair_id)
  {
    auto& pair = segment_mesh_pairs[pair_id];
    if (pair.segment_index0 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index0].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index0].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index0]
          .mesh_pair_indices[segment_properties[pair.segment_index0].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index0].neighbor_indices[segment_properties[pair.segment_index0].neighbor_count]
          = pair.segment_index1; // add neighbor
        segment_properties[pair.segment_index0].neighbor_count++;
      }
    }

    if (pair.segment_index1 != -1)
    {
      // make sure there is space left
      if (segment_properties[pair.segment_index1].neighbor_count >= MeshStructure::SegmentProperties::MAX_NEIGHBORS)
      {
        KINDS_ERROR("Exceeded maximum number of neighbors for segment: "
          << segment_properties[pair.segment_index1].neighbor_count
          << " >= " << MeshStructure::SegmentProperties::MAX_NEIGHBORS);
        // throw std::runtime_error("Exceeded maximum number of neighbors for segment.");
      }
      else
      {
        segment_properties[pair.segment_index1]
          .mesh_pair_indices[segment_properties[pair.segment_index1].neighbor_count] = pair_id; // add mesh pair index
        segment_properties[pair.segment_index1].neighbor_indices[segment_properties[pair.segment_index1].neighbor_count]
          = pair.segment_index0; // add neighbor
        segment_properties[pair.segment_index1].neighbor_count++;
      }
    }
  }
}

void kinDS::SegmentBuilder::splitComponent(
  size_t component_id, const std::vector<std::vector<size_t>>& new_components, double t)
{
  if (new_components.empty())
  {
    return;
  }

  // Update component data
  std::vector<size_t> component_ids(new_components.size(), -1);

  component_ids[0] = component_id;
  size_t new_component_id = kin_del.component_data.components.size();
  for (size_t i = 1; i < new_components.size(); i++)
  {
    component_ids[i] = new_component_id;
    new_component_id++;
  }

  size_t new_size = kin_del.component_data.components.size() + new_components.size() - 1;
  kin_del.component_data.components.resize(new_size);
  kin_del.component_data.component_boundaries.resize(new_size);
  kin_del.component_data.component_centroids.resize(new_size);

  std::vector<bool> he_visited(kin_del.getGraph().getHalfEdges().size(), false);

  for (size_t i = 0; i < new_components.size(); i++)
  {
    size_t cid = component_ids[i];

    for (size_t v : new_components[i])
    {
      kin_del.component_data.component_map[v] = cid;
    }

    kin_del.component_data.components[cid] = new_components[i];
    kin_del.component_data.component_boundaries[cid]
      = kin_del.extractComponentBoundaries(new_components[i], t, he_visited);
    kin_del.component_data.component_centroids[cid]
      = polygonCentroid(kin_del.component_data.component_boundaries[cid][0]);
    kin_del.component_data.component_last_updated[cid] = t;
  }
}

void SegmentBuilder::init()
{
  kin_del.registerEventCallbacks(
    section_callback_.get(), flip_callback_.get(), radius_callback_.get(), crossing_callback_.get());
  kin_del.registerSubdivisionEventCallback(subdivision_callback_.get());

  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.getHalfEdges().size(), -1);
  corner_to_cutoff_mesh_indices.resize(graph.getHalfEdges().size(), -1);

  // Initialize the strand geometries at t = 0.0
  double t = 0.0; // TODO: might be customized later

  // We need a ruled surface for each half-edge in the graph except those having the infinite vertex as
  // origin
  size_t half_edge_count = graph.getHalfEdges().size();

  // initialize segment mesh properties for each strand
  for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
  {
    size_t new_segment_id = segment_properties.size();
    MeshStructure::SegmentProperties properties;
    segment_properties.push_back(properties);
    strand_to_segment_indices[strand_id].push_back(new_segment_id);

    size_t component_index = kin_del.component_data.component_map[strand_id];
    const auto& component = kin_del.component_data.components[component_index];
    const auto& boundary_polygon = kin_del.component_data.component_boundaries[component_index][0];
    const auto& centroid = kin_del.component_data.component_centroids[component_index];
    // create a closing mesh
    size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_polygon, centroid);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[new_segment_id];
    segment_mesh_pair.segment_index0 = -1;
    segment_mesh_pair.segment_index1 = strand_to_segment_indices[strand_id].back();
  }

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    startNewMesh(i, t);
  }

  // initialize boundary mesh
  boundary_mesh_last_left_and_right_vertex.resize(half_edge_count, std::make_pair(-1, -1));
  half_edge_to_boundary_vertex_index.resize(half_edge_count, -1);
  addVoronoiTriangulationToBoundaryMesh(t, false, -0.01);
}

void SegmentBuilder::finalize(double t)
{
  updateBoundaries(t);

  // Finalize the segments by finishing all meshes
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();

  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    auto vertex = graph.getHalfEdges()[i].origin;

    // fall back for infinite vertices
    if (vertex == -1)
    {
      vertex = graph.destination(i);
    }

    size_t component_index = kin_del.component_data.component_map[vertex];
    auto& boundary_points = kin_del.component_data.component_boundaries[component_index][0];

    finishMesh(i, t, boundary_points);
  }

  // finalize closing meshes
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
    // create a closing mesh
    size_t component_index = kin_del.component_data.component_map[strand_id];
    auto& boundary_points = kin_del.component_data.component_boundaries[component_index][0];
    auto& centroid = kin_del.component_data.component_centroids[component_index];
    size_t closing_mesh_index = createClosingMesh(strand_id, t, boundary_points, centroid);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[closing_mesh_index];
    segment_mesh_pair.segment_index0 = strand_to_segment_indices[strand_id].back();
    segment_mesh_pair.segment_index1 = -1;
  }

  accumulateSegmentProperties();

  addVoronoiTriangulationToBoundaryMesh(t, true, 0.01);

  // compute normals
  for (auto& meshlet : meshes)
  {
    meshlet.computeNormals(NormalMode::PerTriangleCorner);
  }

  auto remap1 = boundary_mesh.mergeDuplicateVertices();
  boundary_mesh.removeDegenerateTriangles();
  auto remap2 = boundary_mesh.removeIsolatedVertices();
  boundary_mesh.computeNormals(NormalMode::PerTriangleCorner);

  // Update boundary vertex to strand id mapping
  std::vector<size_t> new_boundary_vertex_to_strand_id;
  new_boundary_vertex_to_strand_id.resize(boundary_mesh.getVertices().size());
  for (size_t old_index = 0; old_index < boundary_vertex_to_strand_id.size(); ++old_index)
  {
    size_t new_index = remap2[remap1[old_index]];
    if (new_index != size_t(-1))
    {
      new_boundary_vertex_to_strand_id[new_index] = boundary_vertex_to_strand_id[old_index];
    }
  }

  boundary_vertex_to_strand_id.swap(new_boundary_vertex_to_strand_id);

  finalized = true; // Set the finalized flag to true
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractMeshes() const { return meshes; }

std::pair<std::vector<VoronoiMesh>, std::vector<std::vector<int>>> kinDS::SegmentBuilder::extractSegmentMeshlets(
  bool merge_by_segment) const
{
  if (!merge_by_segment)
  {
    std::vector<VoronoiMesh> raw_meshlets = meshes;
    std::vector<std::vector<int>> raw_neighbors;
    raw_neighbors.reserve(raw_meshlets.size());

    for (size_t meshlet_index = 0; meshlet_index < raw_meshlets.size(); ++meshlet_index)
    {
      std::vector<int> neighbors_for_meshlet;
      neighbors_for_meshlet.assign(raw_meshlets[meshlet_index].getTriangleCount(), -1);
      raw_neighbors.push_back(std::move(neighbors_for_meshlet));
    }

    return std::make_pair(raw_meshlets, raw_neighbors);
  }

  std::vector<VoronoiMesh> meshlets;
  std::vector<std::vector<int>> neighbor_segments; // accessed as [segment_id][triangle_index]
  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    VoronoiMesh segment_mesh;
    std::vector<int> neighbor_segments_for_meshlet;
    const auto& properties = segment_properties[segment_id];
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      const auto& mesh_pair = segment_mesh_pairs[mesh_pair_index];
      VoronoiMesh mesh = meshes[mesh_pair_index];
      if (segment_mesh_pairs[mesh_pair_index].segment_index0 != segment_id)
      {
        mesh.flipOrientation();
      }
      // Append the mesh to the segment mesh
      neighbor_segments_for_meshlet.insert(
        neighbor_segments_for_meshlet.end(), mesh.getTriangleCount(), properties.neighbor_indices[neighbor_index]);
      segment_mesh += mesh;
    }
    neighbor_segments.push_back(neighbor_segments_for_meshlet);
    segment_mesh.mergeDuplicateVertices(1e-4);
    meshlets.push_back(segment_mesh);
  }

  return std::make_pair(meshlets, neighbor_segments);
}

const VoronoiMesh& kinDS::SegmentBuilder::getBoundaryMesh() const { return boundary_mesh; }

const std::vector<size_t>& kinDS::SegmentBuilder::getBoundaryVertexToStrandId() const
{
  return boundary_vertex_to_strand_id;
}

const std::vector<std::vector<size_t>>& kinDS::SegmentBuilder::getStrandToSegmentIndices() const
{
  return strand_to_segment_indices;
}