#include "SegmentBuilder.hpp"

using namespace kinDS;

[[nodiscard]] Point<3> kinDS::SegmentBuilder::computeVoronoiVertex(size_t half_edge_id, double t, size_t segment_mesh_pair_index) const
{
  const auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[half_edge_id];
  const auto& twin_he = graph.getHalfEdges()[half_edge_id ^ 1];

  // Compute the positions of the Voronoi vertices at t = 0.0
  // First get the two adjacent triangles
  std::array<int, 3> triVertices = graph.adjacentTriangleVertices(half_edge_id);

  // now compute the circumcenters if the triangles are not infinite
  Point<2> circumcenter;

  bool infinite = false;

  std::vector<Point<2>> points;

  if (triVertices[0] != -1)
  {
    points.push_back(splines[triVertices[0]].evaluate(t));
  }

  if (triVertices[1] != -1)
  {
    points.push_back(splines[triVertices[1]].evaluate(t));
  }

  if (triVertices[2] != -1)
  {
    points.push_back(splines[triVertices[2]].evaluate(t));
  }

  if (points.size() == 3)
  {
    circumcenter = graph.circumcenter(points[0], points[1], points[2]);
  }
  else
  {
    infinite = true;

    // For now just take the midpoint of the edge
    circumcenter = (points[0] + points[1]) * 0.5;
  }

  // place circumcenters into the mesh
  return Point<3> { circumcenter[0], circumcenter[1], t };
}

void kinDS::SegmentBuilder::finishMesh(size_t he_id, double t)
{
  size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
  // Get corresponding mesh
  Mesh& mesh = meshes[segment_mesh_pair_index];
  // Insert Voronoi vertex
  size_t new_left_vertex_index = mesh.getVertices().size();
  Point<3> left_vertex = computeVoronoiVertex(he_id & ~1, t, segment_mesh_pair_index);
  mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);
  size_t new_right_vertex_index = mesh.getVertices().size();
  Point<3> right_vertex = computeVoronoiVertex((he_id & ~1) + 1, t, segment_mesh_pair_index);
  mesh.addVertex(right_vertex[0], right_vertex[1], right_vertex[2]);
  // build triangles
  const auto& last_vertices = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  // create two triangles
  mesh.addTriangle(last_vertices.first, last_vertices.second, new_left_vertex_index);
  mesh.addTriangle(new_left_vertex_index, last_vertices.second, new_right_vertex_index);
  // update last vertex indices
  segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index]
    = std::make_pair(new_left_vertex_index, new_right_vertex_index);
}

SegmentBuilder::SegmentBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines, std::vector<std::pair<size_t, double>> subdivisions)
  : kin_del(kin_del)
  , splines(splines)
  , subdivisions(std::move(subdivisions))
{
  // Assert that the subdivisions are sorted by time
  assert(std::is_sorted(this->subdivisions.begin(), this->subdivisions.end(), [](const auto& a, const auto& b)
    { return a.second <= b.second; }));

  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.getHalfEdges().size(), -1);

  // size_t half_edge_count = graph.getHalfEdges().size();
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

  // For now also create a mesh, but this might be changed later
  Mesh mesh;

  Point<3> left_vertex = computeVoronoiVertex(even_id, t, half_edge_index_to_segment_mesh_pair_index[even_id]);
  Point<3> right_vertex = computeVoronoiVertex(odd_id, t, half_edge_index_to_segment_mesh_pair_index[even_id]);

  mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);
  mesh.addVertex(right_vertex[0], right_vertex[1], right_vertex[2]);
  meshes.push_back(mesh);

  // add last vertex indices
  segment_mesh_pair_last_left_and_right_vertex.push_back(
    std::make_pair(mesh.getVertices().size() - 2, mesh.getVertices().size() - 1));

  assert(segment_mesh_pairs.size() == segment_mesh_pair_last_left_and_right_vertex.size());
}

void SegmentBuilder::init()
{
  // Initialize the strand geometries at t = 0.0
  double t = 0.0; // TODO: might be customized later

  // We need a ruled surface for each half-edge in the graph with the exeption of those having the infinite vertex as origin
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();

  // initialize segment mesh properties for each strand
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
    size_t segment_mesh_properties_id = segment_mesh_properties.size();
    MeshStructure::SegmentMeshProperties properties;
    segment_mesh_properties.push_back(properties);
    strand_to_segment_indices[strand_id].push_back(segment_mesh_properties_id);
  }

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    startNewMesh(i, t);
  }
}

void SegmentBuilder::betweenSections(size_t index)
{
  // Check if we need to insert a subdivision before handling this event
  while (subdivision_index < subdivisions.size() && subdivisions[subdivision_index].second <= index)
  {
    insertSubdivision(subdivisions[subdivision_index].first, subdivisions[subdivision_index].second);
    subdivision_index++;
  }

  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    const auto& he = graph.getHalfEdges()[i];
    size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[i];

    // Get corresponding mesh
    Mesh& mesh = meshes[segment_mesh_pair_index];

    // Insert Voronoi vertex
    Point<3> left_vertex = computeVoronoiVertex(i, index, segment_mesh_pair_index);
    size_t new_left_vertex_index = mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);

    Point<3> right_vertex = computeVoronoiVertex(i ^ 1, index, segment_mesh_pair_index);
    size_t new_right_vertex_index = mesh.addVertex(right_vertex[0], right_vertex[1], right_vertex[2]);

    // build triangles
    const auto& last_vertices = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    // create two triangles
    if (last_vertices.first != last_vertices.second) // case of first insertion
    {
      mesh.addTriangle(last_vertices.first, last_vertices.second, new_left_vertex_index);
    }
    mesh.addTriangle(new_left_vertex_index, last_vertices.second, new_right_vertex_index);
    // update last vertex indices
    segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index]
      = std::make_pair(new_left_vertex_index, new_right_vertex_index);
  }
}

void SegmentBuilder::beforeEvent(KineticDelaunay::Event& e)
{
  // Check if we need to insert a subdivision before handling this event
  while (subdivision_index < subdivisions.size() && subdivisions[subdivision_index].second <= e.time)
  {
    insertSubdivision(subdivisions[subdivision_index].first, subdivisions[subdivision_index].second);
    subdivision_index++;
  }

  // Finish the segment mesh pair of the edge being flipped
  Point<3> event_point { e.position[0], e.position[1], e.time };
  size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[e.half_edge_id];
  Mesh& mesh = meshes[segment_mesh_pair_index];
  size_t event_vertex_index = mesh.addVertex(event_point[0], event_point[1], event_point[2]);
  const auto& last_vertices = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
  // create one triangle to the event point
  mesh.addTriangle(last_vertices.first, last_vertices.second, event_vertex_index);
}

void SegmentBuilder::afterEvent(KineticDelaunay::Event& e)
{
  auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[e.half_edge_id];
  const auto& twin_he = graph.getHalfEdges()[e.half_edge_id ^ 1];
  // Create a new segment mesh pair for the two new edges created by the flip
  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pair.segment_index0 = he.origin == -1 ? -1 : strand_to_segment_indices[he.origin].back();
  segment_mesh_pair.segment_index1 = twin_he.origin == -1 ? -1 : strand_to_segment_indices[twin_he.origin].back();

  half_edge_index_to_segment_mesh_pair_index[e.half_edge_id] = segment_mesh_pairs.size();
  half_edge_index_to_segment_mesh_pair_index[e.half_edge_id ^ 1] = segment_mesh_pairs.size();

  segment_mesh_pairs.push_back(segment_mesh_pair);

  // For now also create a mesh, but this might be changed later
  Mesh mesh;
  size_t index = mesh.addVertex(e.position[0], e.position[1], e.time);

  // add last vertex indices
  segment_mesh_pair_last_left_and_right_vertex.push_back(std::make_pair(index, index));

  meshes.push_back(mesh);

  // first get the other half-edges of the quadrilateral
  size_t he0_id = graph.getHalfEdges()[e.half_edge_id].next; // Next half-edge in the quadrilateral
  size_t he1_id = graph.getHalfEdges()[he0_id].next; // Next half-edge in the quadrilateral
  size_t he2_id = graph.getHalfEdges()[e.half_edge_id ^ 1].next; // Next half-edge in the quadrilatera
  size_t he3_id = graph.getHalfEdges()[he2_id].next; // Next half-edge in the quadrilateral

  // for each of them, insert the event vertex with one triangle on one side
  for (size_t he_id : { he0_id, he1_id, he2_id, he3_id })
  {
    size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
    Mesh& mesh = meshes[segment_mesh_pair_index];
    size_t index = mesh.addVertex(e.position[0], e.position[1], e.time);
    const auto& last_vertices = segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index];
    // create one triangle to the event point
    mesh.addTriangle(last_vertices.first, last_vertices.second, index);
    // update last vertex indices
    const MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[segment_mesh_pair_index];
    // Determine whether we have to update the left or right vertex here
    if (segment_mesh_pair.segment_index1 == strand_to_segment_indices[graph.getHalfEdges()[he_id].origin].back())
    {
      segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index] = std::make_pair(last_vertices.first, index);
    }
    else
    {
      segment_mesh_pair_last_left_and_right_vertex[segment_mesh_pair_index] = std::make_pair(index, last_vertices.second);
    }
  }
}

void kinDS::SegmentBuilder::insertSubdivision(size_t strand_id, double t)
{
  logger.log(DEBUG, "Inserting subdivision for strand %zu at t = %f", strand_id, t);
  // Traverse all half-edges around this strand and insert a new vertex into the corresponding segment meshes
  auto& graph = kin_del.getGraph();

  // finish old meshes
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id), end = graph.incidentEdgesEnd(strand_id); it != end; ++it)
  {
    finishMesh(*it, t);
  }

  // Create a new segment mesh property for the new segment
  size_t segment_mesh_properties_id = segment_mesh_properties.size();
  MeshStructure::SegmentMeshProperties properties;
  segment_mesh_properties.push_back(properties);
  strand_to_segment_indices[strand_id].push_back(segment_mesh_properties_id);

  // Start new meshes
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id), end = graph.incidentEdgesEnd(strand_id); it != end; ++it)
  {
    startNewMesh(*it, t);

    // insert vertices into adjacent meshes
    auto& he = graph.getHalfEdges()[*it];

    size_t adjacent_he_id = he.next;
    auto& adjacent_he = graph.getHalfEdges()[adjacent_he_id];
    size_t adjacent_segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[adjacent_he_id];
    auto& adjacent_segment_mesh_pair = segment_mesh_pairs[adjacent_segment_mesh_pair_index];
    Mesh& adjacent_mesh = meshes[adjacent_segment_mesh_pair_index];
    Point<3> vertex = computeVoronoiVertex(adjacent_he_id, t, adjacent_segment_mesh_pair_index);
    size_t new_vertex_index = adjacent_mesh.addVertex(vertex[0], vertex[1], vertex[2]);
    auto& last_vertices = segment_mesh_pair_last_left_and_right_vertex[adjacent_segment_mesh_pair_index];
    adjacent_mesh.addTriangle(last_vertices.first, last_vertices.second, new_vertex_index);

    if (adjacent_he_id % 2 == 0)
    {
      last_vertices.first = new_vertex_index;
    }
    else
    {
      last_vertices.second = new_vertex_index;
    }
  }
}

void SegmentBuilder::finalize(double t)
{
  // Check if we need to insert a subdivision before handling this event
  while (subdivision_index < subdivisions.size() && subdivisions[subdivision_index].second <= t)
  {
    insertSubdivision(subdivisions[subdivision_index].first, subdivisions[subdivision_index].second);
    subdivision_index++;
  }

  // Finalize the mesh by finishing all ruled surfaces
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();

  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    finishMesh(i, t);
  }

  finalized = true; // Set the finalized flag to true
}

std::vector<Mesh> kinDS::SegmentBuilder::extractMeshes() const
{
  return meshes;
}
