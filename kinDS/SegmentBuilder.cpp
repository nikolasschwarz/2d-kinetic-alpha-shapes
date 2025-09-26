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
  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
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
  // split quad differently depending on which side is closer
  if (last_vertices.first == last_vertices.second)
  {
    mesh.addTriangle(new_left_vertex_index, last_vertices.second, new_right_vertex_index);
  }
  else if (mesh.getVertices()[last_vertices.first][2] < mesh.getVertices()[last_vertices.second][2])
  {
    mesh.addTriangle(last_vertices.first, last_vertices.second, new_left_vertex_index);
    mesh.addTriangle(new_left_vertex_index, last_vertices.second, new_right_vertex_index);
  }
  else
  {
    mesh.addTriangle(last_vertices.first, last_vertices.second, new_right_vertex_index);
    mesh.addTriangle(last_vertices.first, new_right_vertex_index, new_left_vertex_index);
  }

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
}

SegmentBuilder::SegmentBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines)
  : kin_del(kin_del)
  , splines(splines)
{
  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.getHalfEdges().size(), -1);
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
  VoronoiMesh mesh;

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

size_t kinDS::SegmentBuilder::createClosingMesh(size_t strand_id, double t)
{
  auto& graph = kin_del.getGraph();

  MeshStructure::SegmentMeshPair segment_mesh_pair;
  segment_mesh_pairs.push_back(segment_mesh_pair);

  VoronoiMesh mesh;

  // we just create a triangle fan because the Voronoi cell is convex
  // iterate over all segment indices of the strand
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id), end = graph.incidentEdgesEnd(strand_id); it != end; ++it)
  {
    Point<3> voronoi_vertex = computeVoronoiVertex(*it, t, half_edge_index_to_segment_mesh_pair_index[*it]);
    mesh.addVertex(voronoi_vertex[0], voronoi_vertex[1], voronoi_vertex[2]);
  }

  // create triangles
  size_t center_index = 0;

  for (size_t voronoi_index = 2; voronoi_index < mesh.getVertices().size(); ++voronoi_index)
  {
    mesh.addTriangle(center_index, voronoi_index - 1, voronoi_index);
  }

  size_t index = meshes.size();
  meshes.push_back(mesh);
  segment_mesh_pair_last_left_and_right_vertex.push_back(
    std::make_pair(-1, -1)); // not needed, so we just set it to -1,-1

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
      assert(segment_properties[pair.segment_index0].neighbor_count < MeshStructure::SegmentProperties::MAX_NEIGHBORS);

      segment_properties[pair.segment_index0].mesh_pair_indices[segment_properties[pair.segment_index0].neighbor_count]
        = pair_id; // add mesh pair index
      segment_properties[pair.segment_index0].neighbor_indices[segment_properties[pair.segment_index0].neighbor_count]
        = pair.segment_index1; // add neighbor
      segment_properties[pair.segment_index0].neighbor_count++;
    }

    if (pair.segment_index1 != -1)
    {
      // make sure there is space left
      assert(segment_properties[pair.segment_index1].neighbor_count < MeshStructure::SegmentProperties::MAX_NEIGHBORS);

      segment_properties[pair.segment_index1].mesh_pair_indices[segment_properties[pair.segment_index1].neighbor_count]
        = pair_id; // add mesh pair index
      segment_properties[pair.segment_index1].neighbor_indices[segment_properties[pair.segment_index1].neighbor_count]
        = pair.segment_index0; // add neighbor
      segment_properties[pair.segment_index1].neighbor_count++;
    }
  }
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
    size_t new_segment_id = segment_properties.size();
    MeshStructure::SegmentProperties properties;
    segment_properties.push_back(properties);
    strand_to_segment_indices[strand_id].push_back(new_segment_id);

    // create a closing mesh
    size_t closing_mesh_index = createClosingMesh(strand_id, t);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[new_segment_id];
    segment_mesh_pair.segment_index0 = -1;
    segment_mesh_pair.segment_index1 = strand_to_segment_indices[strand_id].back();
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
    finishMesh(i, index);
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
  VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
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
  VoronoiMesh mesh;
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
    VoronoiMesh& mesh = meshes[segment_mesh_pair_index];
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

  size_t new_segment_id = segment_properties.size();

  // create a closing mesh
  size_t closing_mesh_index = createClosingMesh(strand_id, t);
  MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[closing_mesh_index];
  segment_mesh_pair.segment_index0 = strand_to_segment_indices[strand_id].back();
  segment_mesh_pair.segment_index1 = new_segment_id;

  // Create a new segment mesh property for the new segment
  MeshStructure::SegmentProperties properties;
  segment_properties.push_back(properties);
  strand_to_segment_indices[strand_id].push_back(new_segment_id);

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
    VoronoiMesh& adjacent_mesh = meshes[adjacent_segment_mesh_pair_index];
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

  // Finalize the segments by finishing all meshes
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();

  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    finishMesh(i, t);
  }

  // finalize closing meshes
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
    // create a closing mesh
    size_t closing_mesh_index = createClosingMesh(strand_id, t);
    MeshStructure::SegmentMeshPair& segment_mesh_pair = segment_mesh_pairs[closing_mesh_index];
    segment_mesh_pair.segment_index0 = strand_to_segment_indices[strand_id].back();
    segment_mesh_pair.segment_index1 = -1;
  }

  accumulateSegmentProperties();

  finalized = true; // Set the finalized flag to true
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractMeshes() const
{
  return meshes;
}

std::vector<VoronoiMesh> kinDS::SegmentBuilder::extractSegmentMeshlets() const
{
  std::vector<VoronoiMesh> meshlets;

  for (size_t segment_id = 0; segment_id < segment_properties.size(); ++segment_id)
  {
    VoronoiMesh segment_mesh;
    const auto& properties = segment_properties[segment_id];
    for (size_t neighbor_index = 0; neighbor_index < properties.neighbor_count; ++neighbor_index)
    {
      size_t mesh_pair_index = properties.mesh_pair_indices[neighbor_index];
      const auto& mesh_pair = segment_mesh_pairs[mesh_pair_index];
      const auto& mesh = meshes[mesh_pair_index];
      // Append the mesh to the segment mesh
      segment_mesh += mesh;
    }
    meshlets.push_back(segment_mesh);
  }

  return meshlets;
}
