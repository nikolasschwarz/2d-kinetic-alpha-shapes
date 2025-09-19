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
    circumcenter = (splines[he.origin].evaluate(t) + splines[twin_he.origin].evaluate(t)) * 0.5;
  }

  // place circumcenters into the mesh
  return Point<3> { circumcenter[0], circumcenter[1], t };
}

SegmentBuilder::SegmentBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines)
  : kin_del(kin_del)
  , splines(splines)
{
  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_to_segment_indices.resize(strand_count);
  half_edge_index_to_segment_mesh_pair_index.resize(graph.getHalfEdges().size(), -1);

  // size_t half_edge_count = graph.getHalfEdges().size();
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
    strand_to_segment_indices[strand_id].push_back(segment_mesh_properties_id);
  }

  // now go through all half-edges and create a segment mesh pair
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    const auto& graph = kin_del.getGraph();
    const auto& he = graph.getHalfEdges()[i];
    const auto& twin_he = graph.getHalfEdges()[i ^ 1];

    MeshStructure::SegmentMeshPair segment_mesh_pair;
    segment_mesh_pair.segment_index0 = strand_to_segment_indices[he.origin].back();
    segment_mesh_pair.segment_index1 = strand_to_segment_indices[twin_he.origin].back();

    half_edge_index_to_segment_mesh_pair_index[i] = segment_mesh_pairs.size();
    half_edge_index_to_segment_mesh_pair_index[i ^ 1] = segment_mesh_pairs.size();

    segment_mesh_pairs.push_back(segment_mesh_pair);

    // For now also create a mesh, but this might be changed later
    Mesh mesh;

    Point<3> left_vertex = computeVoronoiVertex(i, t, half_edge_index_to_segment_mesh_pair_index[i]);
    Point<3> right_vertex = computeVoronoiVertex(i ^ 1, t, half_edge_index_to_segment_mesh_pair_index[i]);

    mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);
    mesh.addVertex(right_vertex[0], right_vertex[1], right_vertex[2]);
    meshes.push_back(mesh);

    // add last vertex indices
    segment_mesh_pair_last_left_and_right_vertex.push_back(
      std::make_pair(mesh.getVertices().size() - 2, mesh.getVertices().size() - 1));

    assert(segment_mesh_pairs.size() == segment_mesh_pair_last_left_and_right_vertex.size());
  }
}

void SegmentBuilder::betweenSections(size_t index)
{
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    const auto& he = graph.getHalfEdges()[i];
    size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[i];

    // Get corresponding mesh
    Mesh& mesh = meshes[segment_mesh_pair_index];

    // Insert Voronoi vertex
    size_t new_left_vertex_index = mesh.getVertices().size();
    Point<3> left_vertex = computeVoronoiVertex(i, index, segment_mesh_pair_index);
    mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);

    size_t new_right_vertex_index = mesh.getVertices().size();
    Point<3> right_vertex = computeVoronoiVertex(i ^ 1, index, segment_mesh_pair_index);
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
}

void SegmentBuilder::beforeEvent(KineticDelaunay::Event& e)
{
  // Add event point to all adjacent segment meshes
}

void SegmentBuilder::afterEvent(KineticDelaunay::Event& e)
{
  auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[e.half_edge_id];
  const auto& twin_he = graph.getHalfEdges()[e.half_edge_id ^ 1];
  // Start a new ruled surface in the point resulting from the edge flip

  // re-use the same ruled surface, just invert the sides

  // first get the other half-edges of the quadrilateral
  size_t he0_id = graph.getHalfEdges()[e.half_edge_id].next; // Next half-edge in the quadrilateral
  size_t he1_id = graph.getHalfEdges()[he0_id].next; // Next half-edge in the quadrilateral
  size_t he2_id = graph.getHalfEdges()[e.half_edge_id ^ 1].next; // Next half-edge in the quadrilatera
  size_t he3_id = graph.getHalfEdges()[he2_id].next; // Next half-edge in the quadrilateral
}

void kinDS::SegmentBuilder::insertSubdivision(size_t strand_id, double t)
{
  // Traverse all half-edges around this strand and insert a new vertex into the corresponding segment meshes
  auto& graph = kin_del.getGraph();
  for (HalfEdgeDelaunayGraph::IncidentEdgeIterator it = graph.incidentEdgesBegin(strand_id), end = graph.incidentEdgesEnd(strand_id); it != end; ++it)
  {
    size_t he_id = *it;
    size_t segment_mesh_pair_index = half_edge_index_to_segment_mesh_pair_index[he_id];
    // Get corresponding mesh
    Mesh& mesh = meshes[segment_mesh_pair_index];
    // Insert Voronoi vertex
    size_t new_left_vertex_index = mesh.getVertices().size();
    Point<3> left_vertex = computeVoronoiVertex(he_id, t, segment_mesh_pair_index);
    mesh.addVertex(left_vertex[0], left_vertex[1], left_vertex[2]);
    size_t new_right_vertex_index = mesh.getVertices().size();
    Point<3> right_vertex = computeVoronoiVertex(he_id ^ 1, t, segment_mesh_pair_index);
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
}

void SegmentBuilder::finalize()
{
  // Finalize the mesh by finishing all ruled surfaces
  auto& graph = kin_del.getGraph();
  for (size_t strand_id = 0; strand_id < graph.getVertexCount(); ++strand_id)
  {
  }

  finalized = true; // Set the finalized flag to true
}
