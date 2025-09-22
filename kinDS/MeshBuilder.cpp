#include "MeshBuilder.hpp"

using namespace kinDS;

std::vector<Mesh> StrandGeometry::extractMesh(const MeshBuilder* mesh_builder, std::vector<std::vector<double>>& strand_subdivisions) const
{
  /* std::vector<Point<3>> vertices;
  std::vector<size_t> indices;
  std::vector<size_t> group_offsets;

  for (auto& [index, inverted] : ruled_surface_indices)
  {
    group_offsets.emplace_back(indices.size() / 3); // Store the current vertex count as a group offset
    mesh_builder->getRuledSurface(index).extractTriangles(vertices, indices, inverted); // Extract triangles from the ruled surface
  }
  Mesh result(vertices, indices); // Return the constructed mesh

  result.setGroupOffsets(group_offsets); // Set the group offsets for the mesh
  return result;*/

  // Re-do this such that we obtain meshlets according to the strand subdivisions
  std::vector<Mesh> meshes;
  for (auto& [index, inverted] : ruled_surface_indices)
  {
    std::vector<Mesh> extracted = mesh_builder->getRuledSurface(index).extractTriangles(inverted, strand_subdivisions, strand_id);
    std::copy(extracted.begin(), extracted.end(), std::back_inserter(meshes));
  }

  return meshes;
}

void kinDS::StrandGeometry::applySubdivision(const std::vector<double>& strand_subdivision, std::vector<RuledSurface>& ruled_surfaces)
{
  // iterate through the ruled surfaces and apply the subdivision
  for (auto& [index, inverted] : ruled_surface_indices)
  {
    RuledSurface& rs = ruled_surfaces[index];
    rs.applySubdivision(strand_subdivision);
  }
}

VoronoiSiteTrajectory MeshBuilder::constructTrajectoryForHalfEdge(size_t half_edge_id, size_t section_index) const
{
  // Get the trajectory for the given half-edge ID
  auto& graph = kin_del.getGraph();
  auto triVertices = graph.adjacentTriangleVertices(half_edge_id);

  int infinite_vertex = -1;
  std::vector<Trajectory<2>> trajs;
  for (int i = 0; i < 3; i++)
  {
    if (triVertices[i] == -1)
    {
      infinite_vertex = i;
    }
    else
    {
      trajs.push_back(splines[triVertices[i]].getPiecePolynomial(section_index)); // Get the trajectory for the vertex
    }
  }

  // No infinite vertex found, construct the circumcenter trajectory
  if (infinite_vertex == -1)
  {
    return VoronoiSiteTrajectory::circumcenterTrajectory(trajs[0], trajs[1], trajs[2]);
  }

  // Else, just let them run along the center of the edge, any offsets will be computed at mesh extraction
  return VoronoiSiteTrajectory { (trajs[0][0] + trajs[1][0]) / 2, (trajs[0][1] + trajs[1][1]) / 2 }; // Return the trajectory for the edge opposite to the infinite vertex
}

void MeshBuilder::insertTrajectoryIntoRuledSurface(size_t half_edge_id, const VoronoiSiteTrajectory& traj, KineticDelaunay::Event& e)
{
  auto& graph = kin_del.getGraph();
  auto& he = graph.getHalfEdges()[half_edge_id];
  auto& strand = strand_geometries[he.origin];
  RuledSurface& ruled_surface = ruled_surfaces[edge_id_to_ruled_surface[half_edge_id]];

  if (half_edge_id % 2 == 0)
  {
    ruled_surface.insertRightEvent(e.position, traj, e.time); // Insert the new trajectory into the ruled surface
    logger.log(DEBUG, "Inserted into right side of strand %i, surface no. %i", he.origin, edge_id_to_ruled_surface[half_edge_id]);
  }
  else
  {
    ruled_surface.insertLeftEvent(e.position, traj, e.time); // Insert the new trajectory into the ruled surface
    logger.log(DEBUG, "Inserted into left side of strand %i, surface no. %i", he.origin, edge_id_to_ruled_surface[half_edge_id]);
  }

  // twin
  /* auto& twin_he = graph.getHalfEdges()[half_edge_id ^ 1];
  auto& twin_strand = strand_geometries[twin_he.origin];
  RuledSurface& twin_ruled_surface = ruled_surfaces[edge_id_to_ruled_surface[half_edge_id ^ 1]];
  twin_ruled_surface.insertLeft(traj, t); // Insert the twin trajectory into the ruled surface
  logger.log(DEBUG, "Inserted into left side of strand %i, surface no. %i", twin_he.origin, edge_id_to_ruled_surface[half_edge_id ^ 1]);

  logger.log(DEBUG, "Inserted trajectory into ruled surface for half-edge %zu at time %f between strands %i and %i", half_edge_id, t, he.origin, twin_he.origin);*/
}

MeshBuilder::MeshBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines)
  : kin_del(kin_del)
  , splines(splines)
{
  auto& graph = kin_del.getGraph();

  size_t strand_count = graph.getVertexCount();
  strand_geometries.reserve(strand_count);

  for (size_t i = 0; i < strand_count; ++i)
  {
    strand_geometries.emplace_back(StrandGeometry { i, {} });
  }

  size_t half_edge_count = graph.getHalfEdges().size();
  edge_id_to_ruled_surface.resize(half_edge_count, -1); // Initialize with -1 to indicate no ruled surface
}

void MeshBuilder::init()
{
  // Initialize the strand geometries at t = 0.0
  double t = 0.0; // TODO: might be customized later

  // We need a ruled surface for each half-edge in the graph with the exeption of those having the infinite vertex as origin
  auto& graph = kin_del.getGraph();

  size_t half_edge_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    const auto& he = graph.getHalfEdges()[i];

    // Build trajectories for each Voronoi vertex
    // TODO: avoid constructing trajectories for the same half-edge twice
    VoronoiSiteTrajectory right_traj = constructTrajectoryForHalfEdge(i);
    VoronoiSiteTrajectory left_traj = constructTrajectoryForHalfEdge(i ^ 1); // Get the twin half-edge trajectory

    RuledSurface ruled_surface(graph.getHalfEdges()[i].face, he.face);
    ruled_surface.init(left_traj, right_traj, t); // Initialize the ruled surface with the trajectories and time
    ruled_surfaces.push_back(ruled_surface);

    edge_id_to_ruled_surface[i] = ruled_surfaces.size() - 1; // Map edge ID to ruled surface index
    if (he.origin != -1) // Skip half-edges with the infinite vertex as origin
    {
      strand_geometries[he.origin].ruled_surface_indices.push_back(std::make_pair(ruled_surfaces.size() - 1, false));
    }

    edge_id_to_ruled_surface[i ^ 1] = ruled_surfaces.size() - 1; // Map edge ID to twin ruled surface index
    if (graph.getHalfEdges()[i ^ 1].origin != -1) // Skip half-edges with the infinite vertex as origin
    {
      strand_geometries[graph.getHalfEdges()[i ^ 1].origin].ruled_surface_indices.push_back(std::make_pair(ruled_surfaces.size() - 1, true));
    }
  }
}

void MeshBuilder::betweenSections(size_t index)
{
  auto& graph = kin_del.getGraph();
  size_t half_edge_count = graph.getHalfEdges().size();
  for (size_t i = 0; i < half_edge_count; i += 2)
  {
    const auto& he = graph.getHalfEdges()[i];

    // Build trajectories for each Voronoi vertex
    // TODO: avoid constructing trajectories for the same half-edge twice
    VoronoiSiteTrajectory right_traj = constructTrajectoryForHalfEdge(i, index);
    VoronoiSiteTrajectory left_traj = constructTrajectoryForHalfEdge(i ^ 1, index); // Get the twin half-edge trajectory

    RuledSurface& ruled_surface = ruled_surfaces[edge_id_to_ruled_surface[i]];

    ruled_surface.insertLeft(left_traj, index); // Insert the new section into the ruled surface
    ruled_surface.insertRight(right_traj, index); // Insert the new section into the twin ruled surface
  }
}

void MeshBuilder::beforeEvent(KineticDelaunay::Event& e)
{
  // Finish the ruled surface in the point resulting from the edge flip and insert it into the adjacent strand geometry
  auto& graph = kin_del.getGraph();
  auto& he = graph.getHalfEdges()[e.half_edge_id];
  RuledSurface& s0 = ruled_surfaces[edge_id_to_ruled_surface[e.half_edge_id]];
  s0.finalizeEvent(e.position, e.time); // Finalize the ruled surface at the event time
}

void MeshBuilder::afterEvent(KineticDelaunay::Event& e)
{
  auto& graph = kin_del.getGraph();
  const auto& he = graph.getHalfEdges()[e.half_edge_id];
  const auto& twin_he = graph.getHalfEdges()[e.half_edge_id ^ 1];
  // Start a new ruled surface in the point resulting from the edge flip
  VoronoiSiteTrajectory traj = constructTrajectoryForHalfEdge(e.half_edge_id);
  VoronoiSiteTrajectory twin_traj = constructTrajectoryForHalfEdge(e.half_edge_id ^ 1); // Get the twin half-edge trajectory

  RuledSurface ruled_surface(twin_he.face, he.face);
  ruled_surface.initEvent(e.position, twin_traj, traj, e.time); // Initialize the ruled surface with the trajectories and time
  auto& strandA = strand_geometries[he.origin];
  ruled_surfaces.push_back(ruled_surface);
  edge_id_to_ruled_surface[e.half_edge_id] = ruled_surfaces.size() - 1; // Map edge ID to ruled surface index
  strand_geometries[he.origin].ruled_surface_indices.push_back(std::make_pair(ruled_surfaces.size() - 1, false));

  // re-use the same ruled surface, just invert the sides
  edge_id_to_ruled_surface[e.half_edge_id ^ 1] = ruled_surfaces.size() - 1; // Map edge ID to twin ruled surface index
  strand_geometries[twin_he.origin].ruled_surface_indices.push_back(std::make_pair(ruled_surfaces.size() - 1, true));

  // first get the other half-edges of the quadrilateral
  size_t he0_id = graph.getHalfEdges()[e.half_edge_id].next; // Next half-edge in the quadrilateral
  size_t he1_id = graph.getHalfEdges()[he0_id].next; // Next half-edge in the quadrilateral
  size_t he2_id = graph.getHalfEdges()[e.half_edge_id ^ 1].next; // Next half-edge in the quadrilatera
  size_t he3_id = graph.getHalfEdges()[he2_id].next; // Next half-edge in the quadrilateral

  // now update the ruled surfaces for the other half-edges
  insertTrajectoryIntoRuledSurface(he0_id, traj, e); // Insert the new trajectory into the ruled surface for the first half-edge
  insertTrajectoryIntoRuledSurface(he1_id, traj, e); // Insert the new trajectory into the ruled surface for the second half-edge
  insertTrajectoryIntoRuledSurface(he2_id, twin_traj, e); // Insert the twin trajectory into the ruled surface for the third half-edge
  insertTrajectoryIntoRuledSurface(he3_id, twin_traj, e); // Insert the twin trajectory into the ruled surface for the fourth half-edge

  // For debugging, iterate through all ruled surfaces and print their debug info
  logger.log(DEBUG, "Debugging ruled surfaces");
  for (const auto& ruled_surface : ruled_surfaces)
  {
    ruled_surface.printDebugInfo(); // Print debug info for each ruled surface
    // visual separator in log
    logger.log(DEBUG, "----------------------------------------");
  }
}

void MeshBuilder::finalize(double t)
{
  // Finalize the mesh by finishing all ruled surfaces
  for (size_t strand_id = 0; strand_id < strand_geometries.size(); ++strand_id)
  {
    auto& strand = strand_geometries[strand_id];

    for (auto& [ruled_surface_index, inverted] : strand.ruled_surface_indices)
    {
      if (!getRuledSurface(ruled_surface_index).isFinalized())
        getRuledSurface(ruled_surface_index).finalize(splines[strand_id].pointCount() - 1); // Finalize the ruled surface with the upper bound of the spline
    }
  }

  finalized = true; // Set the finalized flag to true
}

void kinDS::MeshBuilder::applySubdivisions(const std::vector<std::vector<double>>& strand_subdivisions)
{
  // first of all, we insert the strand subdivisions to all ruled surfaces
  for (size_t strand_id = 0; strand_id < strand_geometries.size(); ++strand_id)
  {
    StrandGeometry& strand = strand_geometries[strand_id];
    strand.applySubdivision(strand_subdivisions[strand_id], ruled_surfaces);
  }
}

std::vector<Mesh> MeshBuilder::extractMeshes(double boundary_offset, double subdivision)
{
  std::vector<std::vector<double>> strand_subdivisions(strand_geometries.size());

  return extractMeshes(boundary_offset, subdivision, strand_subdivisions);
}

std::vector<Mesh> MeshBuilder::extractMeshes(double boundary_offset, double subdivision, std::vector<std::vector<double>>& strand_subdivisions)
{
  assert(strand_subdivisions.size() == strand_geometries.size() && "Strand subdivisions size must match the number of strands.");
  logger.log(INFO, "Extracting meshes with boundary offset %f and subdivision %f", boundary_offset, subdivision);

  if (!finalized)
  {
    throw std::runtime_error("MeshBuilder must be finalized before extracting the mesh.");
  }

  applySubdivisions(strand_subdivisions); // Apply the strand subdivisions to the ruled surfaces

  std::vector<Mesh> meshes;
  for (size_t strand_id = 0; strand_id < strand_geometries.size(); ++strand_id)
  {
    assert(strand_geometries[strand_id].strand_id == strand_id && "Strand ID mismatch.");

    std::vector<Mesh> extracted = strand_geometries[strand_id].extractMesh(this, strand_subdivisions); // Add the mesh for the strand to the list of meshes
    std::copy(extracted.begin(), extracted.end(), std::back_inserter(meshes)); // TODO: we do this twice, this could be really inefficient if the compiler does not optimize it away
  }

  return meshes; // Return the list of extracted meshes
}

const RuledSurface& MeshBuilder::getRuledSurface(size_t index) const
{
  if (index >= ruled_surfaces.size())
  {
    throw std::out_of_range("Index out of range for ruled surfaces.");
  }
  return ruled_surfaces[index]; // Return the ruled surface at the specified index
}

RuledSurface& MeshBuilder::getRuledSurface(size_t index)
{
  if (index >= ruled_surfaces.size())
  {
    throw std::out_of_range("Index out of range for ruled surfaces.");
  }
  return ruled_surfaces[index]; // Return the ruled surface at the specified index
}

void MeshBuilder::printDebugInfo() const
{
  for (size_t i = 0; i < strand_geometries.size(); ++i)
  {
    const auto& strand = strand_geometries[i];
    logger.log(INFO, "======== Strand %zu has %zu ruled surfaces ========", i, strand.ruled_surface_indices.size());
    for (const auto& [ruled_surface_index, inverted] : strand.ruled_surface_indices)
    {
      logger.log(INFO, "  Ruled Surface Index: %zu, Inverted: %s", ruled_surface_index, inverted ? "true" : "false");
      ruled_surfaces[ruled_surface_index].printDebugInfo("    ");
    }
  }
}
