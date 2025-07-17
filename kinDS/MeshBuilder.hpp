#pragma once
#include "KineticDelaunay.hpp"
#include "Mesh.hpp"
#include "RuledSurface.hpp"

namespace kinDS
{

struct StrandGeometry
{
  size_t strand_id; // Unique identifier for the strand
  std::vector<RuledSurface> ruled_surfaces; // List of ruled surfaces associated with the strand

  Mesh extractMesh() const
  {
    std::vector<Point<3>> vertices;
    std::vector<size_t> indices;
    std::vector<size_t> group_offsets;

    for (const auto& ruled_surface : ruled_surfaces)
    {
      group_offsets.emplace_back(indices.size() / 3); // Store the current vertex count as a group offset
      ruled_surface.extractTriangles(vertices, indices); // Extract triangles from the ruled surface
    }
    Mesh result(vertices, indices); // Return the constructed mesh

    result.setGroupOffsets(group_offsets); // Set the group offsets for the mesh
    return result;
  }
};

class MeshBuilder : public KineticDelaunay::EventHandler
{
 private:
  std::vector<StrandGeometry> strand_geometries; // List of strand geometries to be built
  std::vector<int> edge_id_to_ruled_surface; // Maps edge IDs to the corresponding ruled surface index
  const KineticDelaunay& kin_del;
  const std::vector<CubicHermiteSpline<2>>& splines; // Reference to the splines used for the triangulation
  bool finalized = false; // Flag to indicate if the mesh has been finalized

  VoronoiSiteTrajectory constructTrajectoryForHalfEdge(size_t half_edge_id, size_t section_index = 0) const
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

  void insertTrajectoryIntoRuledSurface(size_t half_edge_id, const VoronoiSiteTrajectory& traj, double t)
  {
    auto& graph = kin_del.getGraph();
    auto& he = graph.getHalfEdges()[half_edge_id];
    auto& strand = strand_geometries[he.origin];
    RuledSurface& ruled_surface = strand.ruled_surfaces[edge_id_to_ruled_surface[half_edge_id]];
    ruled_surface.insertRight(traj, t); // Insert the new trajectory into the ruled surface
    logger.log(DEBUG, "Inserted into right side of strand %i, surface no. %i", he.origin, edge_id_to_ruled_surface[half_edge_id]);

    // twin
    auto& twin_he = graph.getHalfEdges()[half_edge_id ^ 1];
    auto& twin_strand = strand_geometries[twin_he.origin];
    RuledSurface& twin_ruled_surface = twin_strand.ruled_surfaces[edge_id_to_ruled_surface[half_edge_id ^ 1]];
    twin_ruled_surface.insertLeft(traj, t); // Insert the twin trajectory into the ruled surface
    logger.log(DEBUG, "Inserted into left side of strand %i, surface no. %i", twin_he.origin, edge_id_to_ruled_surface[half_edge_id ^ 1]);

    logger.log(DEBUG, "Inserted trajectory into ruled surface for half-edge %zu at time %f between strands %i and %i", half_edge_id, t, he.origin, twin_he.origin);
  }

 public:
  MeshBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines)
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

  void init() override
  {
    // Initialize the strand geometries at t = 0.0
    double t = 0.0; // TODO: might be customized later

    // We need a ruled surface for each half-edge in the graph with the exeption of those having the infinite vertex as origin
    auto& graph = kin_del.getGraph();

    size_t half_edge_count = graph.getHalfEdges().size();
    for (size_t i = 0; i < half_edge_count; ++i)
    {
      const auto& he = graph.getHalfEdges()[i];
      if (he.origin != -1) // Skip half-edges with the infinite vertex as origin
      {
        // Build trajectories for each Voronoi vertex
        // TODO: avoid constructing trajectories for the same half-edge twice
        VoronoiSiteTrajectory right_traj = constructTrajectoryForHalfEdge(i);
        VoronoiSiteTrajectory left_traj = constructTrajectoryForHalfEdge(i ^ 1); // Get the twin half-edge trajectory

        RuledSurface ruled_surface(graph.getHalfEdges()[i].face, he.face);
        ruled_surface.init(left_traj, right_traj, t); // Initialize the ruled surface with the trajectories and time
        strand_geometries[he.origin].ruled_surfaces.push_back(ruled_surface);
        edge_id_to_ruled_surface[i] = strand_geometries[he.origin].ruled_surfaces.size() - 1; // Map edge ID to ruled surface index
      }
    }
  }

  void betweenSections(size_t index) override
  {
    auto& graph = kin_del.getGraph();
    size_t half_edge_count = graph.getHalfEdges().size();
    for (size_t i = 0; i < half_edge_count; ++i)
    {
      const auto& he = graph.getHalfEdges()[i];
      if (he.origin != -1) // Skip half-edges with the infinite vertex as origin
      {
        // Build trajectories for each Voronoi vertex
        // TODO: avoid constructing trajectories for the same half-edge twice
        VoronoiSiteTrajectory right_traj = constructTrajectoryForHalfEdge(i, index);
        VoronoiSiteTrajectory left_traj = constructTrajectoryForHalfEdge(i ^ 1, index); // Get the twin half-edge trajectory

        RuledSurface& ruled_surface = strand_geometries[he.origin].ruled_surfaces[edge_id_to_ruled_surface[i]];

        ruled_surface.insertLeft(left_traj, index); // Insert the new section into the ruled surface
        ruled_surface.insertRight(right_traj, index); // Insert the new section into the twin ruled surface
      }
    }
  }

  void beforeEvent(KineticDelaunay::Event& e) override
  {
    // Finish the ruled surface in the point resulting from the edge flip and insert it into the adjacent strand geometry
    auto& graph = kin_del.getGraph();
    auto& he = graph.getHalfEdges()[e.half_edge_id];
    RuledSurface& s0 = strand_geometries[he.origin].ruled_surfaces[edge_id_to_ruled_surface[e.half_edge_id]];
    s0.finalize(e.time); // Finalize the ruled surface at the event time

    // same for twin half-edge
    auto& twin_he = graph.getHalfEdges()[e.half_edge_id ^ 1];
    RuledSurface& s1 = strand_geometries[twin_he.origin].ruled_surfaces[edge_id_to_ruled_surface[e.half_edge_id ^ 1]];
    s1.finalize(e.time); // Finalize the twin ruled surface at the event time
  }

  void afterEvent(KineticDelaunay::Event& e) override
  {
    auto& graph = kin_del.getGraph();
    const auto& he = graph.getHalfEdges()[e.half_edge_id];
    const auto& twin_he = graph.getHalfEdges()[e.half_edge_id ^ 1];
    // Start a new ruled surface in the point resulting from the edge flip
    VoronoiSiteTrajectory traj = constructTrajectoryForHalfEdge(e.half_edge_id);
    VoronoiSiteTrajectory twin_traj = constructTrajectoryForHalfEdge(e.half_edge_id ^ 1); // Get the twin half-edge trajectory

    RuledSurface ruled_surface(twin_he.face, he.face);
    ruled_surface.init(twin_traj, traj, e.time); // Initialize the ruled surface with the trajectories and time
    auto& strandA = strand_geometries[he.origin];
    strandA.ruled_surfaces.push_back(ruled_surface);
    edge_id_to_ruled_surface[e.half_edge_id] = strand_geometries[he.origin].ruled_surfaces.size() - 1; // Map edge ID to ruled surface index

    // now for the twin
    RuledSurface twin_ruled_surface(he.face, twin_he.face);
    twin_ruled_surface.init(traj, twin_traj, e.time); // Initialize the twin ruled surface with the trajectories and time
    auto& strandB = strand_geometries[twin_he.origin];
    strandB.ruled_surfaces.push_back(twin_ruled_surface);
    edge_id_to_ruled_surface[e.half_edge_id ^ 1] = strand_geometries[twin_he.origin].ruled_surfaces.size() - 1; // Map edge ID to twin ruled surface index

    // first get the other half-edges of the quadrilateral
    size_t he0_id = graph.getHalfEdges()[e.half_edge_id].next; // Next half-edge in the quadrilateral
    size_t he1_id = graph.getHalfEdges()[he0_id].next; // Next half-edge in the quadrilateral
    size_t he2_id = graph.getHalfEdges()[e.half_edge_id ^ 1].next; // Next half-edge in the quadrilatera
    size_t he3_id = graph.getHalfEdges()[he2_id].next; // Next half-edge in the quadrilateral

    // now update the ruled surfaces for the other half-edges
    insertTrajectoryIntoRuledSurface(he0_id, traj, e.time); // Insert the new trajectory into the ruled surface for the first half-edge
    insertTrajectoryIntoRuledSurface(he1_id, traj, e.time); // Insert the new trajectory into the ruled surface for the second half-edge
    insertTrajectoryIntoRuledSurface(he2_id, twin_traj, e.time); // Insert the twin trajectory into the ruled surface for the third half-edge
    insertTrajectoryIntoRuledSurface(he3_id, twin_traj, e.time); // Insert the twin trajectory into the ruled surface for the fourth half-edge

    // For debugging, iterate through all ruled surfaces and print their debug info
    for (size_t strand_id = 0; strand_id < strand_geometries.size(); strand_id++)
    {
      logger.log(DEBUG, "Debugging ruled surfaces for strand %zu", strand_id);
      for (const auto& ruled_surface : strand_geometries[strand_id].ruled_surfaces)
      {
        ruled_surface.printDebugInfo(); // Print debug info for each ruled surface
        // visual separator in log
        logger.log(DEBUG, "----------------------------------------");
      }
    }
  }

  void finalize() override
  {
    // Finalize the mesh by finishing all ruled surfaces
    for (size_t strand_id = 0; strand_id < strand_geometries.size(); ++strand_id)
    {
      auto& strand = strand_geometries[strand_id];

      for (auto& ruled_surface : strand.ruled_surfaces)
      {
        if (!ruled_surface.isFinalized())
          ruled_surface.finalize(splines[strand_id].pointCount() - 1); // Finalize the ruled surface with the upper bound of the spline
      }
    }

    finalized = true; // Set the finalized flag to true
  }

  std::vector<Mesh> extractMeshes(double boundary_offset, double subdivision) const
  {
    if (!finalized)
    {
      throw std::runtime_error("MeshBuilder must be finalized before extracting the mesh.");
    }

    std::vector<Mesh> meshes;
    for (const auto& strand : strand_geometries)
    {

      meshes.push_back(strand.extractMesh()); // Add the mesh for the strand to the list of meshes
    }

    return meshes; // Return the list of extracted meshes
  }
};
} // namespace kinDS
