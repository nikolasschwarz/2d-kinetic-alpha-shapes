#pragma once
#include "KineticDelaunay.hpp"
#include "RuledSurface.hpp"

namespace kinDS
{

struct StrandGeometry
{
  size_t strand_id; // Unique identifier for the strand
  std::vector<RuledSurface> ruled_surfaces; // List of ruled surfaces associated with the strand
};

class MeshBuilder : public KineticDelaunay::EventHandler
{
 private:
  std::vector<StrandGeometry> strand_geometries; // List of strand geometries to be built
  std::vector<int> edge_id_to_ruled_surface; // Maps edge IDs to the corresponding ruled surface index
  const KineticDelaunay& kin_del;
  const std::vector<CubicHermiteSpline<2>>& splines; // Reference to the splines used for the triangulation

  VoronoiSiteTrajectory constructTrajectoryForHalfEdge(size_t half_edge_id) const
  {
    // Get the trajectory for the given half-edge ID
    auto& graph = kin_del.getGraph();
    auto triVertices = graph.adjacentTriangleVertices(half_edge_id);

    // TODO: deal with infinite vertex
    // The solution here will probably be to use the two finite vertices and shift outwards by some amount

    Trajectory<2> traj0 = splines[triVertices[0]].getPiecePolynomial(0);
    Trajectory<2> traj1 = splines[triVertices[1]].getPiecePolynomial(0);
    Trajectory<2> traj2 = splines[triVertices[2]].getPiecePolynomial(0);

    return VoronoiSiteTrajectory::circumcenterTrajectory(traj0, traj1, traj2);
  }

  void insertTrajectoryIntoRuledSurface(size_t half_edge_id, const VoronoiSiteTrajectory& traj, double t)
  {
    auto& graph = kin_del.getGraph();
    auto& he = graph.getHalfEdges()[half_edge_id];
    auto& strand = strand_geometries[he.origin];
    strand.ruled_surfaces[edge_id_to_ruled_surface[half_edge_id]].insertRight(traj, t); // Insert the new trajectory into the ruled surface

    // twin
    auto& twin_he = graph.getHalfEdges()[half_edge_id ^ 1];
    auto twin_strand = strand_geometries[twin_he.origin];
    twin_strand.ruled_surfaces[edge_id_to_ruled_surface[half_edge_id ^ 1]].insertLeft(traj, t); // Insert the twin trajectory into the ruled surface
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

        RuledSurface ruled_surface;
        ruled_surface.init(left_traj, right_traj, t); // Initialize the ruled surface with the trajectories and time
        strand_geometries[he.origin].ruled_surfaces.push_back(ruled_surface);
        edge_id_to_ruled_surface[i] = strand_geometries[he.origin].ruled_surfaces.size() - 1; // Map edge ID to ruled surface index
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

    RuledSurface ruled_surface;
    ruled_surface.init(twin_traj, traj, e.time); // Initialize the ruled surface with the trajectories and time
    auto& strandA = strand_geometries[he.origin];
    strandA.ruled_surfaces.push_back(ruled_surface);
    edge_id_to_ruled_surface[e.half_edge_id] = strand_geometries[he.origin].ruled_surfaces.size() - 1; // Map edge ID to ruled surface index

    // now for the twin
    RuledSurface twin_ruled_surface;
    twin_ruled_surface.init(traj, twin_traj, e.time); // Initialize the twin ruled surface with the trajectories and time
    auto& strandB = strand_geometries[twin_he.origin];
    strand_geometries[twin_he.origin].ruled_surfaces.push_back(twin_ruled_surface);
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
  }

  void finalize() override
  {
    // Finalize the mesh by finishing all ruled surfaces
    for (size_t strand_id = 0; strand_id < strand_geometries.size(); ++strand_id)
    {
      auto& strand = strand_geometries[strand_id];

      for (auto& ruled_surface : strand.ruled_surfaces)
      {
        ruled_surface.finalize(splines[strand_id].pointCount() - 1); // Finalize the ruled surface with the upper bound of the spline
      }
    }
  }
};
} // namespace kinDS
