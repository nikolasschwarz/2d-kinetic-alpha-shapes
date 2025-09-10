#pragma once
#include "KineticDelaunay.hpp"
#include "Mesh.hpp"
#include "RuledSurface.hpp"

namespace kinDS
{

// forward declaration of MeshBuilder class
class MeshBuilder;

struct StrandGeometry
{
  size_t strand_id; // Unique identifier for the strand
  // std::vector<RuledSurface> ruled_surfaces; // List of ruled surfaces associated with the strand
  std::vector<std::pair<size_t, bool>> ruled_surface_indices; // Indices of the ruled surfaces associated with the strand

  std::vector<Mesh> extractMesh(const MeshBuilder* mesh_builder, std::vector<std::vector<double>>& strand_subdivisions) const;
  void applySubdivision(const std::vector<double>& strand_subdivision, std::vector<RuledSurface>& ruled_surfaces);
};

class MeshBuilder : public KineticDelaunay::EventHandler
{
 private:
  std::vector<StrandGeometry> strand_geometries; // List of strand geometries to be built
  std::vector<int> edge_id_to_ruled_surface; // Maps edge IDs to the corresponding ruled surface index
  const KineticDelaunay& kin_del;
  const std::vector<CubicHermiteSpline<2>>& splines; // Reference to the splines used for the triangulation
  std::vector<RuledSurface> ruled_surfaces; // list of all ruled surfaces
  bool finalized = false; // Flag to indicate if the mesh has been finalized

  VoronoiSiteTrajectory constructTrajectoryForHalfEdge(size_t half_edge_id, size_t section_index = 0) const;

  void insertTrajectoryIntoRuledSurface(size_t half_edge_id, const VoronoiSiteTrajectory& traj, KineticDelaunay::Event& e);

 public:
  MeshBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines);

  void init() override;

  void betweenSections(size_t index) override;

  void beforeEvent(KineticDelaunay::Event& e) override;

  void afterEvent(KineticDelaunay::Event& e) override;

  void finalize() override;

  void applySubdivisions(const std::vector<std::vector<double>>& strand_subdivisions);

  // TODO: this should be const, but it is very hard to implement it as the vector of subdivisions is passed through a lot of functions
  std::vector<Mesh> extractMeshes(double boundary_offset, double subdivision);

  std::vector<Mesh> extractMeshes(double boundary_offset, double subdivision, std::vector<std::vector<double>>& strand_subdivisions);

  const RuledSurface& getRuledSurface(size_t index) const;

  RuledSurface& getRuledSurface(size_t index);

  void printDebugInfo() const;
};
} // namespace kinDS
