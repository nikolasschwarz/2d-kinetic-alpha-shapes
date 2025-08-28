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
  std::vector<size_t> ruled_surface_indices; // Indices of the ruled surfaces associated with the strand

  Mesh extractMesh(const MeshBuilder* mesh_builder) const;
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

  void insertTrajectoryIntoRuledSurface(size_t half_edge_id, const VoronoiSiteTrajectory& traj, double t);

 public:
  MeshBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines);

  void init() override;

  void betweenSections(size_t index) override;

  void beforeEvent(KineticDelaunay::Event& e) override;

  void afterEvent(KineticDelaunay::Event& e) override;

  void finalize() override;

  std::vector<Mesh> extractMeshes(double boundary_offset, double subdivision) const;

  const RuledSurface& getRuledSurface(size_t index) const;

  RuledSurface& getRuledSurface(size_t index);
};
} // namespace kinDS
