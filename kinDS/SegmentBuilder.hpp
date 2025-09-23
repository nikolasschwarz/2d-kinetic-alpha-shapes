#pragma once
#include "KineticDelaunay.hpp"
#include "Mesh.hpp"
#include "MeshStructure.hpp"
namespace kinDS
{

class SegmentBuilder : public KineticDelaunay::EventHandler
{
 private:
  std::vector<std::vector<size_t>> strand_to_segment_indices; // Maps strand IDs to their corresponding segment indices in correct order
  std::vector<MeshStructure::SegmentMeshProperties> segment_mesh_properties; // Properties for each segment mesh
  std::vector<MeshStructure::SegmentMeshPair> segment_mesh_pairs; // Pairs of segments and their corresponding mesh data
  std::vector<size_t> half_edge_index_to_segment_mesh_pair_index; // Maps edge indices to their corresponding segment mesh pair indices
  std::vector<Mesh> meshes; // List of all generated meshes
  std::vector<std::pair<size_t, size_t>> segment_mesh_pair_last_left_and_right_vertex;

  const KineticDelaunay& kin_del;
  const std::vector<CubicHermiteSpline<2>>& splines; // Reference to the splines used for the triangulation
  bool finalized = false; // Flag to indicate if the mesh has been finalized
  std::vector<std::pair<size_t, double>> subdivisions;
  size_t subdivision_index = 0;

  Point<3> computeVoronoiVertex(size_t half_edge_id, double t, size_t segment_mesh_pair_index) const;

  void finishMesh(size_t half_edge_id, double t);

  void startNewMesh(size_t half_edge_id, double t);

 public:
  SegmentBuilder(const KineticDelaunay& kin_del, std::vector<CubicHermiteSpline<2>>& splines, std::vector<std::pair<size_t, double>> subdivisions);

  void init() override;

  void betweenSections(size_t index) override;

  void beforeEvent(KineticDelaunay::Event& e) override;

  void afterEvent(KineticDelaunay::Event& e) override;

  void insertSubdivision(size_t strand_id, double t);

  void finalize(double t) override;

  std::vector<Mesh> extractMeshes() const;
};
} // namespace kinDS
