#include "VoronoiMesh.hpp"

namespace kinDS
{
// Provides a structure to hold the mesh data for the segments as it will later be on the GPU
struct MeshStructure
{
  struct SegmentMeshPair
  {
    size_t segment_index0;
    size_t segment_index1;
    size_t mesh_start_index; // start index into the index buffer
    size_t mesh_length; // length of the mesh in the index buffer
    size_t connected = 1;
  };

  struct IntersectionMeshPairMetadata
  {
    size_t voronoi_cell_id = static_cast<size_t>(-1);
    size_t owner_segment_id = static_cast<size_t>(-1);
    size_t start_delaunay_edge_id = static_cast<size_t>(-1);
    size_t end_delaunay_edge_id = static_cast<size_t>(-1);
  };

  struct SegmentProperties
  {
    std::vector<size_t> mesh_pair_indices; // every segment mesh pair this segment's meshlet is assembled from
  };

 private:
  VoronoiMesh mesh;

  std::vector<SegmentMeshPair> segment_pairs; // pairs of (start, end) indices into the index buffer for each meshlet
};
} // namespace kinDS