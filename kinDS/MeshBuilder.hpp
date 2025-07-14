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
 public:
  MeshBuilder(std::vector<size_t> strand_ids)
  {
    strand_geometries.reserve(strand_ids.size());
    for (const auto& id : strand_ids)
    {
      strand_geometries.emplace_back(StrandGeometry { id, {} });
    }
  }

  MeshBuilder(size_t strand_count)
  {
    strand_geometries.reserve(strand_count);
    for (size_t i = 0; i < strand_count; ++i)
    {
      strand_geometries.emplace_back(StrandGeometry { i, {} });
    }
  }

  void init()
  {
    // Initialize the strand geometries at t = 0.0
  }

  void beforeEvent(KineticDelaunay::Event& e)
  {
    // Finish the ruled surface in the point resulting from the edge flip and insert it into the adjacent strand geometry
  }

  void afterEvent(KineticDelaunay::Event& e) override
  {
    // Start a new ruled surface in the point resulting from the edge flip
  }

  void finalize()
  {
    // Finalize the mesh by finishing all ruled surfaces
  }
};
} // namespace kinDS
