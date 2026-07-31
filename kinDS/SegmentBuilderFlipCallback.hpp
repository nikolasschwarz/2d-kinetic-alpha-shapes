#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayFlipEvent.hpp"

#include <glm/glm.hpp>

#include <optional>

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderFlipCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderFlipCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;
  /// Canonical (smaller) flip-edge Voronoi vertex id captured in @ref beforeEvent when finite coincidence applies.
  std::optional<size_t> buffered_flip_voronoi_vertex_id_;
  /// Mesh-space and Delaunay XY for that coincidence point, computed once from pre-flip topology/strands.
  std::optional<glm::dvec3> buffered_flip_mesh_position_;
  std::optional<glm::dvec2> buffered_flip_delaunay_xy_;
};
} // namespace kinDS
