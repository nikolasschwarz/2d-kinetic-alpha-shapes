#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunayRadiusEvent.hpp"

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderRadiusCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderRadiusCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;
};
} // namespace kinDS

