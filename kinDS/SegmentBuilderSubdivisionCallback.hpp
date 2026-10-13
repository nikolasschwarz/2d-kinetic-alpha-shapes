#pragma once

#include "KineticDelaunay.hpp"

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderSubdivisionCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderSubdivisionCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;
};
} // namespace kinDS
