#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunaySectionEvent.hpp"

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderSectionCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderSectionCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;
};
} // namespace kinDS

