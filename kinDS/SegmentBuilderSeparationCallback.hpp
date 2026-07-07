#pragma once

#include "KineticDelaunay.hpp"
#include "KineticDelaunaySeparationEvent.hpp"

namespace kinDS
{
class SegmentBuilder;

class SegmentBuilderSeparationCallback final : public KineticDelaunay::EventCallback
{
 public:
  explicit SegmentBuilderSeparationCallback(SegmentBuilder& segment_builder)
    : segment_builder_(segment_builder)
  {
  }

  void beforeEvent(KineticDelaunay::Event& e) override;
  void afterEvent(KineticDelaunay::Event& e) override;

 private:
  SegmentBuilder& segment_builder_;
  void writeSeparationVisualDebugSvg(const KineticDelaunay::SeparationEvent& separation, const char* phase) const;
};
} // namespace kinDS
