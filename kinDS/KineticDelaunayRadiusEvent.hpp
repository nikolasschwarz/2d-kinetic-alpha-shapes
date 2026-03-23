#pragma once

#include "KineticDelaunay.hpp"

#include <stdexcept>

namespace kinDS
{
class KineticDelaunay::RadiusEvent final : public KineticDelaunay::Event
{
public:
  RadiusEvent(
    KineticDelaunay* kd,
    double t,
    size_t he_id,
    double creation_time,
    glm::dvec2 position)
    : KineticDelaunay::Event(kd, t, he_id, creation_time, position, static_cast<size_t>(-1))
  {
  }

  void handleEvent(EventHandler& event_handler) override;
};

class KineticDelaunay::RadiusEventManager final : public KineticDelaunay::EventManager
{
public:
  explicit RadiusEventManager(KineticDelaunay* kd)
    : EventManager(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;
};

inline void KineticDelaunay::RadiusEvent::handleEvent(EventHandler& event_handler)
{
  auto* kd = getKineticDelaunay();
  if (!kd)
  {
    throw std::runtime_error("RadiusEvent has no KineticDelaunay pointer");
  }

  auto& graph = kd->graph;

  // Check if the event is still valid
  size_t face_id = graph.getHalfEdges()[half_edge_id].face;
  if (creation_time < kd->face_last_updated[face_id])
  {
    return;
  }

  event_handler.beforeRadiusEvent(*this);

  kd->setFaceInside(face_id, !kd->face_inside[face_id]);

  event_handler.afterRadiusEvent(*this);
}

} // namespace kinDS

