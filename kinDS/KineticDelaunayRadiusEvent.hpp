#pragma once

#include "KineticDelaunay.hpp"

#include <stdexcept>

namespace kinDS
{
class KineticDelaunay::RadiusEvent final : public KineticDelaunay::Event
{
public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;

  RadiusEvent(
    KineticDelaunay* kd,
    double t,
    size_t he_id,
    double creation_time,
    glm::dvec2 position)
    : KineticDelaunay::Event(t, creation_time)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
  {
  }

  void handleEvent() override;
};

class KineticDelaunay::RadiusEventManager final : public KineticDelaunay::EventManager
{
public:
  explicit RadiusEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id) override;

private:
  KineticDelaunay* kd_;
};

inline void KineticDelaunay::RadiusEvent::handleEvent()
{
  auto* kd = kd_;
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

  auto* event_handler = kd->radius_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  kd->setFaceInside(face_id, !kd->face_inside[face_id]);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }
}

} // namespace kinDS

