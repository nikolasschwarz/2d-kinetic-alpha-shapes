#pragma once

#include "KineticDelaunay.hpp"
#include "Logger.hpp"

#include <stdexcept>

namespace kinDS
{
class KineticDelaunay::RadiusEvent final : public KineticDelaunay::Event
{
public:
  KineticDelaunay* kd_;
  size_t half_edge_id;
  glm::dvec2 position;
  /// Target state implied by the radius predicate sign change. Positive -> negative means the triangle is added
  /// (inside); negative -> positive means it is removed (outside). This makes duplicate events idempotent.
  bool target_inside;

  RadiusEvent(
    KineticDelaunay* kd,
    double t,
    size_t he_id,
    double creation_time,
    glm::dvec2 position,
    bool target_inside)
    : KineticDelaunay::Event(t, creation_time, 30u)
    , kd_(kd)
    , half_edge_id(he_id)
    , position(position)
    , target_inside(target_inside)
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

  // Outdated if the event half-edge or its triangle was tombstoned (e.g. after a branch split).
  if (!graph.isLiveHalfEdge(half_edge_id))
  {
    return;
  }

  size_t face_id = graph.halfEdge(half_edge_id).face;
  if (!graph.isLiveFace(face_id))
  {
    return;
  }

  // Check if the event is still valid
  if (creation_time < kd->face_last_updated[face_id])
  {
    return;
  }

  if (kd->mustRemainInside(face_id, occurrence_time))
  {
    if (!kd->face_inside[face_id])
    {
      kd->setFaceInside(face_id, true, occurrence_time);
    }
    return;
  }

  if (kd->face_inside[face_id] == target_inside)
  {
    KINDS_DEBUG("RadiusEvent no-op: face " << face_id << " is already "
                                           << (target_inside ? "inside" : "outside")
                                           << " at t=" << occurrence_time
                                           << " for half_edge_id=" << half_edge_id);
    return;
  }

  auto* event_handler = kd->radius_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  kd->setFaceInside(face_id, target_inside, occurrence_time);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  kd->validateSitesInsideConvexHull("RadiusEvent:afterEvent", occurrence_time);
}

} // namespace kinDS

