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
  KineticEventType eventType() const override { return KineticEventType::Radius; }
};

class KineticDelaunay::RadiusEventManager final : public KineticDelaunay::EventManager
{
public:
  explicit RadiusEventManager(KineticDelaunay* kd)
    : kd_(kd)
  {
  }

  void computeEvents(double t, size_t event_id,
    std::optional<InfinitesimalComputeContext> infinitesimal = std::nullopt) override;

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
  const double t = occurrence_time.real_time;
  const double infinitesimal_t = occurrence_time.infinitesimal_time;
  const bool is_infinitesimal = infinitesimal_t > 0.0;

  if (is_infinitesimal)
  {
    bool epoch_ok = false;
    for (const auto& entry : kd->pending_branch_splits_.by_parent_)
    {
      if (entry.second.infinitesimal_active && entry.second.infinitesimal_epoch == infinitesimal_epoch_)
      {
        epoch_ok = true;
        break;
      }
    }
    if (!epoch_ok)
    {
      return;
    }
    kd->current_infinitesimal_t_ = infinitesimal_t;
  }

  const auto clear_infinitesimal = [&]()
  {
    if (is_infinitesimal)
    {
      kd->current_infinitesimal_t_ = 0.0;
    }
  };

  // Outdated if the event half-edge or its triangle was tombstoned (e.g. after a branch split).
  if (!graph.isLiveHalfEdge(half_edge_id))
  {
    clear_infinitesimal();
    return;
  }

  size_t face_id = graph.halfEdge(half_edge_id).face;
  if (!graph.isLiveFace(face_id))
  {
    clear_infinitesimal();
    return;
  }

  // Check if the event is still valid
  if (creation_time < kd->face_last_updated[face_id])
  {
    clear_infinitesimal();
    return;
  }

  if (kd->mustRemainInside(face_id, t))
  {
    if (!kd->face_inside[face_id])
    {
      kd->setFaceInside(face_id, true, t);
    }
    clear_infinitesimal();
    return;
  }

  if (kd->face_inside[face_id] == target_inside)
  {
    KINDS_DEBUG("RadiusEvent no-op: face " << face_id << " is already "
                                           << (target_inside ? "inside" : "outside")
                                           << " at t=" << t
                                           << " for half_edge_id=" << half_edge_id);
    clear_infinitesimal();
    return;
  }

  auto* event_handler = kd->radius_event_manager_->getCallback();
  if (event_handler)
  {
    event_handler->beforeEvent(*this);
  }

  kd->setFaceInside(face_id, target_inside, t);

  if (event_handler)
  {
    event_handler->afterEvent(*this);
  }

  kd->validateSitesInsideConvexHull("RadiusEvent:afterEvent", t);

  if (is_infinitesimal)
  {
    size_t parent_component_id = static_cast<size_t>(-1);
    for (const auto& entry : kd->pending_branch_splits_.by_parent_)
    {
      if (entry.second.infinitesimal_active && entry.second.infinitesimal_epoch == infinitesimal_epoch_)
      {
        parent_component_id = entry.first;
        break;
      }
    }
    // Finalize first. Kinetic radius events do not recompute neighbors afterward.
    if (parent_component_id != static_cast<size_t>(-1))
    {
      kd->maybeFinalizeInfinitesimalSeparation(parent_component_id, t);
    }
    kd->current_infinitesimal_t_ = 0.0;
  }
}

} // namespace kinDS

