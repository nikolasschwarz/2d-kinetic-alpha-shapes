#pragma once

#include "KineticDelaunay.hpp"

namespace kinDS
{
class KineticDelaunay::SectionEvent final : public KineticDelaunay::Event
{
public:
    SectionEvent(
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

class KineticDelaunay::SectionEventManager final : public KineticDelaunay::EventManager
{
public:
    explicit SectionEventManager(KineticDelaunay* kd)
    : EventManager(kd)
    {
    }

    void computeEvents(double t, size_t event_id) override;
};

inline void KineticDelaunay::SectionEvent::handleEvent(EventHandler& event_handler)
{
    auto* kd = getKineticDelaunay();
    if (!kd)
    {
        throw std::runtime_error("SectionEvent has no KineticDelaunay pointer");
    }
}
}