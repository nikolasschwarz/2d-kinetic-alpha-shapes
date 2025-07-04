#pragma once
#include "KineticDelaunay.hpp"

namespace MyNamespace
{

class StrandGeometry
{
    size_t strand_id; // Unique identifier for the strand
    std::vector<RuledSurface> ruled_surfaces; // List of ruled surfaces associated with the strand
};

class MeshBuilder
{
public:
    void operator()(KineticDelaunay::Event& e, bool before)
    {
        // Default event handler does nothing
        // This can be overridden by the user to handle events as needed
    }
};
} // namespace MyNamespace
