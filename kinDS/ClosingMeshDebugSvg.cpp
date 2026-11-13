#include "ClosingMeshDebugSvg.hpp"

#include <sstream>

namespace kinDS
{
namespace
{
std::vector<std::string> wrapText(const std::string& text, size_t width)
{
  std::vector<std::string> lines;
  std::istringstream words(text);
  std::string word, line;
  while (words >> word)
  {
    if (!line.empty() && line.size() + 1 + word.size() > width)
    {
      lines.push_back(line);
      line.clear();
    }
    line += (line.empty() ? "" : " ") + word;
  }
  if (!line.empty())
  {
    lines.push_back(line);
  }
  return lines;
}

std::array<int, 3> stripColor(ClosingMeshWalkDebug::StripState state)
{
  switch (state)
  {
  case ClosingMeshWalkDebug::StripState::Unused:
    return { 214, 39, 40 };
  case ClosingMeshWalkDebug::StripState::Walked:
    return { 31, 119, 180 };
  default:
    return { 150, 150, 150 };
  }
}

const char* stripStateName(ClosingMeshWalkDebug::StripState state)
{
  switch (state)
  {
  case ClosingMeshWalkDebug::StripState::Unused:
    return "UNUSED";
  case ClosingMeshWalkDebug::StripState::Walked:
    return "walked";
  default:
    return "closed";
  }
}
} // namespace

HalfEdgeDelaunayGraphToSVG::DebugOverlay buildClosingMeshWalkOverlay(const ClosingMeshWalkDebug& debug)
{
  using Anchor = ClosingMeshWalkDebug::Anchor;
  HalfEdgeDelaunayGraphToSVG::DebugOverlay overlay;

  for (size_t k = 0; k < debug.boundary_edges.size(); ++k)
  {
    const auto& edge = debug.boundary_edges[k];
    overlay.paths.push_back(
      { { Anchor { Anchor::Kind::Site, edge.origin }, Anchor { Anchor::Kind::Site, edge.destination } },
        { 255, 127, 14 }, 0.05, "b" + std::to_string(k) + " he" + std::to_string(edge.half_edge), false });
  }
  for (const auto& strip : debug.strips)
  {
    const bool unused = strip.state == ClosingMeshWalkDebug::StripState::Unused;
    overlay.paths.push_back({ { strip.start, strip.end }, stripColor(strip.state), unused ? 0.04 : 0.03,
      "s" + std::to_string(strip.index), false });
  }
  std::vector<Anchor> walk;
  for (const auto& point : debug.walk_points)
  {
    walk.push_back(point.anchor);
  }
  overlay.paths.push_back({ walk, { 44, 160, 44 }, 0.012, "", true });
  if (debug.stuck.has_value())
  {
    overlay.markers.push_back({ debug.stuck->anchor, "STUCK: " + debug.stuck->label, { 220, 0, 160 } });
  }

  auto& notes = overlay.notes;
  notes.push_back("Closing mesh cap walk FAILED: strand " + std::to_string(debug.strand_id)
    + " t=" + std::to_string(debug.t) + " seed s" + std::to_string(debug.seed_index));
  for (const auto& line : wrapText(debug.reason, 90))
  {
    notes.push_back("  " + line);
  }
  std::string order = "walk order:";
  for (size_t index : debug.walk_strip_order)
  {
    order += " s" + std::to_string(index);
  }
  notes.push_back(order);
  notes.push_back("strips (red UNUSED, blue walked, gray closed earlier), start -> end:");
  for (const auto& strip : debug.strips)
  {
    notes.push_back("  s" + std::to_string(strip.index) + " [" + stripStateName(strip.state) + "] " + strip.label + ": "
      + strip.start_label + " -> " + strip.end_label);
  }
  std::string steps = "walk steps (green):";
  for (size_t i = 0; i < debug.walk_points.size(); ++i)
  {
    steps += " " + std::to_string(i) + "=" + debug.walk_points[i].label;
  }
  for (const auto& line : wrapText(steps, 90))
  {
    notes.push_back(line);
  }
  return overlay;
}
} // namespace kinDS
