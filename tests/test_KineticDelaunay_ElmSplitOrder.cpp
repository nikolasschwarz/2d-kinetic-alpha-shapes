#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/KineticDelaunayCrossingEvent.hpp"
#include "kinDS/KineticDelaunayFlipEvent.hpp"
#include "kinDS/StrandTree.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>

using namespace kinDS;

namespace
{
std::filesystem::path elmFixturePath()
{
  const char* relative = "tests/fixtures/elm_branch_0_sections_14_15.txt";
  for (const std::filesystem::path prefix : { ".", "..", "../.." })
  {
    if (std::filesystem::exists(prefix / relative))
    {
      return prefix / relative;
    }
  }
  throw std::runtime_error(std::string("test fixture not found: ") + relative);
}

/// Mirrors SegmentBuilder's flip callback: restore the order of the flipped quad's crossing lists.
struct FlipResort : public KineticDelaunay::EventCallback
{
  KineticDelaunay& kd;
  explicit FlipResort(KineticDelaunay& kd_ref)
    : kd(kd_ref)
  {
  }

  void afterEvent(KineticDelaunay::Event& event) override
  {
    const auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&event);
    if (flip == nullptr)
    {
      return;
    }
    const HalfEdgeDelaunayGraph& graph = kd.getGraph();
    std::unordered_set<size_t> edges { flip->half_edge_id / 2 };
    for (size_t he_id : graph.getQuadBoundaryHalfEdgeIndices(flip->half_edge_id / 2))
    {
      edges.insert(he_id / 2);
    }
    for (size_t edge_id : edges)
    {
      if (!graph.isInfinite(2 * edge_id))
      {
        kd.refreshAndSortDelaunayEdgeIntersectionParams(edge_id, flip->occurrence_time);
      }
    }
  }
};

/// After each crossing event, the crossed Delaunay edge's crossings must be ordered by position along the edge.
struct CrossedEdgeOrderCheck : public KineticDelaunay::EventCallback
{
  KineticDelaunay& kd;
  std::optional<std::string> first_violation;
  size_t checked_events = 0;
  explicit CrossedEdgeOrderCheck(KineticDelaunay& kd_ref)
    : kd(kd_ref)
  {
  }

  void afterEvent(KineticDelaunay::Event& event) override
  {
    const auto* crossing = dynamic_cast<KineticDelaunay::CrossingEvent*>(&event);
    if (crossing == nullptr || first_violation.has_value())
    {
      return;
    }
    ++checked_events;
    const size_t edge_id = crossing->half_edge_id / 2;
    // Just after the event, when the new crossings have separated along the edge.
    const double t = crossing->occurrence_time + 1e-6;
    double previous = -std::numeric_limits<double>::infinity();
    std::ostringstream list;
    bool sorted = true;
    for (const auto& ref : kd.getCrossingData().delaunay_edge_intersections[edge_id])
    {
      const double d_param = kd.delaunayVoronoiEdgeIntersection(edge_id, ref->voronoi_edge_id, t).first;
      list << " ve" << ref->voronoi_edge_id << '@' << d_param;
      sorted = sorted && d_param >= previous;
      previous = d_param;
    }
    if (!sorted)
    {
      first_violation = "de=" + std::to_string(edge_id) + " t=" + std::to_string(crossing->occurrence_time) + list.str();
    }
  }
};
} // namespace

// Elm tree (TreeBatchApp seed 0 job 3), input branch 0 at global sections 14-15. Two splits on one Delaunay edge
// at t=14.1689 used to be inserted in the wrong order (the Voronoi edges' far endpoints lie at very different
// distances from the edge), which later broke merge_2_to_1 with "could not match middle strip".
TEST_CASE("Elm slice keeps split crossings ordered along the Delaunay edge", "[KineticDelaunay][CrossingData]")
{
  StrandTree tree = StrandTree::loadFromFile(elmFixturePath());
  REQUIRE(tree.getPoints().size() == 544);

  constexpr double cutoff = 10.0;
  KineticDelaunay kd(tree, cutoff, false);
  FlipResort flip_resort(kd);
  CrossedEdgeOrderCheck order_check(kd);
  kd.registerEventCallbacks(nullptr, &flip_resort, nullptr, &order_check);
  kd.init();
  REQUIRE_NOTHROW(kd.compute());

  INFO(order_check.first_violation.value_or("no violation"));
  REQUIRE(order_check.checked_events > 0);
  REQUIRE_FALSE(order_check.first_violation.has_value());
}
