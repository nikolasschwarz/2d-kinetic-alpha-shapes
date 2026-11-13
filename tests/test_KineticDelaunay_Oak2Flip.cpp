#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/KineticDelaunayFlipEvent.hpp"
#include "kinDS/Logger.hpp"
#include "kinDS/SegmentBuilderVisualDebug.hpp"
#include "kinDS/StrandTree.hpp"
#include "kinDS/TreeMesher.hpp"
#include "kinDS/VisualDebugHighlight.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

using namespace kinDS;
using Catch::Matchers::ContainsSubstring;

namespace
{
constexpr double kKnownBadFlipGlobalTime = 15.012518;
constexpr double kKnownBadFlipLocalTime = kKnownBadFlipGlobalTime - 15.0;

static void enable_all_log_levels_for_test()
{
  logger.setLogLevelMask(
    LogLevel::Debug | LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical);
}

static std::filesystem::path repoRootFromFixture(const char* relative_fixture)
{
  std::filesystem::path candidate = std::filesystem::path(relative_fixture);
  if (!std::filesystem::exists(candidate))
  {
    candidate = std::filesystem::path("..") / relative_fixture;
  }
  if (!std::filesystem::exists(candidate))
  {
    candidate = std::filesystem::path("../..") / relative_fixture;
  }
  if (!std::filesystem::exists(candidate))
  {
    throw std::runtime_error(std::string("test fixture not found: ") + relative_fixture);
  }

  std::filesystem::path directory = std::filesystem::absolute(candidate).parent_path();
  while (!directory.empty())
  {
    if (std::filesystem::exists(directory / "CMakeLists.txt") && std::filesystem::exists(directory / "kinDS"))
    {
      return directory;
    }
    directory = directory.parent_path();
  }

  throw std::runtime_error("could not locate kinDS repository root");
}

static std::filesystem::path fixturePath(const char* relative)
{
  return repoRootFromFixture(relative) / relative;
}

static std::filesystem::path oak2TestOutputRoot(const char* test_name)
{
  const std::filesystem::path output_root = repoRootFromFixture("tests/fixtures/test_oak_2_branches_1_2_sections_15_16.txt")
    / "tests" / "output" / test_name;
  // Start clean so assertions on generated files never pass on output left by an earlier run.
  std::filesystem::remove_all(output_root);
  std::filesystem::create_directories(output_root);
  return std::filesystem::absolute(output_root);
}

/// Visual debug limited to [section, section + 0.1] (the polynomial dump is gated on it) with a flip dump target.
static TreeMesher::Settings makeOak2DebugSettings(
  const std::filesystem::path& output_root, double flip_polynomial_dump_target_time)
{
  TreeMesher::Settings settings;
  settings.visual_debug = true;
  settings.visual_debug_output_root = output_root;
  settings.visual_debug_time_lower = std::floor(flip_polynomial_dump_target_time);
  settings.visual_debug_time_upper = std::floor(flip_polynomial_dump_target_time) + 0.1;
  settings.flip_polynomial_dump_target_time = flip_polynomial_dump_target_time;
  return settings;
}

static glm::dvec2 evaluateSitePolynomial(const Trajectory<2>& piece, double frac)
{
  return glm::dvec2(piece[0](frac), piece[1](frac));
}

static void requireFlipPolynomialMatchesGetPointAt(const KineticDelaunay& kd, double t)
{
  const StrandTree& tree = kd.getStrandTree();
  const size_t strand_count = tree.getPoints().size();
  const size_t section = static_cast<size_t>(std::floor(t));
  const double frac = t - static_cast<double>(section);
  if (frac < std::numeric_limits<double>::epsilon())
  {
    return;
  }

  double branch_lookup_time = t;
  if (branch_lookup_time <= static_cast<double>(section) + std::numeric_limits<double>::epsilon())
  {
    branch_lookup_time = static_cast<double>(section) + std::numeric_limits<double>::epsilon();
  }

  for (size_t strand_id = 0; strand_id < strand_count; ++strand_id)
  {
    const glm::dvec2 from_point_at = kd.getPointAt(strand_id, t);
    const size_t reference_branch = kd.getReferenceBranch(strand_id, branch_lookup_time);
    const Trajectory<2> piece = tree.getPiecePolynomial(strand_id, section, reference_branch);
    const glm::dvec2 from_polynomial = evaluateSitePolynomial(piece, frac);
    INFO("strand_id=" << strand_id << " t=" << t << " reference_branch=" << reference_branch);
    REQUIRE(from_point_at.x == Catch::Approx(from_polynomial.x).margin(1e-9));
    REQUIRE(from_point_at.y == Catch::Approx(from_polynomial.y).margin(1e-9));
  }
}

struct FlipSvgDebugHandler : public KineticDelaunay::EventCallback
{
  KineticDelaunay& kd;

  explicit FlipSvgDebugHandler(KineticDelaunay& kd_ref)
    : kd(kd_ref)
  {
  }

  void beforeEvent(KineticDelaunay::Event& event) override
  {
    auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&event);
    if (!flip)
    {
      return;
    }

    writeSegmentBuilderVisualDebugSvg(true, kd, kd.getGraph(), flip->occurrence_time, "before",
      "flip_he" + std::to_string(flip->half_edge_id), VisualDebugHighlight::forFlip(kd.getGraph(), flip->half_edge_id),
      kd.getRuntimeBranchIdForHalfEdge(flip->half_edge_id));
  }
};

static void configureOak2DebugOutput(KineticDelaunay& kd, const std::filesystem::path& output_root,
  double flip_polynomial_dump_target_time, std::optional<size_t> flip_polynomial_dump_target_half_edge = std::nullopt)
{
  kd.setVisualDebugOutputRoot(output_root);
  kd.setFlipPolynomialDumpTargetTime(flip_polynomial_dump_target_time);
  kd.setFlipPolynomialDumpTargetHalfEdge(flip_polynomial_dump_target_half_edge);
}

/// Records the first flip within 1e-4 of the known bad flip time (slice-local time), and the first flip overall.
struct BadFlipRecorder : public KineticDelaunay::EventCallback
{
  std::optional<std::pair<double, size_t>> flip_near_known_bad_time;
  std::optional<std::pair<double, size_t>> first_flip;

  void afterEvent(KineticDelaunay::Event& event) override
  {
    auto* flip = dynamic_cast<KineticDelaunay::FlipEvent*>(&event);
    if (!flip)
    {
      return;
    }
    const auto record = std::make_pair(flip->occurrence_time.real_time, flip->half_edge_id);
    if (!first_flip.has_value())
    {
      first_flip = record;
    }
    if (std::abs(record.first - kKnownBadFlipLocalTime) <= 1e-4 && !flip_near_known_bad_time.has_value())
    {
      flip_near_known_bad_time = record;
    }
  }
};

static StrandTree loadOak2BranchSlice()
{
  return StrandTree::loadFromFile(fixturePath("tests/fixtures/test_oak_2_branches_1_2_sections_15_16.txt"));
}

static void requirePolynomialDumpExists(const std::filesystem::path& output_root)
{
  const std::filesystem::path polynomials_file = output_root / "polynomials.txt";
  INFO("expected polynomials at: " << polynomials_file.string());
  REQUIRE(std::filesystem::exists(polynomials_file));
  std::ifstream in(polynomials_file);
  std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  REQUIRE_THAT(contents, ContainsSubstring("event_trigger(X)"));
  REQUIRE_THAT(contents, ContainsSubstring("inCircle"));
}

static void requireBranchesLogExists(const std::filesystem::path& output_root)
{
  const std::filesystem::path branches_file = output_root / "branches.txt";
  INFO("expected branches log at: " << branches_file.string());
  REQUIRE(std::filesystem::exists(branches_file));
  std::ifstream in(branches_file);
  std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  REQUIRE_THAT(contents, ContainsSubstring("event=create"));
}
} // namespace

TEST_CASE("oak_2 branches 1-2 sections 15-16 fixture loads", "[KineticDelaunay][Oak2Flip]")
{
  enable_all_log_levels_for_test();

  StrandTree tree = loadOak2BranchSlice();
  REQUIRE(tree.getHeight() == 1);
  REQUIRE(tree.getPoints().size() == 16);
}

TEST_CASE("oak_2 slice site polynomials match getPointAt around failing flip", "[KineticDelaunay][Oak2Flip]")
{
  enable_all_log_levels_for_test();

  StrandTree tree = loadOak2BranchSlice();

  constexpr double cutoff = 10.0;
  KineticDelaunay kd(tree, cutoff, false);
  kd.init();

  const std::vector<double> sample_times = {
    0.001,
    0.01,
    kKnownBadFlipLocalTime,
    0.05,
    0.25,
    0.5,
    0.75,
    1.0,
  };

  for (double t : sample_times)
  {
    requireFlipPolynomialMatchesGetPointAt(kd, t);
  }
}

TEST_CASE("oak_2 slice debug output near global t=15.012518", "[KineticDelaunay][Oak2Flip]")
{
  enable_all_log_levels_for_test();

  const std::filesystem::path output_root = oak2TestOutputRoot("oak2_slice");
  StrandTree tree = loadOak2BranchSlice();

  constexpr double cutoff = 10.0;
  KineticDelaunay kd(tree, cutoff, false);
  configureOak2DebugOutput(kd, output_root, kKnownBadFlipLocalTime);
  FlipSvgDebugHandler svg_handler(kd);
  kd.registerEventCallbacks(nullptr, &svg_handler, nullptr, nullptr);
  kd.init();

  REQUIRE_NOTHROW(kd.compute());

  // The bad flip this slice was extracted for no longer occurs (see the regression test below), so no flip debug
  // output (branch folders, SVGs) is written for it any more.
  for (const auto& entry : std::filesystem::recursive_directory_iterator(output_root))
  {
    const std::string filename = entry.path().filename().string();
    CHECK_FALSE((filename.find("flip") != std::string::npos && filename.find("0.0125") != std::string::npos));
  }
  INFO("slice debug output written to: " << output_root.string());
}

// Regression: the inCircle flip near global t=15.012518 was the bad flip this slice was extracted to explore. It must
// no longer be scheduled; fail if it comes back.
TEST_CASE("oak_2 slice has no bad inCircle flip near global t=15.012518", "[KineticDelaunay][Oak2Flip]")
{
  enable_all_log_levels_for_test();

  StrandTree tree = loadOak2BranchSlice();

  constexpr double cutoff = 10.0;
  KineticDelaunay kd(tree, cutoff, false);
  BadFlipRecorder recorder;
  kd.registerEventCallbacks(nullptr, &recorder, nullptr, nullptr);
  kd.init();

  REQUIRE_NOTHROW(kd.compute());
  if (recorder.flip_near_known_bad_time.has_value())
  {
    const auto [flip_time, half_edge_id] = *recorder.flip_near_known_bad_time;
    FAIL("bad flip came back at local t=" << flip_time << " half_edge_id=" << half_edge_id);
  }
}

// The slice has no flips any more, so the dump feature is exercised on the full tree: a first run records a flip that
// actually occurs, a second run targets exactly that flip. Plain KineticDelaunay (no SegmentBuilder) cannot handle the
// branch split, so both runs stop before it.
constexpr size_t kPlainKineticEndSection = 10;
TEST_CASE("test_oak_2 writes the flip polynomial dump for a flip that occurs", "[KineticDelaunay][Oak2Flip]")
{
  enable_all_log_levels_for_test();

  const std::filesystem::path tree_path = fixturePath("input_examples/test_oak_2.txt");
  constexpr double cutoff = 10.0;
  std::optional<std::pair<double, size_t>> target;
  {
    StrandTree tree = StrandTree::loadFromFile(tree_path);
    KineticDelaunay kd(tree, cutoff, false);
    kd.setSectionRange(0, kPlainKineticEndSection);
    BadFlipRecorder recorder;
    kd.registerEventCallbacks(nullptr, &recorder, nullptr, nullptr);
    kd.init();
    REQUIRE_NOTHROW(kd.compute());
    target = recorder.first_flip;
  }
  REQUIRE(target.has_value());

  const std::filesystem::path output_root = oak2TestOutputRoot("oak2_polynomials");
  StrandTree tree = StrandTree::loadFromFile(tree_path);
  KineticDelaunay kd(tree, cutoff, false);
  kd.setSectionRange(0, kPlainKineticEndSection);
  configureOak2DebugOutput(kd, output_root, target->first, target->second);
  kd.setVisualDebugEnabled(true);
  kd.init();
  REQUIRE_NOTHROW(kd.compute());

  INFO("target flip t=" << target->first << " half_edge_id=" << target->second);
  requirePolynomialDumpExists(output_root);
}

// Full-tree counterpart of the slice regression: with visual debug on at section 15 and the polynomial dump targeted
// at the known bad flip time, no dump may be written (that flip must not occur) and meshing must not throw.
TEST_CASE("full test_oak_2 meshing has no bad flip near t=15.012518", "[KineticDelaunay][Oak2Flip][integration]")
{
  enable_all_log_levels_for_test();

  const std::filesystem::path output_root = oak2TestOutputRoot("oak2_integration");
  StrandTree tree = StrandTree::loadFromFile(fixturePath("input_examples/test_oak_2.txt"));
  TreeMesher mesher(tree);
  mesher.setSettings(makeOak2DebugSettings(output_root, kKnownBadFlipGlobalTime));

  REQUIRE_NOTHROW(mesher.runMeshingAlgorithm(true));

  CHECK_FALSE(std::filesystem::exists(output_root / "polynomials.txt"));
  requireBranchesLogExists(output_root);
  INFO("integration debug output written to: " << output_root.string());
}
