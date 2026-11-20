#include "kinDS/KineticDelaunay.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace kinDS;

// A transition registered while another one hands over must start its last hold piece in the frame the previous
// piece ended in, so the section state is identical before and after the section event.
TEST_CASE("Hold pieces start in the handover frame only for the last hold piece", "[KineticDelaunay][FrameTransition]")
{
  PostSplitFrameTransition transition;
  transition.common_reference_branch = 6;
  transition.hold_end_height = 32;
  transition.start_reference_branches.emplace(52, 1);

  CHECK(holdPieceStartReference(transition, 52, 31) == 1);
  CHECK(holdPieceStartReference(transition, 52, 30) == 6);
  CHECK(holdPieceStartReference(transition, 53, 31) == 6);
}
