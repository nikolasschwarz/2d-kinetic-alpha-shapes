#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/HalfEdgeDelaunayGraph.hpp"
#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/StrandTree.hpp"
#include "kinDS/Logger.hpp"

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <bitset>

using namespace kinDS;

static void enable_all_log_levels_for_test()
{
  logger.setLogLevelMask(
    LogLevel::Debug | LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical);
}

// Export current state to SVG before running validation (for debugging failures).
static void exportSvgBeforeAssertion(const KineticDelaunay& kd, double t, const std::string& filename)
{
  std::vector<glm::dvec2> points = kd.getPointsAt(t);
  const HalfEdgeDelaunayGraph& graph = kd.getGraph();
  const size_t face_count = graph.getFaces().size();

  // Build Voronoi vertex -> containing triangle mapping from CrossingData.
  std::vector<size_t> voronoi_vertex_to_tri(face_count);
  constexpr size_t invalid_id = static_cast<size_t>(-1);
  for (size_t vid = 0; vid < face_count; ++vid)
  {
    voronoi_vertex_to_tri[vid] = kd.getCrossingDataContainingTriId(vid);
  }
  auto intersection_debug_data = kd.getCrossingIntersectionDebugData();

  HalfEdgeDelaunayGraphToSVG::write(
    points, graph, filename, 0.1, nullptr, true, &voronoi_vertex_to_tri, &intersection_debug_data); // margin 0.1, draw Voronoi edges in red
}

// Validate CrossingData: each Voronoi vertex's containing tri matches the list, and the vertex lies inside that triangle.
// The 'context' string is used to indicate where in the evolution we are (init, betweenSections, after event, etc.).
static void validateCrossingData(const KineticDelaunay& kd, double t, const std::string& context)
{
  const HalfEdgeDelaunayGraph& graph = kd.getGraph();
  const size_t face_count = graph.getFaces().size();
  const auto& faces = graph.getFaces();

  constexpr size_t invalid_id = static_cast<size_t>(-1);

  auto to_binary_ieee754 = [](double v) -> std::string
  {
    union
    {
      double d;
      std::uint64_t u;
    } u;
    u.d = v;
    std::uint64_t bits = u.u;
    int sign = static_cast<int>((bits >> 63) & 0x1);
    std::uint64_t exponent_bits = (bits >> 52) & 0x7FFu;
    std::uint64_t mantissa_bits = bits & ((1ull << 52) - 1);

    int exponent = static_cast<int>(exponent_bits) - 1023; // unbiased exponent

    std::bitset<11> exp_bits(exponent_bits);
    std::bitset<52> mantissa(mantissa_bits);

    std::ostringstream oss;
    oss << "sign=" << sign
        << ", exp_bits=" << exp_bits
        << ", exp=" << exponent
        << ", mantissa=" << mantissa;
    return oss.str();
  };

  auto edge_function = [](const glm::dvec2& a, const glm::dvec2& b, const glm::dvec2& c)
  {
    // Signed area * 2, orientation-agnostic containment will be done via barycentric weights.
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  };

  // Build convex hull polygon (outer boundary) once, using boundary edge iterator.
  // We keep the boundary half-edge IDs, the vertex indices, and their positions along the hull cycle.
  std::vector<size_t> hull_edge_ids;
  std::vector<size_t> hull_vertex_ids;
  std::vector<glm::dvec2> hull_points;
  {
    auto it = graph.boundaryEdgesBegin();
    auto end = graph.boundaryEdgesEnd();
    for (; it != end; ++it)
    {
      size_t he_id = *it;
      int v = graph.getHalfEdges()[he_id].origin;
      if (v == -1)
      {
        continue;
      }
      hull_edge_ids.push_back(he_id);
      hull_vertex_ids.push_back(static_cast<size_t>(v));
      hull_points.push_back(kd.getPointAt(static_cast<size_t>(v), t));
    }
  }

  auto isInsideConvexHull = [&](const glm::dvec2& p) -> bool
  {
    constexpr double eps = 1e-8;
    if (hull_points.size() < 3)
    {
      return false;
    }

    double sign = 0.0;
    const size_t n = hull_points.size();
    for (size_t i = 0; i < n; ++i)
    {
      const glm::dvec2& a = hull_points[i];
      const glm::dvec2& b = hull_points[(i + 1) % n];
      double cross = edge_function(a, b, p);
      if (std::abs(cross) < eps)
      {
        continue; // on the edge, treat as inside
      }
      if (sign == 0.0)
      {
        sign = cross;
      }
      else if (cross * sign < -eps)
      {
        // Different orientation => outside convex polygon
        return false;
      }
    }
    return true;
  };

  // For a given point p, find the triangle index that contains it geometrically (finite triangles only).
  auto findContainingFaceByGeometry = [&](const glm::dvec2& p) -> size_t
  {
    constexpr double eps = 1e-8;
    for (size_t f = 0; f < face_count; ++f)
    {
      const auto& tri = faces[f];
      auto tri_vertices = graph.adjacentTriangleVertices(tri.half_edges[0]);
      std::vector<glm::dvec2> pts;
      bool has_infinite_vertex = false;
      for (int v : tri_vertices)
      {
        if (v == -1)
        {
          has_infinite_vertex = true;
          break;
        }
        pts.push_back(kd.getPointAt(static_cast<size_t>(v), t));
      }
      if (has_infinite_vertex || pts.size() < 3)
      {
        continue;
      }

      const glm::dvec2& a = pts[0];
      const glm::dvec2& b = pts[1];
      const glm::dvec2& c = pts[2];

      double denom = edge_function(a, b, c);
      if (std::abs(denom) < eps)
      {
        continue; // degenerate
      }

      // Barycentric coordinates (orientation-agnostic).
      double w1 = edge_function(b, c, p) / denom;
      double w2 = edge_function(c, a, p) / denom;
      double w3 = 1.0 - w1 - w2;

      if (w1 >= -eps && w2 >= -eps && w3 >= -eps)
      {
        return f;
      }
    }
    return invalid_id;
  };

  auto format_triangle_and_barycentric = [&](size_t tri_id, const glm::dvec2& p) -> std::string
  {
    constexpr double eps = 1e-8;
    if (tri_id == invalid_id || tri_id >= face_count)
    {
      return " [tri_id invalid for barycentric debug]";
    }

    const auto& tri = faces[tri_id];
    auto tri_vertices = graph.adjacentTriangleVertices(tri.half_edges[0]);
    std::vector<glm::dvec2> pts;

    KINDS_DEBUG("Tri vertices for tri_id " << tri_id << ": " << tri_vertices[0] << ", " << tri_vertices[1] << ", " << tri_vertices[2] << " at t = " << t);
    for (int v : tri_vertices)
    {
      if (v == -1)
      {
        return " [triangle has infinite vertex; barycentric debug skipped]";
      }
      pts.push_back(kd.getPointAt(static_cast<size_t>(v), t));
    }
    if (pts.size() < 3)
    {
      return " [triangle has fewer than 3 finite vertices; barycentric debug skipped]";
    }

    const glm::dvec2& a = pts[0];
    const glm::dvec2& b = pts[1];
    const glm::dvec2& c = pts[2];

    double denom = edge_function(a, b, c);
    if (std::abs(denom) < eps)
    {
      return " [triangle is degenerate (denom≈0); barycentric debug skipped]";
    }

    double w1 = edge_function(b, c, p) / denom;
    double w2 = edge_function(c, a, p) / denom;
    double w3 = 1.0 - w1 - w2;

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(15);
    oss << " Voronoi vertex p=(" << p.x << "," << p.y << ")"
        << " tri_id=" << tri_id
        << " tri_verts=[(" << a.x << "," << a.y << "),("
        << b.x << "," << b.y << "),("
        << c.x << "," << c.y << ")]"
        << " barycentric=(w1=" << w1 << ", w2=" << w2 << ", w3=" << w3 << ")"
        << " | p.x_bin={" << to_binary_ieee754(p.x) << "}"
        << " p.y_bin={" << to_binary_ieee754(p.y) << "}"
        << " w1_bin={" << to_binary_ieee754(w1) << "}"
        << " w2_bin={" << to_binary_ieee754(w2) << "}"
        << " w3_bin={" << to_binary_ieee754(w3) << "}";
    return oss.str();
  };

  KINDS_DEBUG(context << ": validateCrossingData called at t=" << t << " with " << face_count << " faces.");

  for (size_t face_id = 0; face_id < face_count; ++face_id)
  {
    // Skip validation on infinite faces for now (faces with any vertex at infinity).
    const auto& tri_face = faces[face_id];
    auto tri_face_vertices = graph.adjacentTriangleVertices(tri_face.half_edges[0]);
    bool face_has_infinite_vertex = false;
    for (int v : tri_face_vertices)
    {
      if (v == -1)
      {
        face_has_infinite_vertex = true;
        break;
      }
    }
    if (face_has_infinite_vertex)
    {
      continue;
    }

    std::vector<size_t> voronoi_vertices = kd.getCrossingDataVoronoiVerticesInTri(face_id);

    for (size_t voronoi_vertex_id : voronoi_vertices)
    {
      size_t containing_tri = kd.getCrossingDataContainingTriId(voronoi_vertex_id);
      if (containing_tri == invalid_id)
      {
        throw std::runtime_error(
          context + ": CrossingData inconsistency: Voronoi vertex " + std::to_string(voronoi_vertex_id)
          + " has invalid containing triangle id (expected to be in face " + std::to_string(face_id) + ").");
      }
      if (containing_tri != face_id)
      {
        throw std::runtime_error(
          context + ": CrossingData inconsistency for Voronoi vertex " + std::to_string(voronoi_vertex_id)
          + ": tri_id_to_voronoi_vertices lists it in face " + std::to_string(face_id)
          + " but containing tri id is " + std::to_string(containing_tri) + ".");
      }

      glm::dvec3 pos_h = kd.getVoronoiVertexHomogeneous(voronoi_vertex_id, t);
      if (pos_h.z == 0.0)
      {
        continue;
      }

      glm::dvec2 p(pos_h.x / pos_h.z, pos_h.y / pos_h.z);

      // Independently find the triangle that geometrically contains this Voronoi vertex.
      size_t expected_face_by_geometry = findContainingFaceByGeometry(p);

      if (expected_face_by_geometry == invalid_id)
      {
        // No finite triangle contains this vertex. If it lies inside the convex hull, this is an error.
        if (isInsideConvexHull(p))
        {
          std::string msg = context + ": CrossingData geometric inconsistency for Voronoi vertex "
            + std::to_string(voronoi_vertex_id)
            + ": inside convex hull, but no finite triangle found that contains the vertex geometrically at t = "
            + std::to_string(t) + ".";
          msg += format_triangle_and_barycentric(containing_tri, p);
          throw std::runtime_error(msg);
        }

        // Outside hull and no finite containing triangle => considered inside an infinite triangle region.
        // Use the convex hull geometry to choose the infinite face according to the user's rule:
        // for hull vertex B with neighbors A,C, use the line through B orthogonal to AC, and then
        // take the face id from the corresponding boundary half-edge.
        if (hull_points.size() >= 3 && !hull_edge_ids.empty())
        {
          size_t n = hull_points.size();
          size_t best_idx = invalid_id;
          double best_proj = -std::numeric_limits<double>::infinity();
          for (size_t i = 0; i < n; ++i)
          {
            const glm::dvec2& A = hull_points[(i + n - 1) % n];
            const glm::dvec2& B = hull_points[i];
            const glm::dvec2& C = hull_points[(i + 1) % n];
            glm::dvec2 dirAC = C - A;
            double len2 = glm::dot(dirAC, dirAC);
            if (len2 < 1e-16)
            {
              continue;
            }
            glm::dvec2 n_i(-(dirAC.y), dirAC.x); // orthogonal to line through A and C
            n_i /= std::sqrt(len2); // normalize

            double proj = glm::dot(p - B, n_i);
            if (proj > best_proj)
            {
              best_proj = proj;
              best_idx = i;
            }
          }
          if (best_idx != invalid_id)
          {
            size_t he_id = hull_edge_ids[best_idx];
            // The face associated with this boundary half-edge represents the infinite triangle region.
            size_t inf_face = graph.getHalfEdges()[he_id].face;
            if (inf_face != invalid_id)
            {
              expected_face_by_geometry = inf_face;
            }
          }
        }

        // If we still couldn't classify, leave expected_face_by_geometry as invalid_id and
        // accept this as belonging to some infinite region for now.
        if (expected_face_by_geometry == invalid_id)
        {
          continue;
        }
      }

      // Always emit detailed DEBUG info for this Voronoi vertex, now that expected_face_by_geometry
      // is defined (for finite regions; infinite regions are handled by the early-continue above).
      {
        std::string dbg = context + ": CrossingData debug for Voronoi vertex "
          + std::to_string(voronoi_vertex_id)
          + " at t=" + std::to_string(t)
          + " containing_tri=" + std::to_string(containing_tri)
          + " geom_face=" + std::to_string(expected_face_by_geometry);
        dbg += format_triangle_and_barycentric(containing_tri, p);
        dbg += " geom_suggested:" + format_triangle_and_barycentric(expected_face_by_geometry, p);
        KINDS_DEBUG(dbg);
      }

      if (expected_face_by_geometry != containing_tri)
      {
        std::string msg = context + ": CrossingData geometric inconsistency for Voronoi vertex "
          + std::to_string(voronoi_vertex_id) + ": containing_tri=" + std::to_string(containing_tri)
          + " but geometric search suggests face " + std::to_string(expected_face_by_geometry) + " instead.";
        msg += format_triangle_and_barycentric(containing_tri, p);
        msg += " geom_suggested:" + format_triangle_and_barycentric(expected_face_by_geometry, p);
        throw std::runtime_error(msg);
      }
    }
  }
}

// Build a StrandTree similar to the demo data used in kinetic_delaunay_example()
static StrandTree makeDemoStrandTree()
{
  std::vector<glm::dvec2> trajectory_A = {
    { -0.432132, -0.426942 }, { -0.447292, -0.580708 }, { -0.469864, -0.531837 }, { -0.578741, -0.494280 },
    { -0.519044, -0.496727 }, { -0.487418, -0.587100 }, { -0.536664, -0.465019 }, { -0.536664, -0.465019 }
  };

  std::vector<glm::dvec2> trajectory_B = {
    { -0.150887, -0.424968 }, { -0.101774, -0.349936 }, { -0.052661, -0.274904 }, { -0.003548, -0.199872 },
    { 0.045565, -0.124840 },  { 0.094678, -0.049808 },  { 0.143791, 0.025224 },   { 0.143791, 0.025224 }
  };

  std::vector<glm::dvec2> trajectory_C = {
    { -0.048665, -0.333097 }, { -0.197330, -0.266194 }, { -0.345995, -0.199291 }, { -0.494660, -0.132388 },
    { -0.643325, -0.065485 }, { -0.791990, 0.001418 },  { -0.940656, 0.068321 },  { -0.940656, 0.068321 }
  };

  std::vector<glm::dvec2> trajectory_D = {
    { 0.467745, 0.111272 }, { 0.435490, 0.022544 }, { 0.403235, -0.066183 }, { 0.370980, -0.154910 },
    { 0.338725, -0.243637 }, { 0.306470, -0.332364 }, { 0.274215, -0.421091 }, { 0.274215, -0.421091 }
  };

  std::vector<std::vector<glm::dvec2>> support_points { trajectory_A, trajectory_B, trajectory_C, trajectory_D };

  // Same subdivisions pattern as in the demo
  std::vector<std::vector<double>> subdivisions
    = { { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 }, { 0.2, 0.4, 0.6, 0.8 } };

  // No explicit physics segment mapping needed for this test
  std::vector<std::vector<int>> physics_strand_to_segment_indices;

  // Identity transforms per height and branch (single branch)
  const size_t height = support_points.front().size();
  std::vector<std::vector<glm::dmat4>> transforms_by_height_and_branch(height, std::vector<glm::dmat4>(1, glm::dmat4(1.0)));

  // Single branch (0) for all strands at all heights
  std::vector<std::vector<size_t>> branch_indices(
    support_points.size(), std::vector<size_t>(support_points.front().size(), 0));

  // strands_by_branch_id[h][branch_id] = list of strand ids
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch_id;
  strands_by_branch_id.resize(height);
  for (size_t h = 0; h < height; ++h)
  {
    strands_by_branch_id[h].push_back({ 0, 1, 2, 3 }); // single branch containing all strands
  }

  return StrandTree(
    support_points, subdivisions, physics_strand_to_segment_indices, transforms_by_height_and_branch, branch_indices,
    strands_by_branch_id);
}

// Event handler that validates CrossingData during the kinetic Delaunay evolution
struct CrossingDataTestHandler : public KineticDelaunay::EventHandler
{
  KineticDelaunay& kd;

  explicit CrossingDataTestHandler(KineticDelaunay& kd)
    : kd(kd)
  {
  }

  void afterFlipEvent(KineticDelaunay::Event& e) override
  {
    std::string filename
      = "t" + std::to_string(e.time) + "_crossing_data_after_flip_he" + std::to_string(e.half_edge_id) + ".svg";
    exportSvgBeforeAssertion(kd, e.time, filename);
    std::string ctx
      = "afterFlipEvent(time=" + std::to_string(e.time) + ", half_edge_id=" + std::to_string(e.half_edge_id) + ")";
    validateCrossingData(kd, e.time, ctx);
  }

  void afterCrossingEvent(KineticDelaunay::Event& e) override
  {
    const auto& graph = kd.getGraph();
    size_t old_tri = graph.getHalfEdges()[e.half_edge_id].face;
    size_t new_tri = graph.getHalfEdges()[e.half_edge_id ^ 1].face;
    std::string filename = "t" + std::to_string(e.time) + "_crossing_data_after_crossing_v"
                           + std::to_string(e.voronoi_vertex_id) + "_"
                           + std::to_string(old_tri) + "_to_" + std::to_string(new_tri) + ".svg";
    exportSvgBeforeAssertion(kd, e.time, filename);
    std::string ctx = "afterCrossingEvent(time=" + std::to_string(e.time)
      + ", voronoi_vertex_id=" + std::to_string(e.voronoi_vertex_id) + ")";
    validateCrossingData(kd, e.time, ctx);
  }

  void betweenSections(size_t index) override
  {
    std::string filename = "t" + std::to_string(index) + "_crossing_data_between_sections_" + std::to_string(index)
      + ".svg";
    exportSvgBeforeAssertion(kd, static_cast<double>(index), filename);
    std::string ctx = "betweenSections(index=" + std::to_string(index) + ")";
    validateCrossingData(kd, static_cast<double>(index), ctx);
  }
};

TEST_CASE("KineticDelaunay CrossingData consistency on demo data", "[KineticDelaunay][CrossingData]")
{
  enable_all_log_levels_for_test();
  KINDS_DEBUG("Debug test: Starting test");
  KINDS_INFO("Info test: Starting test");

  StrandTree tree = makeDemoStrandTree();

  double cutoff = 10.0;
  bool add_dummy_boundary = false;

  KineticDelaunay kd(tree, cutoff, add_dummy_boundary);
  CrossingDataTestHandler handler(kd);

  kd.init();

  // Debug: dump initial CrossingData mapping voronoi_vertex_id -> containing triangle.
  const HalfEdgeDelaunayGraph& graph_after_init = kd.getGraph();
  const size_t face_count_after_init = graph_after_init.getFaces().size();
  for (size_t face_id = 0; face_id < face_count_after_init; ++face_id)
  {
    std::vector<size_t> voronoi_vertices = kd.getCrossingDataVoronoiVerticesInTri(face_id);
    for (size_t voronoi_vertex_id : voronoi_vertices)
    {
      size_t containing_tri = kd.getCrossingDataContainingTriId(voronoi_vertex_id);
      KINDS_DEBUG("Init CrossingData: voronoi_vertex " << voronoi_vertex_id << " listed in tri " << face_id
                                                       << ", containing_tri=" << containing_tri);
    }
  }

  // Export SVG and validate immediately after initialization at t = 0.
  exportSvgBeforeAssertion(kd, 0.0, "t0_after_init.svg");
  REQUIRE_NOTHROW(validateCrossingData(kd, 0.0, "after init"));

  // The test passes if validateCrossingData never throws during the full evolution
  REQUIRE_NOTHROW(kd.compute(handler));
}

