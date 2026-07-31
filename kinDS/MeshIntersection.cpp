#include "MeshIntersection.hpp"

#include "Logger.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef USE_CGAL
#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/intersections.h> // triangle–triangle intersection
#include <CGAL/version.h>
#endif

using namespace kinDS;

namespace
{
#ifdef USE_CGAL
// Only fill small cracks/holes on meshlets; large openings stay open and fail closedness checks.
constexpr std::size_t kMaxHoleHalfedgesToFill = 32;
// Boundary solid should be closed: fill any hole that triangulate_hole can handle.
constexpr std::size_t kMaxHoleHalfedgesToFillBoundary = std::numeric_limits<std::size_t>::max();
constexpr std::size_t kMaxLoggedBoundaryEdges = 1;

struct MeshReadinessReport
{
  std::size_t vertex_count = 0;
  std::size_t face_count = 0;
  std::size_t edge_count = 0;
  std::size_t boundary_edge_count = 0;
  std::size_t non_manifold_vertex_count = 0;
  std::size_t add_face_failures = 0;
  bool is_closed = false;
  bool self_intersects = false;
  bool self_intersect_checked = false;

  bool readyForBoolean() const
  {
    return is_closed && non_manifold_vertex_count == 0;
  }

  std::string summary() const
  {
    std::ostringstream oss;
    oss << "closed=" << (is_closed ? "true" : "false") << ", V=" << vertex_count << ", F=" << face_count
        << ", E=" << edge_count << ", boundary_edges=" << boundary_edge_count
        << ", non_manifold_vertices=" << non_manifold_vertex_count;
    if (add_face_failures > 0)
    {
      oss << ", add_face_failures=" << add_face_failures;
    }
    if (self_intersect_checked)
    {
      oss << ", self_intersects=" << (self_intersects ? "true" : "false");
    }
    return oss.str();
  }

  std::string failureReason() const
  {
    if (boundary_edge_count > 0 || !is_closed)
    {
      return "Input mesh is not closed (" + std::to_string(boundary_edge_count) + " boundary edges).";
    }
    if (non_manifold_vertex_count > 0)
    {
      return "Input mesh has non-manifold vertices (" + std::to_string(non_manifold_vertex_count) + ").";
    }
    return "Input mesh is not suitable for boolean intersection.";
  }
};

std::size_t countBoundaryEdges(const MeshCGAL_internal& mesh)
{
  std::size_t count = 0;
  for (const auto e : edges(mesh))
  {
    if (CGAL::is_border(e, mesh))
    {
      ++count;
    }
  }
  return count;
}

std::size_t countNonManifoldVertices(const MeshCGAL_internal& mesh)
{
  std::size_t count = 0;
  for (const auto v : vertices(mesh))
  {
    if (PMP::is_non_manifold_vertex(v, mesh))
    {
      ++count;
    }
  }
  return count;
}

MeshReadinessReport diagnoseMesh(
  const MeshCGAL_internal& mesh, std::size_t add_face_failures = 0, bool check_self_intersection = false)
{
  MeshReadinessReport report;
  report.vertex_count = mesh.number_of_vertices();
  report.face_count = mesh.number_of_faces();
  report.edge_count = mesh.number_of_edges();
  report.add_face_failures = add_face_failures;
  report.boundary_edge_count = countBoundaryEdges(mesh);
  report.non_manifold_vertex_count = countNonManifoldVertices(mesh);
  report.is_closed = CGAL::is_closed(mesh) && report.boundary_edge_count == 0;
  if (check_self_intersection)
  {
    report.self_intersect_checked = true;
    try
    {
      report.self_intersects = PMP::does_self_intersect(mesh);
    }
    catch (...)
    {
      report.self_intersects = true;
    }
  }
  return report;
}

std::string formatCorefineFailureDiag(const MeshReadinessReport& input_pre, const MeshReadinessReport& boundary_pre,
  const MeshReadinessReport& input_post, const MeshReadinessReport& boundary_post, const MeshReadinessReport& output_post)
{
  std::ostringstream oss;
  oss << "corefine_and_compute_intersection returned false (CGAL: output is not a manifold volume; "
         "not necessarily that inputs were open). "
      << "input_pre[" << input_pre.summary() << "], boundary_pre[" << boundary_pre.summary() << "], "
      << "input_post[" << input_post.summary() << "], boundary_post[" << boundary_post.summary() << "], "
      << "output[" << output_post.summary() << "]";
  return oss.str();
}

void logSampleBoundaryEdges(const MeshCGAL_internal& mesh, const std::string& suffix)
{
  const std::size_t total = countBoundaryEdges(mesh);
  std::size_t logged = 0;
  for (const auto e : edges(mesh))
  {
    if (!CGAL::is_border(e, mesh))
    {
      continue;
    }

    const auto h = halfedge(e, mesh);
    const auto v0 = source(h, mesh);
    const auto v1 = target(h, mesh);
    const auto& p0 = mesh.point(v0);
    const auto& p1 = mesh.point(v1);
    //KINDS_WARNING("  boundary edge " << logged << ": v" << v0 << " (" << p0 << ") -- v" << v1 << " (" << p1 << ")"
    //                                 << suffix);
    if (++logged >= kMaxLoggedBoundaryEdges)
    {
      if (total > logged)
      {
        //KINDS_WARNING("  ... and " << (total - logged) << " more boundary edge(s)" << suffix);
      }
      break;
    }
  }
}

std::size_t collectBoundaryCycle(
  MeshCGAL_internal::Halfedge_index start, const MeshCGAL_internal& mesh, std::vector<MeshCGAL_internal::Halfedge_index>& cycle)
{
  cycle.clear();
  auto h = start;
  do
  {
    cycle.push_back(h);
    h = next(h, mesh);
    // Guard against corrupt cycles.
    if (cycle.size() > mesh.number_of_halfedges())
    {
      break;
    }
  } while (h != start);
  return cycle.size();
}

/// Best-effort repair so corefinement can run. Returns false only on exception.
bool tryRepairMeshForIntersection(
  MeshCGAL_internal& mesh, std::string* log, std::size_t max_hole_halfedges = kMaxHoleHalfedgesToFill)
{
  std::ostringstream oss;
  try
  {
    const std::size_t removed_deg = PMP::remove_degenerate_faces(mesh);
    PMP::remove_isolated_vertices(mesh);
    const std::size_t stitched = PMP::stitch_borders(mesh);
    const std::size_t duplicated_nm = PMP::duplicate_non_manifold_vertices(mesh);

    std::size_t holes_filled = 0;
    std::size_t faces_from_holes = 0;
    std::size_t holes_skipped_too_large = 0;
    std::unordered_set<MeshCGAL_internal::Halfedge_index> visited;
    for (const auto h : halfedges(mesh))
    {
      if (!mesh.is_border(h) || visited.count(h) != 0)
      {
        continue;
      }

      std::vector<MeshCGAL_internal::Halfedge_index> cycle;
      collectBoundaryCycle(h, mesh, cycle);
      for (const auto hh : cycle)
      {
        visited.insert(hh);
      }

      if (cycle.size() < 3)
      {
        continue;
      }
      if (cycle.size() > max_hole_halfedges)
      {
        ++holes_skipped_too_large;
        continue;
      }

      const std::size_t before_faces = mesh.number_of_faces();
      // Return type is an output iterator in this CGAL version (not a face count).
      PMP::triangulate_hole(mesh, h, PMP::parameters::use_delaunay_triangulation(true));
      const std::size_t added = mesh.number_of_faces() - before_faces;
      if (added > 0)
      {
        ++holes_filled;
        faces_from_holes += added;
      }
    }

    oss << "repair: removed_degenerate_faces=" << removed_deg << ", stitched_borders=" << stitched
        << ", duplicated_non_manifold_vertices=" << duplicated_nm << ", holes_filled=" << holes_filled
        << " (+" << faces_from_holes << " faces)";
    if (holes_skipped_too_large > 0)
    {
      oss << ", holes_skipped_too_large=" << holes_skipped_too_large << " (max_halfedges=" << max_hole_halfedges
          << ")";
    }
    if (log)
    {
      *log = oss.str();
    }
    return true;
  }
  catch (const std::exception& e)
  {
    oss << "repair threw: " << e.what();
    if (log)
    {
      *log = oss.str();
    }
    return false;
  }
  catch (...)
  {
    oss << "repair threw an unknown exception";
    if (log)
    {
      *log = oss.str();
    }
    return false;
  }
}
#endif
} // namespace

#ifdef USE_CGAL
using SideTest = CGAL::Side_of_triangle_mesh<MeshCGAL_internal, Kernel>;

void debug_print_face_index_map(const MeshCGAL<size_t>& mesh, const std::string& name)
{
  using FaceDescriptor = MeshCGAL_internal::Face_index;

  std::cout << "=== Face index map dump: " << name << " ===\n";

  std::size_t count = 0;
  for (FaceDescriptor f : faces(mesh.mesh))
  {
    std::cout << "  Face " << count << " --> fmap[f] = " << mesh.fidx[f] << "\n";
    ++count;
  }

  std::cout << "Total faces: " << count << "\n";
  std::cout << "=== End dump ===\n";
}

void debug_print_face_origin_map(const MeshCGAL<Origin>& mesh, const std::string& name)
{
  using FaceDescriptor = MeshCGAL_internal::Face_index;

  std::cout << "=== Face origin map dump: " << name << " ===\n";

  std::size_t count = 0;
  for (FaceDescriptor f : faces(mesh.mesh))
  {
    auto val = mesh.fidx[f];
    // val is std::pair<int, size_t>

    std::cout << "  Face " << count << " --> origin = { mesh_id = " << val.mesh_index << ", face_idx = " << val.face_id
              << " }\n";

    ++count;
  }

  std::cout << "Total faces: " << count << "\n";
  std::cout << "=== End dump ===\n";
}

// Convert from std::vector-based mesh data to CGAL Surface_mesh
static MeshCGAL_internal vectorToCgalMesh(
  const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<size_t, 3>>& triangles)
{
  MeshCGAL_internal mesh;
  std::vector<MeshCGAL_internal::Vertex_index> vmap(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    vmap[i] = mesh.add_vertex(Point_CGAL(vertices[i][0], vertices[i][1], vertices[i][2]));
  }

  for (const auto& t : triangles)
    mesh.add_face(vmap[t[0]], vmap[t[1]], vmap[t[2]]);

  return mesh;
}

static MeshCGAL<Origin> voronoiMeshToCgalMesh(const VoronoiMesh& input_mesh,
  [[maybe_unused]] const std::vector<int>& neighbor_segments, int mesh_id = -1, std::size_t* add_face_failures = nullptr)
{
  MeshCGAL<Origin> output_mesh("f:origin", Origin { -1, 0 });
  if (add_face_failures)
  {
    *add_face_failures = 0;
  }

  auto& vertices = input_mesh.getVertices();
  std::vector<MeshCGAL_internal::Vertex_index> vmap(vertices.size());
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    vmap[i] = output_mesh.mesh.add_vertex(Point_CGAL(vertices[i][0], vertices[i][1], vertices[i][2]));
  }

  auto& triangles = input_mesh.getTriangles();
  for (size_t i = 0; i < triangles.size(); i += 3)
  {
    auto face_index = output_mesh.mesh.add_face(vmap[triangles[i]], vmap[triangles[i + 1]], vmap[triangles[i + 2]]);

    if (output_mesh.mesh.is_valid(face_index))
    {
      // Always record origin so hole-fill defaults (-1) stay distinguishable from real faces.
      output_mesh.fidx[face_index] = { mesh_id, i / 3 };
    }
    else
    {
      if (add_face_failures)
      {
        ++(*add_face_failures);
        if (*add_face_failures <= 5)
        {
          //KINDS_WARNING("Adding triangle no. " << (i / 3) << " failed (likely degenerate or non-manifold), ignored.");
        }
      }
    }
  }

  return output_mesh;
}

static double eps = 1e-12;

kinDS::MeshIntersection::MeshIntersection(const VoronoiMesh& static_mesh)
  : boundary_mesh_voronoi(static_mesh)
{
  // All neighbor segments are -1
  std::vector<int> neighbor_segments(boundary_mesh_voronoi.getTriangleCount(), -1);

  std::size_t add_face_failures = 0;
  boundary_mesh = voronoiMeshToCgalMesh(boundary_mesh_voronoi, neighbor_segments, 0, &add_face_failures);

  MeshReadinessReport readiness = diagnoseMesh(boundary_mesh.mesh, add_face_failures, true);
  if (add_face_failures > 0)
  {
    //KINDS_WARNING("Boundary CGAL conversion dropped " << add_face_failures << " triangle(s); " << readiness.summary());
  }

  if (!readiness.readyForBoolean())
  {
    KINDS_WARNING("Boundary mesh not ready for intersection (" << readiness.summary() << "). Attempting repair.");
    if (readiness.boundary_edge_count > 0)
    {
      logSampleBoundaryEdges(boundary_mesh.mesh, {});
    }

    std::string repair_log;
    const bool repaired =
      tryRepairMeshForIntersection(boundary_mesh.mesh, &repair_log, kMaxHoleHalfedgesToFillBoundary);
    KINDS_WARNING(repair_log);

    readiness = diagnoseMesh(boundary_mesh.mesh, add_face_failures, true);
    if (!repaired || !readiness.readyForBoolean())
    {
      if (readiness.boundary_edge_count > 0)
      {
        logSampleBoundaryEdges(boundary_mesh.mesh, {});
      }
      KINDS_ERROR("Boundary mesh still not closed/manifold after repair [" << readiness.summary()
                                                                          << "]. Classification and intersection may fail.");
    }
    else
    {
      KINDS_WARNING("Boundary repair succeeded [" << readiness.summary() << "].");
    }
  }
  else
  {
    KINDS_DEBUG("Boundary mesh ready for intersection [" << readiness.summary() << "].");
  }

  try
  {
    PMP::orient_to_bound_a_volume(boundary_mesh.mesh);
  }
  catch (const std::exception& e)
  {
    KINDS_ERROR("orient_to_bound_a_volume failed on boundary mesh: " << e.what());
  }
  catch (...)
  {
    KINDS_ERROR("orient_to_bound_a_volume failed on boundary mesh (unknown exception).");
  }

  // Build AABB tree only after repair/orientation so ClassifyMeshRelation uses the fixed mesh.
  tree = TreeCGAL(faces(boundary_mesh.mesh).first, faces(boundary_mesh.mesh).second, boundary_mesh.mesh);
  tree.build();
  tree.accelerate_distance_queries();
}
#else
kinDS::MeshIntersection::MeshIntersection(const VoronoiMesh& static_mesh)
{
  KINDS_ERROR("CGAL was not found, acceleration data structure could not be constructed!");
}
#endif

void interpolateProperties(
  const VoronoiMesh& original_mesh, VoronoiMesh& new_mesh, size_t original_face_id, std::array<size_t, 3> new_tri)
{
  std::array<size_t, 3> old_tri;
  std::array<glm::dvec3, 3> old_normals;

  bool interpolate_uv = true;
  for (size_t i = 0; i < 3; i++)
  {
    old_tri[i] = original_mesh.getTriangles()[3 * original_face_id + i];
    old_normals[i] = original_mesh.getNormal(3 * original_face_id + i);
    // KINDS_DEBUG("Old tri id " << i << ": " << old_tri[i]);
    if (!original_mesh.hasValidUVIndex(3 * original_face_id + i))
    {
      interpolate_uv = false;
    }
  }

  // KINDS_DEBUG("interpolate_uv: " << interpolate_uv);
  if (interpolate_uv)
  {
    new_mesh.getUVIndices().resize(new_mesh.getTriangleCount() * 3, -1);
  }

  // compute the barycentric coordinates of the new face with regard to the new one so we can interpolate properties
  for (size_t i = 0; i < 3; i++)
  {
    auto barycentric_coords
      = original_mesh.computeBarycentricCoordinates(original_face_id, new_mesh.getVertices()[new_tri[i]]);

    // compute new normals and UVs by interpolating from old mesh
    glm::dvec3 interpolated_normal = barycentric_coords[0] * old_normals[0] + barycentric_coords[1] * old_normals[1]
      + barycentric_coords[2] * old_normals[2];
    size_t normal_index = new_mesh.addNormal(interpolated_normal);

    if (interpolate_uv)
    {
      auto interpolated_uv
        = barycentric_coords[0] * original_mesh.getUVs()[original_mesh.getUVIndices()[3 * original_face_id]]
        + barycentric_coords[1] * original_mesh.getUVs()[original_mesh.getUVIndices()[3 * original_face_id + 1]]
        + barycentric_coords[2] * original_mesh.getUVs()[original_mesh.getUVIndices()[3 * original_face_id + 2]];

      size_t index = new_mesh.addUV(interpolated_uv);
      new_mesh.getUVIndices()[new_mesh.getTriangleCount() * 3 - 3 + i] = index;
    }
  }
}

std::pair<VoronoiMesh, std::vector<int>> MeshIntersection::Intersect(const VoronoiMesh& mesh,
  const std::vector<int>& neighbor_segments, std::optional<size_t> meshlet_index, bool* failed)
{
  std::pair<VoronoiMesh, std::vector<int>> ret_val {
    VoronoiMesh(std::vector<std::string> {}, NormalMode::PerTriangleCorner), {}
  };
  auto& [intersection_mesh, out_neighbor_segments] = ret_val;
  if (failed)
  {
    *failed = false;
  }

  const auto meshlet_suffix = [&]() -> std::string
  {
    if (!meshlet_index.has_value())
    {
      return {};
    }
    return " (meshlet index " + std::to_string(*meshlet_index) + ")";
  };

  auto mark_failed = [&](const std::string& message)
  {
    if (failed)
    {
      *failed = true;
    }
    KINDS_ERROR(message << meshlet_suffix());
  };

  if (mesh.getTriangles().empty())
  {
    KINDS_WARNING("The input is empty. Returning empty intersection mesh." << meshlet_suffix());
    return ret_val; // empty mesh
  }

#ifdef USE_CGAL

  std::size_t add_face_failures = 0;
  MeshCGAL<Origin> input_mesh = voronoiMeshToCgalMesh(mesh, neighbor_segments, 1, &add_face_failures);
  // we need a copy because the corefinement is destructive
  MeshCGAL<Origin> boundary_mesh_copy = boundary_mesh;
  MeshCGAL<Origin> output_mesh("f:origin", { -1, 0 });

  MeshReadinessReport readiness = diagnoseMesh(input_mesh.mesh, add_face_failures);
  if (add_face_failures > 0)
  {
    //KINDS_WARNING("CGAL conversion dropped " << add_face_failures << " triangle(s); " << readiness.summary()
    //                                         << meshlet_suffix());
  }

  if (!readiness.readyForBoolean())
  {
    KINDS_WARNING("Input mesh not ready for intersection (" << readiness.summary() << "). Attempting repair."
                                                           << meshlet_suffix());
    if (readiness.boundary_edge_count > 0)
    {
      logSampleBoundaryEdges(input_mesh.mesh, meshlet_suffix());
    }

    std::string repair_log;
    const bool repaired = tryRepairMeshForIntersection(input_mesh.mesh, &repair_log);
    KINDS_WARNING(repair_log << meshlet_suffix());

    readiness = diagnoseMesh(input_mesh.mesh, add_face_failures);
    if (!repaired || !readiness.readyForBoolean())
    {
      if (readiness.boundary_edge_count > 0)
      {
        logSampleBoundaryEdges(input_mesh.mesh, meshlet_suffix());
      }
      mark_failed("Intersection failed - " + readiness.failureReason() + " [" + readiness.summary() + "]");
      return ret_val; // empty mesh
    }

    KINDS_WARNING("Repair succeeded for intersection (" << readiness.summary() << ")." << meshlet_suffix());
  }

  // Capture pre-corefine state: corefine mutates both input meshes.
  MeshReadinessReport input_pre = diagnoseMesh(input_mesh.mesh, add_face_failures, true);
  const MeshReadinessReport boundary_pre = diagnoseMesh(boundary_mesh_copy.mesh, 0, true);

  if (input_pre.self_intersects)
  {
    mark_failed("Intersection failed - Mesh self-intersects. input[" + input_pre.summary() + "], boundary["
      + boundary_pre.summary() + "]");
    return ret_val; // empty mesh
  }

  RecordingVisitor visitor(boundary_mesh_copy, input_mesh, output_mesh);

  bool success = false;
  try
  {
    success = PMP::corefine_and_compute_intersection(
      boundary_mesh_copy.mesh, input_mesh.mesh, output_mesh.mesh, PMP::parameters::visitor(visitor));
  }
  catch (const std::exception& e)
  {
    const MeshReadinessReport input_post = diagnoseMesh(input_mesh.mesh, add_face_failures, true);
    const MeshReadinessReport boundary_post = diagnoseMesh(boundary_mesh_copy.mesh, 0, true);
    const MeshReadinessReport output_post = diagnoseMesh(output_mesh.mesh, 0, true);
    mark_failed(std::string("Intersection failed - exception: ") + e.what() + ". "
      + formatCorefineFailureDiag(input_pre, boundary_pre, input_post, boundary_post, output_post));
    return ret_val; // empty mesh
  }
  catch (...)
  {
    const MeshReadinessReport input_post = diagnoseMesh(input_mesh.mesh, add_face_failures, true);
    const MeshReadinessReport boundary_post = diagnoseMesh(boundary_mesh_copy.mesh, 0, true);
    const MeshReadinessReport output_post = diagnoseMesh(output_mesh.mesh, 0, true);
    mark_failed("Intersection failed - An exception was thrown. "
      + formatCorefineFailureDiag(input_pre, boundary_pre, input_post, boundary_post, output_post));
    return ret_val; // empty mesh
  }

  if (!success)
  {
    const MeshReadinessReport input_post = diagnoseMesh(input_mesh.mesh, add_face_failures, true);
    const MeshReadinessReport boundary_post = diagnoseMesh(boundary_mesh_copy.mesh, 0, true);
    const MeshReadinessReport output_post = diagnoseMesh(output_mesh.mesh, 0, true);
    mark_failed(
      "Intersection failed - " + formatCorefineFailureDiag(input_pre, boundary_pre, input_post, boundary_post, output_post));
    return ret_val; // empty mesh
  }

  for (const auto& v : output_mesh.mesh.vertices())
  {
    const Point_CGAL& p = output_mesh.mesh.point(v);
    intersection_mesh.addVertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
  }

  for (const auto& f : output_mesh.mesh.faces())
  {
    std::array<size_t, 3> tri;
    size_t idx = 0;
    for (const auto& v : CGAL::vertices_around_face(output_mesh.mesh.halfedge(f), output_mesh.mesh))
    {
      tri[idx++] = static_cast<size_t>(v);
    }
    size_t new_triangle_index = intersection_mesh.addTriangle(tri[0], tri[1], tri[2]);

    // Neighbor segment material tags:
    //  >= 0  shared face with another segment
    //  -1    open / cut / intersection face → interior material
    //  -2    exterior bark surface
    // Default to interior so faces from the intersection boundary (and synthetic repair faces)
    // use inner wood rather than bark.
    int neighbor_segment = -1;

    auto origin = output_mesh.fidx[f];
    if (origin.mesh_index == 1)
    {
      if (origin.face_id < neighbor_segments.size())
      {
        neighbor_segment = neighbor_segments[origin.face_id];
      }
      else
      {
        KINDS_WARNING("Invalid origin triangle index: " << origin.face_id << meshlet_suffix());
      }

      interpolateProperties(mesh, intersection_mesh, origin.face_id, tri);
    }
    else if (origin.mesh_index == 0)
    {
      // Boundary/intersection faces: keep neighbor_segment = -1 (interior material).
      interpolateProperties(boundary_mesh_voronoi, intersection_mesh, origin.face_id, tri);
    }
    else if (origin.mesh_index < 0)
    {
      // Synthetic faces from repair (e.g. small hole fills): interior material, no property interp.
    }
    else
    {
      KINDS_WARNING("Invalid origin mesh index: " << origin.mesh_index << meshlet_suffix());
    }

    out_neighbor_segments.push_back(neighbor_segment);
  }
#else
  mark_failed(
    "CGAL is required for the intersection computation but was not found. Returning empty intersection mesh.");
#endif

  return ret_val;
}

kinDS::MeshIntersection::MeshRelation kinDS::MeshIntersection::ClassifyMeshRelation(
  const VoronoiMesh& mesh, bool assume_inside)
{
#ifdef USE_CGAL
  // --- Side-of-mesh using the same tree (no rebuild) ---
  SideTest side(tree);

  bool any_inside = false;
  bool any_outside = false;

  // ---------- 1. Triangle-M0 intersection using the tree ----------

  for (size_t i = 0; i < mesh.getTriangles().size(); i += 3)
  {
    const glm::dvec3& p0 = mesh.getVertices()[mesh.getTriangles()[i]];
    Point_CGAL p0_cgal(p0[0], p0[1], p0[2]);
    const glm::dvec3& p1 = mesh.getVertices()[mesh.getTriangles()[i + 1]];
    Point_CGAL p1_cgal(p1[0], p1[1], p1[2]);
    const glm::dvec3& p2 = mesh.getVertices()[mesh.getTriangles()[i + 2]];
    Point_CGAL p2_cgal(p2[0], p2[1], p2[2]);

    CGAL::Triangle_3<Kernel> tri(p0_cgal, p1_cgal, p2_cgal);
    // This performs triangle vs. all boundary_mesh triangles intersection test via the tree
    if (tree.do_intersect(tri))
    {
      return MeshRelation::INTERSECTING;
    }
  }

  if (assume_inside)
  {
    return MeshRelation::INSIDE;
  }

  // ---------- 2. No intersections --> classify inside/outside ----------
  for (size_t i = 0; i < mesh.getTriangles().size(); i += 3)
  {
    const glm::dvec3& p = mesh.getVertices()[mesh.getTriangles()[i]];
    Point_CGAL p_cgal(p[0], p[1], p[2]);

    CGAL::Bounded_side bs = side(p_cgal);

    if (bs == CGAL::ON_BOUNDARY)
      return MeshRelation::INTERSECTING; // touching boundary = intersecting

    if (bs == CGAL::ON_BOUNDED_SIDE)
      any_inside = true;
    else
      any_outside = true;

    // If both occur, surface must cross the boundary (even without intersection)
    if (any_inside && any_outside)
      return MeshRelation::INTERSECTING;
  }

  if (any_inside)
    return MeshRelation::INSIDE;
  return MeshRelation::OUTSIDE;
#else
  KINDS_ERROR("CGAL was not found, cannot determine mesh relation!");
  return MeshRelation::UNDEFINED;
#endif
}

MatchResult MeshIntersection::MatchPointOnSurface(const glm::dvec3& p, double epsilon) const
{
#ifdef USE_CGAL
  const CGAL::Surface_mesh<Point_CGAL>& mesh = boundary_mesh.mesh;
  Point_CGAL query(p[0], p[1], p[2]);

  auto result = tree.closest_point_and_primitive(query);
  const Point_CGAL& closest = result.first;
  CGAL::Surface_mesh<Point_CGAL>::Face_index f = result.second;
  auto origin = boundary_mesh.fidx[f];

  double dist2 = CGAL::squared_distance(query, closest);
  if (dist2 > epsilon * epsilon)
    return {};

  // Extract triangle vertices from face
  auto h = halfedge(f, mesh);
  auto v0 = target(h, mesh);
  h = next(h, mesh);
  auto v1 = target(h, mesh);
  h = next(h, mesh);
  auto v2 = target(h, mesh);

  const Point_CGAL& a = mesh.point(v0);
  const Point_CGAL& b = mesh.point(v1);
  const Point_CGAL& c = mesh.point(v2);

  // Barycentric coordinates
  Kernel::Vector_3 v0v = b - a;
  Kernel::Vector_3 v1v = c - a;
  Kernel::Vector_3 v2v = query - a;

  double d00 = v0v * v0v;
  double d01 = v0v * v1v;
  double d11 = v1v * v1v;
  double d20 = v2v * v0v;
  double d21 = v2v * v1v;

  double denom = d00 * d11 - d01 * d01;
  if (std::abs(denom) < 1e-14)
    return {};

  double v = (d11 * d20 - d01 * d21) / denom;
  double w = (d00 * d21 - d01 * d20) / denom;
  double u = 1.0 - v - w;

  return { true, origin.face_id, u, v, w };
#endif
  return {};
}
