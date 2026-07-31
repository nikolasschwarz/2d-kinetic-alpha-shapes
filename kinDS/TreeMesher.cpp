#include "TreeMesher.hpp"
#include "IndexIterator.hpp"
#include "KineticDelaunay.hpp"
#include "Logger.hpp"
#include "ObjExporter.hpp"
#include "SegmentBuilder.hpp"
#include "TreeMesher.hpp"
#include "Validator.hpp"
#include <algorithm>
#include <execution>
#include <filesystem>
#include <mutex>
#include <sstream>
#include <stdexcept>

#ifdef _DEBUG
#pragma message("kinDS: TreeMesher.cpp built in DEBUG")
#endif

#ifdef _RELEASE
#pragma message("kinDS: TreeMesher.cpp built in RELEASE")
#endif

using namespace kinDS;

static std::vector<std::pair<size_t, double>> MergeSortedVectors(const std::vector<std::vector<double>>& inputs)
{
  using Entry = std::pair<size_t, double>; // (index of input vector, value)
  std::vector<std::pair<size_t, double>> result;

  struct HeapNode
  {
    size_t vec_idx; // which input vector
    size_t elem_idx; // index inside that vector
    double value; // value itself

    bool operator>(const HeapNode& other) const
    {
      return value > other.value; // for min-heap
    }
  };

  std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<>> min_heap;

  // Initialize heap with the first element of each vector
  for (size_t i = 0; i < inputs.size(); ++i)
  {
    if (!inputs[i].empty())
    {
      min_heap.push({ i, 0, inputs[i][0] });
    }
  }

  while (!min_heap.empty())
  {
    auto node = min_heap.top();
    min_heap.pop();

    // record (vector index, value)
    result.emplace_back(node.vec_idx, node.value);

    // advance in that vector
    if (node.elem_idx + 1 < inputs[node.vec_idx].size())
    {
      min_heap.push({ node.vec_idx, node.elem_idx + 1, inputs[node.vec_idx][node.elem_idx + 1] });
    }
  }

  return result;
}

kinDS::TreeMesher::TreeMesher(StrandTree& strand_tree)
  : strand_tree(strand_tree)
  , parallel_for(
      [&](size_t count, std::function<void(size_t)> func)
      {
        std::for_each(
          std::execution::par_unseq, IndexIterator { 0 }, IndexIterator { count }, [&](size_t i) { func(i); });
      })
{
}

TreeMesher::TreeMesher(StrandTree& strand_tree, std::function<void(size_t, std::function<void(size_t)>)> parallel_for)
  : strand_tree(strand_tree)
  , parallel_for(parallel_for)
{
}

void TreeMesher::exportMeshlets(MeshletExportMode export_mode, const std::filesystem::path& export_path,
  std::optional<bool> transformed, std::optional<size_t> max_exports) const
{
  if (!mesh_builder)
  {
    throw std::runtime_error("exportMeshlets: meshing has not been run yet.");
  }

  const bool apply_export_transform = transformed.value_or(!settings.transform_mesh_at_construction);
  // Explicit profile-space export (@c transformed == false), e.g. CLI --untransformed — not merely
  // "skip export transform because meshlets were already transformed at construction".
  const bool untransformed_export = transformed.has_value() && !transformed.value();
  const bool include_metadata = mesh_builder->store_mesh_metadata;

  std::vector<VoronoiMesh> meshlets_to_export;
  std::vector<std::string> export_suffixes;
  if (export_mode == MeshletExportMode::Raw)
  {
    std::vector<std::vector<int>> neighbors_unused;
    std::tie(meshlets_to_export, neighbors_unused) = mesh_builder->extractSegmentMeshlets(false);
    export_suffixes = mesh_builder->extractSegmentMeshletExportSuffixes(false);
  }
  else
  {
    meshlets_to_export = segment_meshlets;
    export_suffixes = segment_meshlet_export_suffixes;
  }

  const size_t export_count = max_exports.has_value()
    ? std::min(max_exports.value(), meshlets_to_export.size())
    : meshlets_to_export.size();

  auto meshlet_for_export = [&](size_t i) -> VoronoiMesh
  {
    VoronoiMesh mesh = meshlets_to_export[i];
    if (apply_export_transform)
    {
      const size_t strand_id = export_mode == MeshletExportMode::Raw
        ? mesh_builder->strandIdForRawMeshlet(i)
        : mesh_builder->strandIdForSegment(i);
      transformToWorldSpace(mesh, strand_id);
    }
    return mesh;
  };

  auto write_meshlet_obj = [&](const VoronoiMesh& mesh, const std::filesystem::path& path)
  {
    const bool include_vertex_colors = Validator::meshUsesValidationErrorMaterial(mesh);
    kinDS::ObjExporter::writeMesh(mesh, path, 1.0, 1.0, {}, include_metadata, include_vertex_colors,
      settings.alternate_section_shading);
  };

  const std::string untransformed_suffix = untransformed_export ? "_untransformed" : "";

  auto meshlet_filename = [&](size_t i) -> std::filesystem::path
  {
    const std::string suffix = (i < export_suffixes.size()) ? export_suffixes[i] : "";
    return std::filesystem::path("meshlet" + std::to_string(i) + suffix
      + meshlets_to_export[i].creationKineticTimeFilenameSuffix() + untransformed_suffix + ".obj");
  };

  auto with_untransformed_path_suffix = [&](std::filesystem::path path) -> std::filesystem::path
  {
    if (!untransformed_export)
    {
      return path;
    }
    const std::filesystem::path ext = path.extension();
    path.replace_extension();
    path += untransformed_suffix;
    if (!ext.empty())
    {
      path += ext;
    }
    return path;
  };

  switch (export_mode)
  {
  case MeshletExportMode::PerSegment:
  case MeshletExportMode::Raw:
  {
    std::filesystem::create_directories(export_path);
    for (size_t i = 0; i < export_count; ++i)
    {
      write_meshlet_obj(meshlet_for_export(i), export_path / meshlet_filename(i));
    }
    KINDS_DEBUG("Exported " << export_count << " meshlet OBJ file(s) (mode=" << meshletExportModeToString(export_mode)
                            << ", transformed=" << apply_export_transform << ") to " << export_path.string() << ".");
    return;
  }

  case MeshletExportMode::Combined:
    break;
  }

  std::filesystem::path obj_path = with_untransformed_path_suffix(export_path);
  if (obj_path.extension().empty())
  {
    obj_path.replace_extension(".obj");
  }
  if (obj_path.has_parent_path())
  {
    std::filesystem::create_directories(obj_path.parent_path());
  }

  kinDS::VoronoiMesh combined_mesh;
  bool combined_mesh_initialized = false;
  // Empty meshlets (e.g. unused segments) keep default NoNormals; lock combined mode from the first
  // meshlet that actually has geometry so the first += does not mismatch.
  NormalMode combined_normal_mode = NormalMode::NoNormals;
  for (size_t i = 0; i < export_count; ++i)
  {
    const VoronoiMesh& candidate = meshlets_to_export[i];
    if (candidate.getVertexCount() > 0 || candidate.getTriangleCount() > 0)
    {
      combined_normal_mode = candidate.getNormalMode();
      break;
    }
    if (i == 0)
    {
      combined_normal_mode = candidate.getNormalMode();
    }
  }
  for (size_t i = 0; i < export_count; ++i)
  {
    VoronoiMesh mesh = meshlet_for_export(i);
    const bool mesh_has_geometry = mesh.getVertexCount() > 0 || mesh.getTriangleCount() > 0;

    if (!combined_mesh_initialized)
    {
      combined_mesh = kinDS::VoronoiMesh(SegmentBuilder::MeshletExportMaterialNames, combined_normal_mode);
      // One OBJ group per meshlet index; boundaries are managed here, not via per-meshlet group_offsets.
      combined_mesh.setGroupOffsets({ 0 });
      combined_mesh.setGroupNames({ "meshlet_0" });
      combined_mesh_initialized = true;
    }
    else
    {
      combined_mesh.startNewGroup("meshlet_" + std::to_string(i));
    }

    if (mesh_has_geometry)
    {
      // operator+= merges the rhs group_offsets, which would shift/absorb prior meshlet groups.
      mesh.setGroupOffsets({});
      mesh.setGroupNames({});
      combined_mesh += mesh;
    }
  }

  if (!combined_mesh_initialized)
  {
    KINDS_DEBUG("exportMeshlets: no meshlets to export.");
    return;
  }

  kinDS::ObjExporter::writeMesh(combined_mesh, obj_path, 1.0, 1.0, {}, include_metadata,
    Validator::meshUsesValidationErrorMaterial(combined_mesh), settings.alternate_section_shading);
  const size_t group_count = combined_mesh.getGroupOffsets().size();
  KINDS_DEBUG("Exported combined mesh (" << group_count << " group(s) from " << export_count
                                       << " meshlet(s), mode=combined, transformed=" << apply_export_transform
                                       << ") to " << obj_path.string() << ".");
}

void TreeMesher::truncateToBoundary(const VoronoiMesh& boundary_mesh)
{
  // intersect all meshes with the boundary mesh and save the result
  // Build an AABB-tree of the boundary-mesh to prefilter
  kinDS::MeshIntersection boundary_intersector(boundary_mesh);

  /*ProgressBar intersection_progress_bar(
    0, segment_meshlets.size(), "Computing Mesh Intersections", ProgressBar::Display::Absolute, 50);

  std::atomic<int> progress_counter { 0 };

  std::atomic<int> threads { 0 };
  thread_local bool counted = false;*/

  std::mutex failed_mutex;
  std::vector<std::pair<size_t, VoronoiMesh>> failed_meshlets;

  parallel_for(segment_meshlets.size(),
    [&](auto mesh_index)
    {
      // progress_counter.fetch_add(1, std::memory_order_relaxed);
      // intersection_progress_bar.Update(progress_counter);
      // Do not assume inside: classify fully so outside meshlets are cleared correctly.
      auto intersect_relation = boundary_intersector.ClassifyMeshRelation(segment_meshlets[mesh_index], false);

      switch (intersect_relation)
      {
      case kinDS::MeshIntersection::MeshRelation::INSIDE:
        // do nothing
        break;

      case kinDS::MeshIntersection::MeshRelation::INTERSECTING:
        if (settings.debug_export_meshes && mesh_index < settings.max_meshlet_export)
        {
          const std::string suffix
            = (mesh_index < segment_meshlet_export_suffixes.size()) ? segment_meshlet_export_suffixes[mesh_index] : "";
          kinDS::ObjExporter::writeMesh(segment_meshlets[mesh_index],
            "meshlet" + std::to_string(mesh_index) + suffix
              + segment_meshlets[mesh_index].creationKineticTimeFilenameSuffix() + "_raw.obj",
            1.0, 1.0, {}, mesh_builder ? mesh_builder->store_mesh_metadata : settings.store_mesh_metadata);
        }
        {
          bool intersection_failed = false;
          auto [clipped_mesh, clipped_neighbors] = boundary_intersector.Intersect(
            segment_meshlets[mesh_index], meshing_neighbor_indices[mesh_index], mesh_index, &intersection_failed);
          if (intersection_failed)
          {
            std::lock_guard<std::mutex> lock(failed_mutex);
            failed_meshlets.emplace_back(mesh_index, segment_meshlets[mesh_index]);
            if (!settings.keep_original_on_intersection_failure)
            {
              segment_meshlets[mesh_index] = kinDS::VoronoiMesh();
              meshing_neighbor_indices[mesh_index] = {};
            }
            // else keep the uncut meshlet and its neighbor indices.
          }
          else
          {
            segment_meshlets[mesh_index] = std::move(clipped_mesh);
            meshing_neighbor_indices[mesh_index] = std::move(clipped_neighbors);
          }
        }
        break;

      case kinDS::MeshIntersection::MeshRelation::OUTSIDE:
        // fully outside, result is empty mesh
        segment_meshlets[mesh_index] = kinDS::VoronoiMesh();
        meshing_neighbor_indices[mesh_index] = {};
        break;

      case kinDS::MeshIntersection::MeshRelation::UNDEFINED:
        KINDS_ERROR("Mesh relation returned UNDEFINED (meshlet index " << mesh_index << ")");
        break;

      default:
        KINDS_ERROR("Unknown return value of mesh relation (meshlet index " << mesh_index << ")");
        break;
      }
    });

  if (!failed_meshlets.empty())
  {
    std::sort(failed_meshlets.begin(), failed_meshlets.end(),
      [](const auto& a, const auto& b) { return a.first < b.first; });

    std::ostringstream index_list;
    for (size_t i = 0; i < failed_meshlets.size(); ++i)
    {
      if (i > 0)
      {
        index_list << ", ";
      }
      index_list << failed_meshlets[i].first;
    }
    KINDS_ERROR("Intersection failed for " << failed_meshlets.size()
                                           << " meshlet(s). Indices: [" << index_list.str() << "]. "
                                           << "Dumping failed meshlets as OBJ.");

    const bool include_metadata = mesh_builder ? mesh_builder->store_mesh_metadata : settings.store_mesh_metadata;
    for (const auto& [mesh_index, failed_mesh] : failed_meshlets)
    {
      const std::string suffix
        = (mesh_index < segment_meshlet_export_suffixes.size()) ? segment_meshlet_export_suffixes[mesh_index] : "";
      const std::filesystem::path obj_path = "failed_meshlet_" + std::to_string(mesh_index) + suffix
        + failed_mesh.creationKineticTimeFilenameSuffix() + ".obj";
      kinDS::ObjExporter::writeMesh(failed_mesh, obj_path, 1.0, 1.0, {}, include_metadata);
      KINDS_ERROR("Wrote failed meshlet " << mesh_index << " to " << obj_path.string());
    }
  }

  // intersection_progress_bar.Finish();
  // std::cout << "Threads used: " << threads.load() << "\n";
  //  Find empty meshes and try to fix them by expanding neighboring meshes
  if (settings.fix_missing_meshes)
  {
    fixFailedSegments(boundary_intersector);
  }
}

void TreeMesher::fixFailedSegments(const MeshIntersection& boundary_intersector)
{
  std::vector<size_t> empty_mesh_indices;
  for (size_t mesh_index = 0; mesh_index < segment_meshlets.size(); mesh_index++)
  {
    if (segment_meshlets[mesh_index].getVertexCount() == 0)
    {
      empty_mesh_indices.push_back(mesh_index);
    }
  }

  // go through all triangles of all meshes and check their neighbor indices
  // if a neighbor index corresponds to an empty mesh, try to copy the triangle into the corresponding empty mesh
  for (size_t mesh_index = 0; mesh_index < segment_meshlets.size(); mesh_index++)
  {
    if (std::binary_search(empty_mesh_indices.begin(), empty_mesh_indices.end(), mesh_index))
    {
      continue;
    }

    auto& mesh = segment_meshlets[mesh_index];
    auto& neighbor_indices = meshing_neighbor_indices[mesh_index];
    const auto& triangles = mesh.getTriangles();

    for (size_t triangle_index = 0; triangle_index < triangles.size(); triangle_index += 3)
    {
      int neighbor_mesh_index = neighbor_indices[triangle_index / 3];
      if (neighbor_mesh_index >= 0)
      {
        // indices are sorted, so we can use binary search
        if (std::binary_search(empty_mesh_indices.begin(), empty_mesh_indices.end(), neighbor_mesh_index))
        {
          kinDS::VoronoiMesh& neighbor_mesh = segment_meshlets[neighbor_mesh_index];

          size_t v0_index = triangles[triangle_index];
          size_t v1_index = triangles[triangle_index + 1];
          size_t v2_index = triangles[triangle_index + 2];

          glm::dvec3 v0 = mesh.getVertices()[v0_index];
          glm::dvec3 v1 = mesh.getVertices()[v1_index];
          glm::dvec3 v2 = mesh.getVertices()[v2_index];

          glm::dvec3 t0 = mesh.getUV(triangle_index);
          glm::dvec3 t1 = mesh.getUV(triangle_index + 1);
          glm::dvec3 t2 = mesh.getUV(triangle_index + 2);

          // TODO: normals, uvs, other
          size_t new_v0_index = neighbor_mesh.addVertex(v0);
          size_t new_v1_index = neighbor_mesh.addVertex(v1);
          size_t new_v2_index = neighbor_mesh.addVertex(v2);

          size_t new_t0_index = neighbor_mesh.addUV(t0);
          size_t new_t1_index = neighbor_mesh.addUV(t1);
          size_t new_t2_index = neighbor_mesh.addUV(t2);

          // invert order to maintain consistent orientation
          neighbor_mesh.addTriangle(new_v1_index, new_v0_index, new_v2_index, new_t1_index, new_t0_index, new_t2_index);
          meshing_neighbor_indices[neighbor_mesh_index].push_back(mesh_index);
        }
      }
    }
  }

  // merge vertices and check how many could be (partially) fixed
  const auto& boundary_mesh = getBoundaryMesh();
  size_t fixed_mesh_count = 0;
  for (size_t mesh_index : empty_mesh_indices)
  {
    auto& mesh = segment_meshlets[mesh_index];

    if (mesh.getVertexCount() == 0)
    {
      continue;
    }

    mesh.mergeDuplicateVertices(1e-6);
    mesh.removeDegenerateTriangles();
    mesh.removeIsolatedVertices();
    mesh.computeNormals(kinDS::PerTriangleCorner);

    mesh.patchHoles(
      [&](size_t tri_index)
      {
        meshing_neighbor_indices[mesh_index].push_back(-2);
        // TODO: copy from neighbor mesh
      },
      [&](size_t v_index, size_t corner_index)
      {
        auto& v = mesh.getVertices()[v_index];
        // query point in boundary mesh
        auto match = boundary_intersector.MatchPointOnSurface(v);
        if (!match.hit)
        {
          //   mark as not bark
          //   meshing_neighbor_indices[mesh_index][tri_index] = -1;

          // get properties from neighbor mesh
          std::vector<size_t> corner_indices = mesh.findTriangleCorners(v_index);
          if (corner_indices.empty())
          {
            KINDS_ERROR("Could not find neighboring vertex!");
            mesh.setUV({ 0, 0, 0 }, corner_index);
            mesh.setNormal({ 0, 0, 0 }, corner_index);
            return;
          }

          size_t neighbor_corner_index = corner_indices[0];
          auto uv = mesh.getUV(neighbor_corner_index);

          // the UV we got here is not in polar coordinates that are suitable for the bark
          glm::dvec2 coords { uv[0], uv[1] };
          double angle = std::atan2(coords[1] - 0.5, coords[0] - 0.5);

          glm::dvec3 new_uv { angle / (2 * glm::pi<double>()), uv[2], uv[2] };

          mesh.setUV(new_uv, corner_index);
          // compute normal from the triangle, we don't have better information here
          size_t triangle_index = neighbor_corner_index / 3;
          size_t t0_index = mesh.getTriangles()[triangle_index * 3];
          size_t t1_index = mesh.getTriangles()[triangle_index * 3 + 1];
          size_t t2_index = mesh.getTriangles()[triangle_index * 3 + 2];
          auto p0 = mesh.getVertices()[t0_index];
          auto p1 = mesh.getVertices()[t1_index];
          auto p2 = mesh.getVertices()[t2_index];
          auto n = glm::normalize(-glm::cross(p1 - p0, p2 - p0));

          mesh.setNormal(n, corner_index);
        }
        else
        {
          size_t matched_triangle_index = match.triangle_index;

          // use barycentric coordinates to interpolate normal and uv from boundary mesh

          // first get boundary triangle vertices
          const auto& boundary_triangles = boundary_mesh.getTriangles();
          size_t bt0_index = matched_triangle_index * 3;
          size_t bt1_index = matched_triangle_index * 3 + 1;
          size_t bt2_index = matched_triangle_index * 3 + 2;

          // get normals
          auto n0 = boundary_mesh.getNormal(bt0_index);
          auto n1 = boundary_mesh.getNormal(bt1_index);
          auto n2 = boundary_mesh.getNormal(bt2_index);

          // get uvs
          auto uv0 = boundary_mesh.getUV(bt0_index);
          auto uv1 = boundary_mesh.getUV(bt1_index);
          auto uv2 = boundary_mesh.getUV(bt2_index);

          // interpolate
          auto n = n0 * match.u + n1 * match.v + n2 * match.w;
          n = glm::normalize(n);

          auto uv = uv0 * match.u + uv1 * match.v + uv2 * match.w;
          mesh.setUV(uv, corner_index);
          mesh.setNormal(n, corner_index);
        }
      });
    fixed_mesh_count++;
  }

  KINDS_DEBUG("Fixed " << fixed_mesh_count << " out of " << empty_mesh_indices.size() << " meshes.");
}

std::pair<std::vector<float>, std::vector<float>> TreeMesher::computeTopAndBottomBoundaryDistances(
  const std::vector<float>& boundary_distance_by_segment_id)
{
  std::vector<float> bottom_boundary_distances_by_strand_id(strand_tree.getPhysicsStrandToSegmentIndices().size());
  std::vector<float> top_boundary_distances_by_strand_id(strand_tree.getPhysicsStrandToSegmentIndices().size());

  for (size_t strand_id = 0; strand_id < strand_tree.getPhysicsStrandToSegmentIndices().size(); strand_id++)
  {
    int bottom_segment_id = strand_tree.getPhysicsStrandToSegmentIndices()[strand_id].front();
    int top_segment_id = strand_tree.getPhysicsStrandToSegmentIndices()[strand_id].back();

    bottom_boundary_distances_by_strand_id[strand_id] = boundary_distance_by_segment_id[bottom_segment_id];
    top_boundary_distances_by_strand_id[strand_id] = boundary_distance_by_segment_id[top_segment_id];
  }

  return std::make_pair(bottom_boundary_distances_by_strand_id, top_boundary_distances_by_strand_id);
}

void TreeMesher::runKineticDelaunay(bool visual_debug)
{
  KINDS_DEBUG("Starting Kinetic Delaunay Voronoi Meshing...");
  // sort subdivisions into a single array
  std::vector<std::pair<size_t, double>> subdivisions = MergeSortedVectors(strand_tree.getSubdivisionsByStrand());

  kinetic_delaunay = std::make_shared<KineticDelaunay>(strand_tree, settings.alpha_cutoff, false);
  // Debug/error file dumps key off the visual-debug output root; default to cwd when either is enabled.
  const bool enable_debug_output_root = visual_debug || settings.error_files;
  if (enable_debug_output_root && !settings.visual_debug_output_root.has_value())
  {
    kinetic_delaunay->setVisualDebugOutputRoot(std::filesystem::current_path());
  }
  else if (settings.visual_debug_output_root.has_value())
  {
    kinetic_delaunay->setVisualDebugOutputRoot(*settings.visual_debug_output_root);
  }
  kinetic_delaunay->setVisualDebugEnabled(visual_debug);
  kinetic_delaunay->setErrorFilesEnabled(settings.error_files || visual_debug);
  kinetic_delaunay->setVisualDebugSeparatePendingSplits(settings.visual_debug_separate_pending_splits);
  if (settings.flip_polynomial_dump_target_time.has_value())
  {
    kinetic_delaunay->setFlipPolynomialDumpTargetTime(settings.flip_polynomial_dump_target_time);
  }
  if (settings.flip_polynomial_dump_target_half_edge.has_value())
  {
    kinetic_delaunay->setFlipPolynomialDumpTargetHalfEdge(settings.flip_polynomial_dump_target_half_edge);
  }
  kinetic_delaunay->setComponentSplitPolicy(settings.retriangulate_on_component_split
      ? KineticDelaunay::ComponentSplitPolicy::Retriangulate
      : KineticDelaunay::ComponentSplitPolicy::InPlaceCut);
  kinetic_delaunay->setSectionRange(settings.start_section, settings.end_section);

  const bool transform_mesh_at_construction = settings.transform_mesh_at_construction;

  mesh_builder = std::make_shared<SegmentBuilder>(*kinetic_delaunay, subdivisions, transform_mesh_at_construction,
    visual_debug, parallel_for);
  mesh_builder->setErrorFiles(settings.error_files || visual_debug);
  // Validate needs source metadata; otherwise honor Settings::store_mesh_metadata (CLI-controllable).
  if (settings.validate_mesh_vertex_sources && !settings.store_mesh_metadata)
  {
    KINDS_WARNING("validate_mesh_vertex_sources requires store_mesh_metadata; enabling metadata storage");
  }
  mesh_builder->store_mesh_metadata = settings.store_mesh_metadata || settings.validate_mesh_vertex_sources;
  mesh_builder->validate_mesh_vertex_sources = settings.validate_mesh_vertex_sources;
  if (settings.validate_mesh_vertex_sources)
  {
    Validator::setLogFile(settings.validate_mesh_vertex_sources_log_path);
  }
  mesh_builder->diagnostics = settings.diagnostics;
  mesh_builder->mesh_cap_at_start = settings.mesh_cap_at_start;
  mesh_builder->mesh_cap_at_end = settings.mesh_cap_at_end;
  kinetic_delaunay->setDiagnosticsEnabled(settings.diagnostics);
  kinetic_delaunay->setSitesInsideConvexHullCheckEnabled(settings.check_sites_inside_convex_hull);

  KINDS_INFO("Starting Kinetic Delaunay Voronoi Meshing with settings: alpha_cutoff=" << settings.alpha_cutoff
                                                                                      << ", visual_debug=" << visual_debug
                                                                                      << ", error_files=" << mesh_builder->shouldDumpErrorFiles()
                                                                                      << ", transform_mesh_at_construction="
                                                                                      << transform_mesh_at_construction
                                                                                      << ", debug_export_meshes="
                                                                                      << settings.debug_export_meshes
                                                                                      << ", max_meshlet_export="
                                                                                      << settings.max_meshlet_export
                                                                                      << ", fix_missing_meshes="
                                                                                      << settings.fix_missing_meshes
                                                                                      << ", store_mesh_metadata="
                                                                                      << mesh_builder->store_mesh_metadata
                                                                                      << ", diagnostics="
                                                                                      << settings.diagnostics
                                                                                      << ", mesh_cap_at_start="
                                                                                      << settings.mesh_cap_at_start
                                                                                      << ", mesh_cap_at_end="
                                                                                      << settings.mesh_cap_at_end
                                                                                      << ", check_sites_inside_convex_hull="
                                                                                      << settings.check_sites_inside_convex_hull
                                                                                      << ", radius_vertex_shift_enabled="
                                                                                      << mesh_builder
                                                                                           ->radius_boundary_transition_shift_enabled);
  kinetic_delaunay->init(mesh_builder.get());
  kinetic_delaunay->compute();
  KINDS_INFO("Kinetic Delaunay Voronoi Meshing finished.");

  std::tie(segment_meshlets, meshing_neighbor_indices) = mesh_builder->extractSegmentMeshlets(true);
  segment_meshlet_export_suffixes = mesh_builder->extractSegmentMeshletExportSuffixes(true);
  size_t flipped_meshlet_count = 0;
  for (VoronoiMesh& meshlet : segment_meshlets)
  {
    if (meshlet.orientFacesAwayFromCentroid())
    {
      ++flipped_meshlet_count;
    }
  }
  KINDS_INFO(flipped_meshlet_count << "/" << segment_meshlets.size()
                                   << " meshlets flipped to orient faces away from centroid.");
}

void TreeMesher::mapMeshingToPhysicsSegmentIndices()
{
  const auto& meshing_strand_to_segment_indices = mesh_builder->getStrandToSegmentIndices();

  size_t max_meshing_id = 0;
  for (const auto& meshing_strand_to_segment_indice : meshing_strand_to_segment_indices)
  {
    for (size_t segment_no = 0; segment_no < meshing_strand_to_segment_indice.size(); ++segment_no)
    {
      size_t meshing_segment_id = meshing_strand_to_segment_indice[segment_no];
      max_meshing_id = std::max(max_meshing_id, meshing_segment_id);
    }
  }

  meshing_to_physics_segment_indices.clear();
  meshing_to_physics_segment_indices.resize(max_meshing_id + 1, -1);
  for (size_t strand_id = 0; strand_id < strand_tree.getPhysicsStrandToSegmentIndices().size(); ++strand_id)
  {
    for (size_t segment_no = 0; segment_no < meshing_strand_to_segment_indices[strand_id].size(); ++segment_no)
    {
      size_t meshing_segment_id = meshing_strand_to_segment_indices[strand_id][segment_no];
      int physics_segment_id = strand_tree.getPhysicsStrandToSegmentIndices()[strand_id][segment_no];
      meshing_to_physics_segment_indices[meshing_segment_id] = physics_segment_id;
    }
  }
}

const std::vector<VoronoiMesh>& kinDS::TreeMesher::runMeshingAlgorithm(bool visual_debug)
{
  runKineticDelaunay(visual_debug);

  auto& boundary_mesh = getBoundaryMesh();

  if (settings.debug_export_meshes)
  {
    kinDS::ObjExporter::writeMesh(
      boundary_mesh, "boundary_mesh.obj", 1.0, 1.0, {}, mesh_builder ? mesh_builder->store_mesh_metadata : settings.store_mesh_metadata);
  }

  // truncateToBoundary(boundary_mesh);

  mapMeshingToPhysicsSegmentIndices();

  if (mesh_builder && mesh_builder->store_mesh_metadata && settings.validate_mesh_vertex_sources)
  {
    mesh_vertex_source_validation_passed_ = Validator::validateAndReportMeshVertexSources(
      segment_meshlets, "segment meshlets", &segment_meshlet_export_suffixes);
  }
  else
  {
    mesh_vertex_source_validation_passed_ = true;
  }

  return segment_meshlets;
}

const VoronoiMesh& TreeMesher::getBoundaryMesh() const
{
  if (mesh_builder)
  {
    return mesh_builder->getBoundaryMesh();
  }

  throw std::runtime_error("Boundary mesh is not available before running the meshing algorithm.");
}

const std::vector<std::vector<int>>& kinDS::TreeMesher::getMeshingNeighborIndices() const
{
  return meshing_neighbor_indices;
}

const std::vector<size_t>& kinDS::TreeMesher::getMeshingToPhysicsSegmentIndices() const
{
  return meshing_to_physics_segment_indices;
}

const std::vector<std::vector<size_t>>& kinDS::TreeMesher::getMeshingStrandToSegmentIndices() const
{
  if (mesh_builder)
  {
    return mesh_builder->getStrandToSegmentIndices();
  }

  throw std::runtime_error("Strand to segment indices are not available before running the meshing algorithm.");
}
const std::vector<size_t>& TreeMesher::getBoundaryVertexToStrandId() const
{
  return mesh_builder->getBoundaryVertexToStrandId();
}

static glm::dvec3 ProfileToModelCoordinates(const std::vector<std::vector<glm::dmat4>>& profile_to_model_transforms,
  glm::dvec3 point, double t, const std::vector<size_t>& branch_indices, double w = 1.0f)
{
  size_t lower_section_index = static_cast<size_t>(std::max(0.0, glm::floor(t)));

  size_t upper_section_index = std::min(profile_to_model_transforms.size() - 1, static_cast<size_t>(glm::ceil(t)));

  // check range
  auto coord_str = std::to_string(t);
  if (lower_section_index >= profile_to_model_transforms.size())
  {
    std::cout << ("ProfileToModelCoordinates: lower bound of point z-coordinate out of range: " + coord_str).c_str()
              << std::endl;
  }
  if (upper_section_index >= profile_to_model_transforms.size())
  {
    std::cout << ("ProfileToModelCoordinates: upper bound of point z-coordinate out of range: " + coord_str).c_str()
              << std::endl;
  }

  // only set second coordinate to 0 for points, not for normal vectors
  // TODO: I actually wanted to get rid of this coordinate swap at some point
  glm::dvec4 local_pos(point[0], (1.0f - w) * point[2], point[1], w);
  size_t lower_branch_index = branch_indices[lower_section_index];
  glm::dvec4 global_pos = profile_to_model_transforms[lower_section_index][lower_branch_index] * local_pos;

  if (upper_section_index != lower_section_index)
  {
    size_t upper_branch_index = branch_indices[upper_section_index];
    glm::dvec4 upper_global_pos = profile_to_model_transforms[upper_section_index][upper_branch_index] * local_pos;
    float frac = static_cast<float>(t - static_cast<double>(lower_section_index));
    global_pos = glm::mix(global_pos, upper_global_pos, frac);
  }

  if (w == 0.0f)
  {
    global_pos = glm::normalize(global_pos);
  }

  return glm::dvec3(global_pos);
}

static const std::vector<size_t>& transformBranchIndicesForStrand(const StrandTree& strand_tree, size_t strand_id)
{
  return strand_tree.getBranchIndices(strand_id);
}

void TreeMesher::transformBoundaryMesh(kinDS::VoronoiMesh& boundary_mesh, const glm::dmat4& root_transform)
{
  auto& vertices = boundary_mesh.getVertices();
  std::vector<double> kinetic_times(vertices.size(), 0.0);
  for (size_t i = 0; i < vertices.size(); i++)
  {
    size_t strand_id = getBoundaryVertexToStrandId()[i];
    auto& v = vertices[i];
    kinetic_times[i] = v[2];
    const std::vector<size_t>& branch_indices = transformBranchIndicesForStrand(strand_tree, strand_id);
    // v is a relative position in 2D, we need to convert it to 3D
    v = glm::dvec3(root_transform
      * glm::dvec4(
        ProfileToModelCoordinates(strand_tree.getTransformsByHeightAndBranch(), v, kinetic_times[i], branch_indices),
        1.0));
  }

  const auto& triangles = boundary_mesh.getTriangles();
  for (size_t triangle_vertex_index = 0; triangle_vertex_index < triangles.size(); triangle_vertex_index += 3)
  {
    size_t material_id = boundary_mesh.getMaterialIDs()[triangle_vertex_index / 3];

    for (size_t j = 0; j < 3; j++)
    {
      auto source_tri_vertex_index = boundary_mesh.getTriangles()[triangle_vertex_index + j];

      size_t strand_id = getBoundaryVertexToStrandId()[source_tri_vertex_index];
      size_t normal_index = triangle_vertex_index + j;
      const glm::vec3& old_normal = boundary_mesh.getNormal(normal_index);
      const double kinetic_time = kinetic_times[source_tri_vertex_index];
      const std::vector<size_t>& branch_indices = transformBranchIndicesForStrand(strand_tree, strand_id);
      glm::vec3 transformed_normal = root_transform
        * glm::vec4(ProfileToModelCoordinates(strand_tree.getNormalTransformsByHeightAndBranch(), old_normal,
                      kinetic_time, branch_indices, 0.0f),
          0.0f);
      boundary_mesh.replaceNormal(normal_index, transformed_normal);
    }
  }
}

void TreeMesher::transformToWorldSpace(VoronoiMesh& mesh, size_t strand_id, const glm::dmat4& root_transform) const
{
  auto& vertices = mesh.getVertices();
  std::vector<double> kinetic_times(vertices.size(), 0.0);

  // transform points
  for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index)
  {
    auto& v = vertices[vertex_index];
    kinetic_times[vertex_index] = v.z;
    const std::vector<size_t>& branch_indices = transformBranchIndicesForStrand(strand_tree, strand_id);
    v = root_transform
      * glm::dvec4(ProfileToModelCoordinates(
                     strand_tree.getTransformsByHeightAndBranch(), v, kinetic_times[vertex_index], branch_indices),
        1.0);
  }

  if (mesh.getNormalMode() == NormalMode::NoNormals)
  {
    return;
  }

  auto& normal_transforms = strand_tree.getNormalTransformsByHeightAndBranch();

  if (mesh.getNormalMode() == NormalMode::PerTriangleCorner)
  {
    const auto& triangles = mesh.getTriangles();
    for (size_t triangle_vertex_index = 0; triangle_vertex_index < triangles.size(); ++triangle_vertex_index)
    {
      const size_t vertex_index = triangles[triangle_vertex_index];
      const glm::dvec3& old_normal = mesh.getNormal(triangle_vertex_index);
      const double kinetic_time = kinetic_times[vertex_index];
      const std::vector<size_t>& branch_indices = transformBranchIndicesForStrand(strand_tree, strand_id);
      const glm::dvec3 transformed_normal = root_transform
        * glm::dvec4(ProfileToModelCoordinates(normal_transforms, old_normal, kinetic_time, branch_indices, 0.0f),
          0.0);
      mesh.setNormal(transformed_normal, triangle_vertex_index);
    }
  }
  else
  {
    auto& normals = mesh.getNormals();
    for (size_t vertex_index = 0; vertex_index < normals.size(); ++vertex_index)
    {
      glm::dvec3& n = normals[vertex_index];
      const double kinetic_time = kinetic_times[vertex_index];
      const std::vector<size_t>& branch_indices = transformBranchIndicesForStrand(strand_tree, strand_id);
      n = root_transform
        * glm::dvec4(ProfileToModelCoordinates(normal_transforms, n, kinetic_time, branch_indices, 0.0f),
          0.0);
    }
  }
}
