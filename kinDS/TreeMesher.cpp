#include "TreeMesher.hpp"

#include "IndexIterator.hpp"
#include "KineticDelaunay.hpp"
#include "MeshIntersection.hpp"
#include "ObjExporter.hpp"
#include "SegmentBuilder.hpp"

#include <execution>
#include <set>

using namespace kinDS;

std::vector<std::pair<size_t, double>> MergeSortedVectors(const std::vector<std::vector<double>>& inputs)
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

void TreeMesher::RunMeshingAlgorithm(const std::vector<std::vector<glm::dvec2>>& support_points,
  std::vector<std::vector<double>>& subdivisions_by_strand,
  std::vector<std::vector<int>>& physics_strand_to_segment_indices,
  const std::vector<std::vector<glm::mat4>>& transforms_by_height_and_branch,
  const std::vector<std::vector<size_t>>& branch_indices,
  std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id,
  std::vector<float>& boundary_distance_by_segment_id, double alpha_cutoff)
{
  bool recompute_segment_pairs = false; // TODO: expose as option?

  std::vector<std::vector<glm::mat4>> normal_transforms_by_height_and_branch(transforms_by_height_and_branch.size());

  for (size_t i = 0; i < transforms_by_height_and_branch.size(); i++)
  {
    normal_transforms_by_height_and_branch[i].resize(transforms_by_height_and_branch[i].size());
    for (size_t j = 0; j < normal_transforms_by_height_and_branch[i].size(); j++)
    {
      normal_transforms_by_height_and_branch[i][j]
        = glm::transpose(glm::inverse(transforms_by_height_and_branch[i][j]));
      normal_transforms_by_height_and_branch[i][j][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }
  }

  std::vector<float> bottom_boundary_distances_by_strand_id(physics_strand_to_segment_indices.size());
  std::vector<float> top_boundary_distances_by_strand_id(physics_strand_to_segment_indices.size());

  for (size_t strand_id = 0; strand_id < physics_strand_to_segment_indices.size(); strand_id++)
  {
    int bottom_segment_id = physics_strand_to_segment_indices[strand_id].front();
    int top_segment_id = physics_strand_to_segment_indices[strand_id].back();

    bottom_boundary_distances_by_strand_id[strand_id] = boundary_distance_by_segment_id[bottom_segment_id];
    top_boundary_distances_by_strand_id[strand_id] = boundary_distance_by_segment_id[top_segment_id];
  }

  // sort subdivisions into a single array
  std::vector<std::pair<size_t, double>> subdivisions = MergeSortedVectors(subdivisions_by_strand);

  KINDS_INFO("Starting Kinetic Delaunay Voronoi Meshing...");
  kinDS::KineticDelaunay kinetic_delaunay(
    kinDS::BranchTrajectories(support_points, transforms_by_height_and_branch, branch_indices, strands_by_branch_id),
    alpha_cutoff, false, branch_indices, strands_by_branch_id);

  bool transform_mesh_at_construction = false;

  kinetic_delaunay.init();
  kinDS::SegmentBuilder mesh_builder(kinetic_delaunay, subdivisions, transform_mesh_at_construction);
  mesh_builder.init();
  // auto points = kinetic_delaunay.getPointsAt(0.0);

  size_t section_count = kinetic_delaunay.getSectionCount();

  ProgressBar section_progress_bar(
    0, section_count, "Computing Kinetic Voronoi Sections", ProgressBar::Display::Absolute);
  for (size_t i = 0; i < section_count; ++i)
  {
    section_progress_bar.Update(i);
    if (i != 0)
      mesh_builder.betweenSections(i);
    kinetic_delaunay.advanceOneSection(mesh_builder);

    // points = kinetic_delaunay.getPointsAt(static_cast<double>(i + 1));
  }
  section_progress_bar.Finish();

  KINDS_INFO("Finalizing Kinetic Delaunay Voronoi Meshing...");
  mesh_builder.finalize(section_count);

  auto [meshes, meshing_neighbor_indices] = mesh_builder.extractSegmentMeshlets();

  auto& boundary_mesh = mesh_builder.getBoundaryMesh();
  auto& boundary_vertex_to_strand_id = mesh_builder.getBoundaryVertexToStrandId();

  bool debug_export_meshes = false;

  if (debug_export_meshes)
  {
    kinDS::ObjExporter::writeMesh(boundary_mesh, "boundary_mesh.obj");
  }
  size_t max_meshlet_export = 50000000;
  // intersect all meshes with the boundary mesh and save the result
  // Build an AABB-tree of the boundary-mesh to prefilter
  kinDS::MeshIntersection boundary_intersector(boundary_mesh);

  // ProgressBar intersection_progress_bar(0, meshes.size(), "Computing Mesh Intersections",
  //                                       ProgressBar::Display::Absolute, 50);

  std::atomic<int> progress_counter { 0 };
  std::for_each(std::execution::par, IndexIterator { 0 }, IndexIterator { meshes.size() },
    [&](auto mesh_index)
    {
      progress_counter.fetch_add(1, std::memory_order_relaxed);
      // intersection_progress_bar.Update(progress_counter);
      auto intersect_relation = boundary_intersector.ClassifyMeshRelation(meshes[mesh_index], true);

      switch (intersect_relation)
      {
      case kinDS::MeshIntersection::MeshRelation::INSIDE:
        // do nothing
        break;

      case kinDS::MeshIntersection::MeshRelation::INTERSECTING:
        if (debug_export_meshes && mesh_index < max_meshlet_export)
        {
          kinDS::ObjExporter::writeMesh(meshes[mesh_index], "meshlet" + std::to_string(mesh_index) + "_raw.obj");
        }
        std::tie(meshes[mesh_index], meshing_neighbor_indices[mesh_index])
          = boundary_intersector.Intersect(meshes[mesh_index], meshing_neighbor_indices[mesh_index]);
        break;

      case kinDS::MeshIntersection::MeshRelation::OUTSIDE:
        // fully outside, result is empty mesh
        meshes[mesh_index] = kinDS::VoronoiMesh();
        meshing_neighbor_indices[mesh_index] = {};
        break;

      case kinDS::MeshIntersection::MeshRelation::UNDEFINED:
        KINDS_ERROR("Mesh relation returned UNDEFINED");
        break;

      default:
        KINDS_ERROR("Unknown return value of mesh relation");
        break;
      }
    });

  // intersection_progress_bar.Finish();

  // Find empty meshes and try to fix them by expanding neighboring meshes
  bool fix_missing_meshes = false;
  if (fix_missing_meshes)
  {
    std::vector<size_t> empty_mesh_indices;
    for (size_t mesh_index = 0; mesh_index < meshes.size(); mesh_index++)
    {
      if (meshes[mesh_index].getVertexCount() == 0)
      {
        empty_mesh_indices.push_back(mesh_index);
      }
    }

    // go through all triangles of all meshes and check their neighbor indices
    // if a neighbor index corresponds to an empty mesh, try to copy the triangle into the corresponding empty mesh
    for (size_t mesh_index = 0; mesh_index < meshes.size(); mesh_index++)
    {
      if (std::binary_search(empty_mesh_indices.begin(), empty_mesh_indices.end(), mesh_index))
      {
        continue;
      }

      auto& mesh = meshes[mesh_index];
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
            kinDS::VoronoiMesh& neighbor_mesh = meshes[neighbor_mesh_index];

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
            neighbor_mesh.addTriangle(
              new_v1_index, new_v0_index, new_v2_index, new_t0_index, new_t1_index, new_t2_index);
            meshing_neighbor_indices[neighbor_mesh_index].push_back(mesh_index);
          }
        }
      }
    }

    // merge vertices and check how many could be (partially) fixed
    size_t fixed_mesh_count = 0;
    for (size_t mesh_index : empty_mesh_indices)
    {
      auto& mesh = meshes[mesh_index];

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

    KINDS_INFO("Fixed " << fixed_mesh_count << " out of " << empty_mesh_indices.size() << " meshes.");
  }

  const auto& meshing_strand_to_segment_indices = mesh_builder.getStrandToSegmentIndices();

  size_t max_meshing_id = 0;
  for (size_t strand_id = 0; strand_id < meshing_strand_to_segment_indices.size(); ++strand_id)
  {
    for (size_t segment_no = 0; segment_no < meshing_strand_to_segment_indices[strand_id].size(); ++segment_no)
    {
      size_t meshing_segment_id = meshing_strand_to_segment_indices[strand_id][segment_no];
      max_meshing_id = std::max(max_meshing_id, meshing_segment_id);
    }
  }

  std::vector<size_t> meshing_to_physics_segment_indices(max_meshing_id + 1, -1);
  for (size_t strand_id = 0; strand_id < physics_strand_to_segment_indices.size(); ++strand_id)
  {
    for (size_t segment_no = 0; segment_no < meshing_strand_to_segment_indices[strand_id].size(); ++segment_no)
    {
      size_t meshing_segment_id = meshing_strand_to_segment_indices[strand_id][segment_no];
      int physics_segment_id = physics_strand_to_segment_indices[strand_id][segment_no];
      meshing_to_physics_segment_indices[meshing_segment_id] = physics_segment_id;
    }
  }

  // for now, just combine all meshes into one
  kinDS::VoronoiMesh combined_mesh;
  for (const auto& mesh : meshes)
  {
    combined_mesh += mesh;
  }

  // also export some meshlets:
  for (size_t i = 0; i < std::min(max_meshlet_export, meshes.size()); i++)
  {
    kinDS::ObjExporter::writeMesh(meshes[i], "meshlet" + std::to_string(i) + ".obj");
  }

  // export combined mesh
  // combined_mesh.mergeDuplicateVertices(0.0001);
  // kinDS::ObjExporter::writeMesh(transformed_mesh, "transformed_mesh.obj");
  // kinDS::ObjExporter::writeMesh(transformed_boundary_mesh, "transformed_boundary_mesh.obj");
  // kinDS::ObjExporter::writeMesh(combined_mesh, "meshtest_subdivided.obj");
  KINDS_INFO("Kinetic Delaunay Voronoi Meshes exported.");
}