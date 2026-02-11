#include "TreeMesher.hpp"
#include "IndexIterator.hpp"
#include "KineticDelaunay.hpp"
#include "ObjExporter.hpp"
#include "SegmentBuilder.hpp"
#include <execution>

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

void TreeMesher::exportCombinedMesh() const
{
  // for now, just combine all meshes into one
  kinDS::VoronoiMesh combined_mesh;
  for (const auto& mesh : segment_meshlets)
  {
    combined_mesh += mesh;
  }

  // also export some meshlets:
  for (size_t i = 0; i < std::min(settings.max_meshlet_export, segment_meshlets.size()); i++)
  {
    kinDS::ObjExporter::writeMesh(segment_meshlets[i], "meshlet" + std::to_string(i) + ".obj");
  }
  KINDS_INFO("Kinetic Delaunay Voronoi Meshes exported.");
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

  parallel_for(segment_meshlets.size(),
    [&](auto mesh_index)
    {
      // progress_counter.fetch_add(1, std::memory_order_relaxed);
      // intersection_progress_bar.Update(progress_counter);
      auto intersect_relation = boundary_intersector.ClassifyMeshRelation(segment_meshlets[mesh_index], true);

      switch (intersect_relation)
      {
      case kinDS::MeshIntersection::MeshRelation::INSIDE:
        // do nothing
        break;

      case kinDS::MeshIntersection::MeshRelation::INTERSECTING:
        if (settings.debug_export_meshes && mesh_index < settings.max_meshlet_export)
        {
          kinDS::ObjExporter::writeMesh(
            segment_meshlets[mesh_index], "meshlet" + std::to_string(mesh_index) + "_raw.obj");
        }
        std::tie(segment_meshlets[mesh_index], meshing_neighbor_indices[mesh_index])
          = boundary_intersector.Intersect(segment_meshlets[mesh_index], meshing_neighbor_indices[mesh_index]);
        break;

      case kinDS::MeshIntersection::MeshRelation::OUTSIDE:
        // fully outside, result is empty mesh
        segment_meshlets[mesh_index] = kinDS::VoronoiMesh();
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
          neighbor_mesh.addTriangle(new_v1_index, new_v0_index, new_v2_index, new_t0_index, new_t1_index, new_t2_index);
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

  KINDS_INFO("Fixed " << fixed_mesh_count << " out of " << empty_mesh_indices.size() << " meshes.");
}

std::vector<std::vector<glm::mat4>> TreeMesher::computeNormalTransforms()
{
  std::vector<std::vector<glm::mat4>> normal_transforms_by_height_and_branch(
    strand_tree.getTransformsByHeightAndBranch().size());

  for (size_t i = 0; i < strand_tree.getTransformsByHeightAndBranch().size(); i++)
  {
    normal_transforms_by_height_and_branch[i].resize(strand_tree.getTransformsByHeightAndBranch()[i].size());
    for (size_t j = 0; j < normal_transforms_by_height_and_branch[i].size(); j++)
    {
      normal_transforms_by_height_and_branch[i][j]
        = glm::transpose(glm::inverse(strand_tree.getTransformsByHeightAndBranch()[i][j]));
      normal_transforms_by_height_and_branch[i][j][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }
  }

  return normal_transforms_by_height_and_branch;
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

void TreeMesher::runKineticDelaunay()
{
  KINDS_INFO("Starting Kinetic Delaunay Voronoi Meshing...");
  // sort subdivisions into a single array
  std::vector<std::pair<size_t, double>> subdivisions = MergeSortedVectors(strand_tree.getSubdivisionsByStrand());

  kinetic_delaunay = std::make_shared<KineticDelaunay>(strand_tree, settings.alpha_cutoff, false);

  bool transform_mesh_at_construction = false;

  kinetic_delaunay->init();
  mesh_builder = std::make_shared<SegmentBuilder>(*kinetic_delaunay, subdivisions, transform_mesh_at_construction);
  mesh_builder->init();
  // auto points = kinetic_delaunay.getPointsAt(0.0);

  size_t section_count = kinetic_delaunay->getSectionCount();

  ProgressBar section_progress_bar(
    0, section_count, "Computing Kinetic Voronoi Sections", ProgressBar::Display::Absolute);
  for (size_t i = 0; i < section_count; ++i)
  {
    section_progress_bar.Update(i);
    if (i != 0)
      mesh_builder->betweenSections(i);
    kinetic_delaunay->advanceOneSection(*mesh_builder);

    // points = kinetic_delaunay.getPointsAt(static_cast<double>(i + 1));
  }
  section_progress_bar.Finish();

  KINDS_INFO("Finalizing Kinetic Delaunay Voronoi Meshing...");
  mesh_builder->finalize(section_count);

  std::tie(segment_meshlets, meshing_neighbor_indices) = mesh_builder->extractSegmentMeshlets();
}

void TreeMesher::mapMeshingToPhysicsSegmentIndices()
{
  const auto& meshing_strand_to_segment_indices = mesh_builder->getStrandToSegmentIndices();

  size_t max_meshing_id = 0;
  for (size_t strand_id = 0; strand_id < meshing_strand_to_segment_indices.size(); ++strand_id)
  {
    for (size_t segment_no = 0; segment_no < meshing_strand_to_segment_indices[strand_id].size(); ++segment_no)
    {
      size_t meshing_segment_id = meshing_strand_to_segment_indices[strand_id][segment_no];
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

const std::vector<VoronoiMesh>& kinDS::TreeMesher::runMeshingAlgorithm()
{
  runKineticDelaunay();

  auto& boundary_mesh = getBoundaryMesh();

  if (settings.debug_export_meshes)
  {
    kinDS::ObjExporter::writeMesh(boundary_mesh, "boundary_mesh.obj");
  }

  truncateToBoundary(boundary_mesh);

  mapMeshingToPhysicsSegmentIndices();

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
