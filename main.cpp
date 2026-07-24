#include "kinDS/HalfEdgeDelaunayGraphToSVG.hpp"
#include "kinDS/KineticDelaunay.hpp"
#include "kinDS/Logger.hpp"
#include "kinDS/ObjExporter.hpp"
#include "kinDS/Polynomial.hpp"
#include "kinDS/SegmentBuilder.hpp"
#include "kinDS/StrandTree.hpp"
#include "kinDS/TreeMesher.hpp"
#include "kinDS/Validator.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <utility> // for std::pair
#include <vector>

std::vector<std::pair<size_t, double>> merge_sorted_vectors(const std::vector<std::vector<double>>& inputs)
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

static bool kinetic_delaunay_example(bool validate_mesh_vertex_sources, const std::string& validate_log_path,
  bool alternate_section_shading)
{
#define TEST_TRAJECTORIES
#ifdef TEST_TRAJECTORIES
  std::vector<glm::dvec2> trajectory_A = {
    //{ -0.420113, -0.558875 },
    { -0.432132, -0.426942 },
    { -0.447292, -0.580708 },
    { -0.469864, -0.531837 },
    { -0.578741, -0.494280 },
    { -0.519044, -0.496727 },
    { -0.487418, -0.587100 },
    { -0.536664, -0.465019 },
  };

  std::vector<glm::dvec2> trajectory_B = {
    //{ -0.200000, -0.500000 },
    { -0.150887, -0.424968 },
    { -0.101774, -0.349936 },
    { -0.052661, -0.274904 },
    { -0.003548, -0.199872 },
    { 0.045565, -0.124840 },
    { 0.094678, -0.049808 },
    { 0.143791, 0.025224 },
  };

  std::vector<glm::dvec2> trajectory_C = {
    //{ 0.100000, -0.400000 },
    { -0.048665, -0.333097 },
    { -0.197330, -0.266194 },
    { -0.345995, -0.199291 },
    { -0.494660, -0.132388 },
    { -0.643325, -0.065485 },
    { -0.791990, 0.001418 },
    { -0.940656, 0.068321 },
  };

  std::vector<glm::dvec2> trajectory_D = {
    //{ 0.500000, 0.200000 },
    { 0.467745, 0.111272 },
    { 0.435490, 0.022544 },
    { 0.403235, -0.066183 },
    { 0.370980, -0.154910 },
    { 0.338725, -0.243637 },
    { 0.306470, -0.332364 },
    { 0.274215, -0.421091 },
  };

  /* std::vector<glm::dvec2> trajectory_B = {
    { 0.556188, -0.476838 },
    { 0.572065, -0.428676 },
    { 0.478597, -0.594223 },
    { 0.520113, -0.439088 },
    { 0.443697, -0.404739 },
    { 0.490520, -0.404842 },
    { 0.457843, -0.405508 },
    { 0.430738, -0.451790 },
  };

  std::vector<kinDS::Point<2>> trajectory_C = {
    { 0.541622, 0.548149 },
    { 0.429803, 0.496248 },
    { 0.413522, 0.498909 },
    { 0.580101, 0.499790 },
    { 0.479949, 0.562262 },
    { 0.593821, 0.428951 },
    { 0.595407, 0.407437 },
    { 0.447700, 0.444491 },
  };

  std::vector<glm::dvec2> trajectory_D = {
    { -0.559656, 0.421442 },
    { -0.464952, 0.428177 },
    { -0.578349, 0.568869 },
    { -0.453764, 0.424407 },
    { -0.420543, 0.594949 },
    { -0.489771, 0.566934 },
    { -0.431853, 0.415965 },
    { -0.428599, 0.554348 },
  };*/

  std::vector<std::vector<glm::dvec2>> support_points { trajectory_A, trajectory_B, trajectory_C, trajectory_D };

#else
  std::vector<std::vector<glm::dvec2>> strand_guide_points = {
    { kinDS::Point<2> { 4.12703, -1.57022 }, kinDS::Point<2> { 4.12703, -1.57022 },
      kinDS::Point<2> { 4.12703, -1.57022 } },
    { kinDS::Point<2> { 4.3178, -3.52816 }, kinDS::Point<2> { 4.3178, -3.52816 },
      kinDS::Point<2> { 4.3178, -3.52816 } },
    { kinDS::Point<2> { -2.97276, 4.03951 }, kinDS::Point<2> { -2.97276, 4.03951 },
      kinDS::Point<2> { -2.97276, 4.03951 } },
    { kinDS::Point<2> { 3.88387, 2.44586 }, kinDS::Point<2> { 3.88387, 2.44586 },
      kinDS::Point<2> { 3.88387, 2.44586 } },
    { kinDS::Point<2> { -2.87601, 0.0555282 }, kinDS::Point<2> { -2.87601, 0.0555282 },
      kinDS::Point<2> { -2.87601, 0.0555282 } },
    { kinDS::Point<2> { -4.58249, 0.899374 }, kinDS::Point<2> { -4.58249, 0.899374 },
      kinDS::Point<2> { -4.58249, 0.899374 } },
    { kinDS::Point<2> { 2.64728, -2.65377 }, kinDS::Point<2> { 2.64728, -2.65377 },
      kinDS::Point<2> { 2.64728, -2.65377 } },
    { kinDS::Point<2> { -0.927691, -0.936193 }, kinDS::Point<2> { -0.927691, -0.936193 },
      kinDS::Point<2> { -0.927691, -0.936193 } },
    { kinDS::Point<2> { 4.04879, 0.409984 }, kinDS::Point<2> { 4.04879, 0.409984 },
      kinDS::Point<2> { 4.04879, 0.409984 } },
    { kinDS::Point<2> { 2.49795, -0.656439 }, kinDS::Point<2> { 2.49795, -0.656439 },
      kinDS::Point<2> { 2.49795, -0.656439 } },
    { kinDS::Point<2> { 0.861075, -3.77961 }, kinDS::Point<2> { 0.861075, -3.77961 },
      kinDS::Point<2> { 0.861075, -3.77961 } },
    { kinDS::Point<2> { -2.54759, -3.88424 }, kinDS::Point<2> { -2.54759, -3.88424 },
      kinDS::Point<2> { -2.54759, -3.88423 } },
    { kinDS::Point<2> { -4.39657, -0.990622 }, kinDS::Point<2> { -4.39657, -0.990622 },
      kinDS::Point<2> { -4.39657, -0.990622 } },
    { kinDS::Point<2> { -0.764416, -4.80726 }, kinDS::Point<2> { -0.764416, -4.80726 },
      kinDS::Point<2> { -0.764416, -4.80726 } },
    { kinDS::Point<2> { -1.018, 1.12591 }, kinDS::Point<2> { -1.018, 1.12591 }, kinDS::Point<2> { -1.018, 1.12591 } },
    { kinDS::Point<2> { 0.724303, -1.71089 }, kinDS::Point<2> { 0.724303, -1.71089 },
      kinDS::Point<2> { 0.724303, -1.71089 } },
    { kinDS::Point<2> { 0.429867, 2.15062 }, kinDS::Point<2> { 0.429867, 2.15062 },
      kinDS::Point<2> { 0.429867, 2.15062 } },
    { kinDS::Point<2> { 2.19038, 3.23518 }, kinDS::Point<2> { 2.19038, 3.23517 },
      kinDS::Point<2> { 2.19038, 3.23517 } },
    { kinDS::Point<2> { -1.13884, 3.0388 }, kinDS::Point<2> { -1.13884, 3.0388 },
      kinDS::Point<2> { -1.13884, 3.0388 } },
    { kinDS::Point<2> { -4.52197, 2.9438 }, kinDS::Point<2> { -4.52197, 2.9438 },
      kinDS::Point<2> { -4.52197, 2.9438 } },
    { kinDS::Point<2> { 5.8964, -0.557665 }, kinDS::Point<2> { 5.8964, -0.557665 },
      kinDS::Point<2> { 5.8964, -0.557665 } },
    { kinDS::Point<2> { 2.44269, 1.24945 }, kinDS::Point<2> { 2.44269, 1.24945 },
      kinDS::Point<2> { 2.44269, 1.24945 } },
    { kinDS::Point<2> { 2.03784, 5.19061 }, kinDS::Point<2> { 2.03784, 5.19061 },
      kinDS::Point<2> { 2.03784, 5.19061 } },
    { kinDS::Point<2> { -0.901509, -2.7655 }, kinDS::Point<2> { -0.901509, -2.7655 },
      kinDS::Point<2> { -0.901509, -2.7655 } },
    { kinDS::Point<2> { 0.289699, 4.10103 }, kinDS::Point<2> { 0.289699, 4.10103 },
      kinDS::Point<2> { 0.289699, 4.10103 } },
    { kinDS::Point<2> { 0.59763, 0.123849 }, kinDS::Point<2> { 0.59763, 0.123849 },
      kinDS::Point<2> { 0.59763, 0.123849 } },
    { kinDS::Point<2> { -4.27136, -3.04379 }, kinDS::Point<2> { -4.27136, -3.04379 },
      kinDS::Point<2> { -4.27136, -3.04379 } },
    { kinDS::Point<2> { 5.76047, 1.54084 }, kinDS::Point<2> { 5.76047, 1.54084 },
      kinDS::Point<2> { 5.76047, 1.54084 } },
    { kinDS::Point<2> { -2.88041, 1.87929 }, kinDS::Point<2> { -2.88041, 1.87929 },
      kinDS::Point<2> { -2.88041, 1.87929 } },
    { kinDS::Point<2> { -2.59, -2.00589 }, kinDS::Point<2> { -2.59, -2.00589 }, kinDS::Point<2> { -2.59, -2.00589 } }
  };

  std::vector<kinDS::CubicHermiteSpline<2>> splines;

  for (auto& points : strand_guide_points)
  {
    splines.emplace_back(points);
  }
#endif

  std::vector<std::vector<double>> subdivisions = { { 0.4 }, { 0.8 }, { 0.2 }, { 0.6 } };

  // Sort subdivisions into pairs of strand index and parameter
  std::vector<std::pair<size_t, double>> sorted_subdivisions = merge_sorted_vectors(subdivisions);

  // Create a branch index lookup using branch_indices[strand_id][h] = branch_id
  const std::vector<std::vector<size_t>> branch_indices(
    support_points.size(), std::vector<size_t>(support_points.front().size(), 0));

  // Maintain the branches as strands_by_branch_id[h][branch_id][strand_no] = strand_id
  std::vector<std::vector<std::vector<size_t>>> strands_by_branch_id;

  // build strands_by_branch_id from branch_indices
  for (size_t h = 0; h < support_points.front().size(); ++h)
  {
    std::map<size_t, std::vector<size_t>> branch_to_strands;
    for (size_t strand_id = 0; strand_id < support_points.size(); ++strand_id)
    {
      size_t branch_id = branch_indices[strand_id][h];
      branch_to_strands[branch_id].push_back(strand_id);
    }
    std::vector<std::vector<size_t>> strands_in_branches;
    for (auto& [branch_id, strands] : branch_to_strands)
    {
      strands_in_branches.push_back(strands);
    }
    strands_by_branch_id.push_back(strands_in_branches);
  }

  std::vector<std::vector<glm::dmat4>> transforms_by_height_and_branch;
  for (size_t h = 0; h < support_points.front().size(); ++h)
  {
    std::vector<glm::dmat4> transforms_at_height;
    size_t branch_count = strands_by_branch_id[h].size();
    for (size_t b = 0; b < branch_count; ++b)
    {
      transforms_at_height.push_back(glm::dmat4(1.0)); // identity
    }
    transforms_by_height_and_branch.push_back(transforms_at_height);
  }

  /*std::vector<std::vector<glm::dmat4>> normal_transforms_by_height_and_branch(transforms_by_height_and_branch.size());

  for (size_t i = 0; i < transforms_by_height_and_branch.size(); i++)
  {
    normal_transforms_by_height_and_branch[i].resize(transforms_by_height_and_branch[i].size());
    for (size_t j = 0; j < normal_transforms_by_height_and_branch[i].size(); j++)
    {
      normal_transforms_by_height_and_branch[i][j]
        = glm::transpose(glm::inverse(transforms_by_height_and_branch[i][j]));
      normal_transforms_by_height_and_branch[i][j][3] = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }
  }*/

  kinDS::KineticDelaunay kinetic_delaunay(kinDS::StrandTree(support_points, subdivisions, {},
                                            transforms_by_height_and_branch, branch_indices, strands_by_branch_id),
    10.0, false);

  kinDS::SegmentBuilder mesh_builder(kinetic_delaunay, sorted_subdivisions, false);
  mesh_builder.radius_boundary_transition_shift_enabled = true;
  mesh_builder.store_mesh_metadata = validate_mesh_vertex_sources;
  mesh_builder.validate_mesh_vertex_sources = validate_mesh_vertex_sources;
  kinetic_delaunay.init(&mesh_builder);
  auto points = kinetic_delaunay.getPointsAt(0.0);

  // Build Voronoi vertex -> containing triangle mapping from CrossingData for SVG labeling.
  const kinDS::HalfEdgeDelaunayGraph& demo_graph = kinetic_delaunay.getGraph();
  const size_t demo_face_count = demo_graph.faceSlotCount();
  std::vector<size_t> demo_voronoi_vertex_to_tri(demo_face_count);
  constexpr size_t demo_invalid_id = static_cast<size_t>(-1);
  for (size_t vid = 0; vid < demo_face_count; ++vid)
  {
    demo_voronoi_vertex_to_tri[vid] = kinetic_delaunay.getCrossingDataContainingTriId(vid);
  }

  // t=0 snapshot
  auto demo_intersection_debug_data = kinetic_delaunay.getCrossingIntersectionDebugData();
  kinDS::HalfEdgeDelaunayGraphToSVG::write(points, demo_graph, "t0.000000_demo.svg", 0.1, nullptr, true,
    &demo_voronoi_vertex_to_tri, &demo_intersection_debug_data);

  kinetic_delaunay.compute();

  // mesh_builder.printDebugInfo();
  bool merge_by_segment = true;
  // auto meshes = mesh_builder.extractMeshes();
  auto meshes = mesh_builder.extractSegmentMeshlets(merge_by_segment);
  auto mesh_suffixes = mesh_builder.extractSegmentMeshletExportSuffixes(merge_by_segment);
  //(0.1, 0.01, subdivisions);

  bool validation_passed = true;
  if (validate_mesh_vertex_sources)
  {
    kinDS::Validator::setLogFile(validate_log_path);
    validation_passed
      = kinDS::Validator::validateAndReportMeshVertexSources(meshes.first, "demo segment meshlets", &mesh_suffixes);
    std::cout << "Validation report written to: " << kinDS::Validator::logFilePath() << std::endl;
    if (!validation_passed)
    {
      std::cerr << "Mesh vertex source validation failed. See " << kinDS::Validator::logFilePath() << std::endl;
    }
  }

  for (size_t i = 0; i < meshes.first.size(); ++i)
  {
    const std::string suffix = (i < mesh_suffixes.size()) ? mesh_suffixes[i] : "";
    std::string filename = "mesh_" + std::to_string(i) + suffix + meshes.first[i].creationKineticTimeFilenameSuffix()
      + ".obj";
    const bool include_vertex_colors = kinDS::Validator::meshUsesValidationErrorMaterial(meshes.first[i]);
    kinDS::ObjExporter::writeMesh(meshes.first[i], filename, 1.0, 1.0, {}, true, include_vertex_colors,
      alternate_section_shading);
    std::cout << "Mesh saved to " << filename << std::endl;
  }

  auto boundary_interval_meshes = mesh_builder.extractBoundaryIntervalMeshlets();
  auto boundary_interval_suffixes = mesh_builder.extractBoundaryIntervalMeshletExportSuffixes();
  for (size_t i = 0; i < boundary_interval_meshes.size(); ++i)
  {
    const std::string suffix = (i < boundary_interval_suffixes.size()) ? boundary_interval_suffixes[i] : "";
    std::string filename = "mesh_bound_" + std::to_string(i) + suffix
      + boundary_interval_meshes[i].creationKineticTimeFilenameSuffix() + ".obj";
    kinDS::ObjExporter::writeMesh(boundary_interval_meshes[i], filename, 1.0, 1.0, {}, true, false,
      alternate_section_shading);
    std::cout << "Boundary interval mesh saved to " << filename << std::endl;
  }

  auto& boundary_mesh = mesh_builder.getBoundaryMesh();
  kinDS::ObjExporter::writeMesh(boundary_mesh, "boundary_mesh.obj", 1.0, 1.0, {}, true);

  boundary_mesh.checkForDegenerateTriangles();

  if (validate_mesh_vertex_sources)
  {
    return validation_passed;
  }
  return true;
}

// Helper: enable/disable specific log levels (relative to current)
static bool modify_log_level(const std::string& levels_str, bool enable)
{
  using namespace kinDS;

  // Parse comma-separated list
  std::stringstream ss(levels_str);
  std::string level;
  bool found_any = false;

  while (std::getline(ss, level, ','))
  {
    // Trim whitespace
    level.erase(0, level.find_first_not_of(" \t"));
    level.erase(level.find_last_not_of(" \t") + 1);

    // Convert to lowercase for case-insensitive matching
    std::string level_lower = level;
    std::transform(level_lower.begin(), level_lower.end(), level_lower.begin(), ::tolower);

    if (level_lower == "debug")
    {
      logger.setLogLevel(LogLevel::Debug, enable);
      found_any = true;
    }
    else if (level_lower == "monitor")
    {
      logger.setLogLevel(LogLevel::Monitor, enable);
      found_any = true;
    }
    else if (level_lower == "info")
    {
      logger.setLogLevel(LogLevel::Info, enable);
      found_any = true;
    }
    else if (level_lower == "warning" || level_lower == "warn")
    {
      logger.setLogLevel(LogLevel::Warning, enable);
      found_any = true;
    }
    else if (level_lower == "error")
    {
      logger.setLogLevel(LogLevel::Error, enable);
      found_any = true;
    }
    else if (level_lower == "critical")
    {
      logger.setLogLevel(LogLevel::Critical, enable);
      found_any = true;
    }
    else
    {
      std::cerr << "Warning: Unknown log level: " << level << std::endl;
    }
  }

  return found_any;
}

// Helper: describe currently enabled log levels as a comma-separated string
static std::string get_enabled_log_levels()
{
  using namespace kinDS;
  LogLevel mask = logger.getLogLevel();
  std::vector<std::string> names;

  auto add_if = [&](LogLevel lvl, const char* name)
  {
    if (static_cast<unsigned>(mask & lvl))
    {
      names.emplace_back(name);
    }
  };

  add_if(LogLevel::Debug, "debug");
  add_if(LogLevel::Monitor, "monitor");
  add_if(LogLevel::Info, "info");
  add_if(LogLevel::Warning, "warning");
  add_if(LogLevel::Error, "error");
  add_if(LogLevel::Critical, "critical");

  if (names.empty())
  {
    return "none";
  }

  std::string result;
  for (size_t i = 0; i < names.size(); ++i)
  {
    if (i > 0)
      result += ",";
    result += names[i];
  }
  return result;
}

// Absolute setter: replace current levels with exactly the ones given
static void set_log_level(const std::string& levels_str)
{
  using namespace kinDS;

  // First disable all known levels
  logger.setLogLevel(
    LogLevel::Debug | LogLevel::Monitor | LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical, false);

  bool ok = modify_log_level(levels_str, true);

  if (!ok)
  {
    std::cerr << "Warning: No valid log levels specified. Using default (info,warning,error,critical)." << std::endl;
    // Restore default: everything except Debug
    logger.setLogLevel(LogLevel::Info | LogLevel::Warning | LogLevel::Error | LogLevel::Critical, true);
  }

  std::cout << "Set log levels to: " << get_enabled_log_levels() << std::endl;
}

static void print_usage(const char* program_name)
{
  std::cout << "Usage: " << program_name << " [OPTIONS] COMMAND\n"
            << "\n"
            << "Options:\n"
            << "  --log-level <levels>      Set log levels (comma-separated: debug,monitor,info,warning,error,critical)\n"
            << "                            Replaces current levels. Default: info,warning,error,critical (no debug/monitor)\n"
            << "  --log-add <levels>        Enable additional log levels (relative to current)\n"
            << "  --log-remove <levels>     Disable specific log levels (relative to current)\n"
            << "  --log-file <path>         Write logs to file (default: no log file, console only)\n"
            << "  --export-mode <mode>      Also export combined mesh when set to combined (default: raw)\n"
            << "  --export-path <path>      Base output directory; writes segment_meshlets/ and raw_meshlets/ beneath it\n"
            << "  --untransformed           Profile-space mesh (no object-space transform, no kinetic separation);\n"
            << "                            applies Y/Z swap for viewing alongside SVG (default: world space)\n"
            << "  --transform-at-construction  Store vertices in object space at add time (default)\n"
            << "  --transform-at-export        Keep vertices in profile space until OBJ export\n"
            << "  --validate                After meshing, check that vertices with the same metadata source agree in world space\n"
            << "  --validate-log <path>     Validation report file (default: mesh_vertex_validation.log)\n"
            << "  --section-shading         Alternate light/dark green and brown materials by even/odd section\n"
            << "  --start <section>         Start kinetic meshing at this section index (default: 0)\n"
            << "  --end <section>           Exclusive stop/finalize time (default: tree height).\n"
            << "                            Section events run on [start, end); events with t >= end are skipped.\n"
            << "  --cutoff <value>          Alpha / radius-event circumradius cutoff (default: 10)\n"
            << "  --debug-files [path]      Write full debug SVGs/TXTs (segmentbuilder snapshots, branch-split dumps,\n"
            << "                            and error dumps). Optional output directory (default: cwd).\n"
            << "                            Off by default.\n"
            << "  --svg-separate-pending-splits\n"
            << "                            With --debug-files: write pending split-off child branches into their own\n"
            << "                            SVG folders (own strand subsets) as soon as a radius event notes the split,\n"
            << "                            not only after the graph cut. Off by default.\n"
            << "  --error-files [path]      Write failure SVGs/TXTs (ring-walk / triangulate FAIL, etc.) plus the\n"
            << "                            common event-style kinetic SVG with affected sites/edges/VVs highlighted,\n"
            << "                            without full visual debug. Optional output directory (default: cwd).\n"
            << "                            Off by default; also enabled by --debug-files.\n"
            << "\n"
            << "Commands:\n"
            << "  --demo                    Run the kinetic Delaunay example\n"
            << "  --mesh <file>             Load StrandTree from file and run TreeMesher\n"
            << "  --help, -h                Show this help message\n"
            << "\n"
            << "Examples:\n"
            << "  " << program_name << " --demo\n"
            << "  " << program_name << " --mesh strandtree.txt\n"
            << "  " << program_name << " --validate --mesh strandtree.txt\n"
            << "  " << program_name << " --log-file output.log --demo\n"
            << "  " << program_name << " --log-level debug,info --log-file debug.log --demo\n"
            << "  " << program_name << " --log-add debug --log-file mesh.log --mesh strandtree.txt\n"
            << "  " << program_name << " --export-mode combined --mesh strandtree.txt combined.obj\n";
}

static const std::vector<std::string>& known_cli_flags()
{
  static const std::vector<std::string> flags = {
    "--log-level",
    "--log-add",
    "--log-remove",
    "--log-file",
    "--export-mode",
    "--export-path",
    "--untransformed",
    "--transform-at-construction",
    "--transform-at-export",
    "--validate",
    "--validate-log",
    "--section-shading",
    "--start",
    "--end",
    "--cutoff",
    "--debug-files",
    "--svg-separate-pending-splits",
    "--error-files",
    "--help",
    "-h",
    "--demo",
    "--mesh",
  };
  return flags;
}

/// Normalize user typos like "debug-files" / "-debug-files" to the canonical "--..." form when unique.
static std::optional<std::string> suggest_cli_flag(const std::string& arg)
{
  if (arg.empty())
  {
    return std::nullopt;
  }

  std::string stripped = arg;
  while (!stripped.empty() && stripped.front() == '-')
  {
    stripped.erase(stripped.begin());
  }
  if (stripped.empty())
  {
    return std::nullopt;
  }

  std::vector<std::string> exact;
  std::vector<std::string> prefix;
  for (const std::string& flag : known_cli_flags())
  {
    std::string flag_body = flag;
    while (!flag_body.empty() && flag_body.front() == '-')
    {
      flag_body.erase(flag_body.begin());
    }
    if (flag_body == stripped)
    {
      exact.push_back(flag);
    }
    else if (flag_body.size() > stripped.size() && flag_body.compare(0, stripped.size(), stripped) == 0)
    {
      prefix.push_back(flag);
    }
  }

  if (exact.size() == 1)
  {
    return exact.front();
  }
  if (exact.empty() && prefix.size() == 1)
  {
    return prefix.front();
  }
  return std::nullopt;
}

static void report_unknown_cli_argument(const std::string& arg)
{
  std::cerr << "Error: Unknown option or command: " << arg << std::endl;
  if (!arg.empty() && arg.front() != '-')
  {
    std::cerr << "Note: options and commands require a leading \"--\" (e.g. --" << arg << ")." << std::endl;
  }
  if (const std::optional<std::string> suggestion = suggest_cli_flag(arg); suggestion.has_value())
  {
    std::cerr << "Did you mean " << suggestion.value() << "?" << std::endl;
  }
}

static bool parse_meshlet_export_mode(const std::string& value, kinDS::MeshletExportMode& out)
{
  if (value == "combined")
  {
    out = kinDS::MeshletExportMode::Combined;
    return true;
  }
  if (value == "segments" || value == "per_segment" || value == "per-segment")
  {
    out = kinDS::MeshletExportMode::PerSegment;
    return true;
  }
  if (value == "raw")
  {
    out = kinDS::MeshletExportMode::Raw;
    return true;
  }
  return false;
}

static std::filesystem::path mesh_export_base_directory(const std::optional<std::filesystem::path>& export_path)
{
  if (!export_path.has_value())
  {
    return std::filesystem::current_path();
  }

  const std::filesystem::path& path = *export_path;
  if (path.has_extension() && path.extension() == ".obj")
  {
    if (path.has_parent_path() && !path.parent_path().empty())
    {
      return path.parent_path();
    }
    return std::filesystem::current_path();
  }
  return path;
}

static bool mesh_from_file(const std::string& filename, kinDS::MeshletExportMode export_mode,
  const std::optional<std::filesystem::path>& export_path, bool profile_space_export,
  bool transform_mesh_at_construction, bool validate_mesh_vertex_sources,
  const std::string& validate_log_path, bool alternate_section_shading, size_t start_section,
  std::optional<size_t> end_section, double alpha_cutoff, bool visual_debug, bool error_files,
  bool visual_debug_separate_pending_splits, const std::optional<std::filesystem::path>& visual_debug_output_root)
{
  std::cout << "Loading StrandTree from: " << filename << std::endl;

  try
  {
    kinDS::StrandTree strand_tree = kinDS::StrandTree::loadFromFile(filename);
    std::cout << "StrandTree loaded successfully. Height: " << strand_tree.getHeight()
              << ", Number of strands: " << strand_tree.getPoints().size();

    const auto& branches_by_height = strand_tree.getStrandsByBranchId();
    if (!branches_by_height.empty())
    {
      size_t min_branches = branches_by_height.front().size();
      size_t max_branches = min_branches;
      for (const auto& at_height : branches_by_height)
      {
        min_branches = std::min(min_branches, at_height.size());
        max_branches = std::max(max_branches, at_height.size());
      }
      std::cout << ", Number of branches: ";
      if (min_branches == max_branches)
      {
        std::cout << min_branches;
      }
      else
      {
        std::cout << min_branches << "-" << max_branches << " (varies by height)";
      }
    }
    std::cout << std::endl;

    std::cout << "Running TreeMesher..." << std::endl;
    kinDS::TreeMesher mesher(strand_tree);
    mesher.getSettings().transform_mesh_at_construction = transform_mesh_at_construction;
    mesher.getSettings().validate_mesh_vertex_sources = validate_mesh_vertex_sources;
    mesher.getSettings().validate_mesh_vertex_sources_log_path = validate_log_path;
    mesher.getSettings().alternate_section_shading = alternate_section_shading;
    mesher.getSettings().start_section = start_section;
    mesher.getSettings().end_section = end_section;
    mesher.getSettings().alpha_cutoff = alpha_cutoff;
    mesher.getSettings().visual_debug = visual_debug;
    mesher.getSettings().error_files = error_files || visual_debug;
    mesher.getSettings().visual_debug_separate_pending_splits = visual_debug_separate_pending_splits;
    if (visual_debug_output_root.has_value())
    {
      mesher.getSettings().visual_debug_output_root = visual_debug_output_root;
    }
    std::cout << "Alpha cutoff: " << alpha_cutoff << std::endl;
    if (visual_debug || error_files)
    {
      const std::filesystem::path debug_root
        = visual_debug_output_root.value_or(std::filesystem::current_path());
      if (visual_debug)
      {
        std::cout << "Debug files: enabled (root=" << debug_root.string() << ")" << std::endl;
        if (visual_debug_separate_pending_splits)
        {
          std::cout << "SVG separate pending splits: enabled" << std::endl;
        }
      }
      else
      {
        std::cout << "Error files: enabled (root=" << debug_root.string() << ")" << std::endl;
      }
    }
    if (start_section != 0 || end_section.has_value())
    {
      std::cout << "Section range: start=" << start_section;
      if (end_section.has_value())
      {
        std::cout << " end=" << end_section.value();
      }
      else
      {
        std::cout << " end=<last>";
      }
      std::cout << std::endl;
    }
    const auto& meshes = mesher.runMeshingAlgorithm(visual_debug);

    std::cout << "Meshing completed. Generated " << meshes.size() << " meshlets." << std::endl;

    const std::filesystem::path export_base = mesh_export_base_directory(export_path);
    const std::filesystem::path segment_export_dir = export_base / "segment_meshlets";
    const std::filesystem::path raw_export_dir = export_base / "raw_meshlets";
    const std::optional<bool> export_transform
      = profile_space_export ? std::optional<bool>(false) : std::nullopt;
    const bool apply_export_transform
      = export_transform.value_or(!mesher.meshletsTransformedAtConstruction());
    const char* space_label = profile_space_export ? "profile space" : "world space";
    const char* transform_label = mesher.meshletsTransformedAtConstruction() && !profile_space_export
      ? "transformed at construction"
      : apply_export_transform ? "transformed at export" : "no export transform";

    std::cout << "Exporting segment meshlets to: " << segment_export_dir.string() << " (" << space_label << ", "
              << transform_label << (profile_space_export ? ", _untransformed suffix" : "") << ")" << std::endl;
    mesher.exportMeshlets(kinDS::MeshletExportMode::PerSegment, segment_export_dir, export_transform);

    std::cout << "Exporting raw meshlets to: " << raw_export_dir.string() << " (" << space_label << ", "
              << transform_label << (profile_space_export ? ", _untransformed suffix" : "") << ")" << std::endl;
    mesher.exportMeshlets(kinDS::MeshletExportMode::Raw, raw_export_dir, export_transform);

    if (export_mode == kinDS::MeshletExportMode::Combined)
    {
      const std::filesystem::path combined_export_path = export_path.has_value() && export_path->has_extension()
        && export_path->extension() == ".obj"
        ? *export_path
        : export_base / "combined_mesh.obj";
      std::cout << "Exporting combined mesh to: " << combined_export_path.string()
                << (profile_space_export ? " (with _untransformed suffix)" : "") << " (" << space_label << ", "
                << transform_label << ")" << std::endl;
      mesher.exportMeshlets(kinDS::MeshletExportMode::Combined, combined_export_path, export_transform);
    }

    // Export boundary mesh
    kinDS::VoronoiMesh boundary_mesh = mesher.getBoundaryMesh();
    if (apply_export_transform)
    {
      mesher.transformBoundaryMesh(boundary_mesh);
    }
    const std::filesystem::path boundary_export_path
      = profile_space_export ? "boundary_mesh_untransformed.obj" : "boundary_mesh.obj";
    kinDS::ObjExporter::writeMesh(boundary_mesh, boundary_export_path, 1.0, 1.0, {}, true);
    std::cout << "Boundary mesh exported to: " << boundary_export_path.string() << std::endl;

    std::cout << "Meshing complete!" << std::endl;

    if (validate_mesh_vertex_sources)
    {
      std::cout << "Validation report written to: " << kinDS::Validator::logFilePath() << std::endl;
      if (!mesher.meshVertexSourceValidationPassed())
      {
        std::cerr << "Mesh vertex source validation failed. See " << kinDS::Validator::logFilePath() << std::endl;
        return false;
      }
    }
    return true;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return false;
  }
}

int main(int argc, char* argv[])
{
  // Sanity check: print all command line arguments
  std::cout << "Arguments (" << argc << "):";
  for (int i = 0; i < argc; ++i)
  {
    std::cout << " [" << i << "]=\"" << argv[i] << "\"";
  }
  std::cout << std::endl;

  if (argc < 2)
  {
    print_usage(argv[0]);
    return 1;
  }

  // First pass: parse all options, remember the command (order-independent)
  std::string command; // "demo", "mesh", "help"
  std::string mesh_file; // for --mesh
  kinDS::MeshletExportMode mesh_export_mode = kinDS::MeshletExportMode::Raw;
  std::optional<std::filesystem::path> mesh_export_path;
  bool mesh_export_profile_space = false;
  bool mesh_transform_at_construction = true;
  bool mesh_validate_vertex_sources = false;
  std::string mesh_validate_log_path = kinDS::Validator::defaultLogFilePath();
  bool mesh_alternate_section_shading = false;
  size_t mesh_start_section = 0;
  std::optional<size_t> mesh_end_section;
  double mesh_alpha_cutoff = kinDS::TreeMesher::Settings {}.alpha_cutoff;
  bool mesh_visual_debug = false;
  bool mesh_error_files = false;
  bool mesh_visual_debug_separate_pending_splits = false;
  std::optional<std::filesystem::path> mesh_visual_debug_output_root;

  int arg_idx = 1;
  while (arg_idx < argc)
  {
    std::string arg = argv[arg_idx];

    if (arg == "--log-level")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --log-level requires a value." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      set_log_level(argv[arg_idx + 1]);
      arg_idx += 2; // Skip both --log-level and its value
    }
    else if (arg == "--log-add")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --log-add requires a value." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      if (!modify_log_level(argv[arg_idx + 1], true))
      {
        std::cerr << "Warning: --log-add did not recognize any valid log levels." << std::endl;
      }
      else
      {
        std::cout << "Added log levels: " << argv[arg_idx + 1] << " -> now enabled: " << get_enabled_log_levels()
                  << std::endl;
      }
      arg_idx += 2;
    }
    else if (arg == "--log-remove")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --log-remove requires a value." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      if (!modify_log_level(argv[arg_idx + 1], false))
      {
        std::cerr << "Warning: --log-remove did not recognize any valid log levels." << std::endl;
      }
      else
      {
        std::cout << "Removed log levels: " << argv[arg_idx + 1] << " -> now enabled: " << get_enabled_log_levels()
                  << std::endl;
      }
      arg_idx += 2;
    }
    else if (arg == "--log-file")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --log-file requires a file path." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      kinDS::logger.setLogFile(argv[arg_idx + 1]);
      std::cout << "Log file set to: " << argv[arg_idx + 1] << std::endl;
      arg_idx += 2;
    }
    else if (arg == "--export-mode")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --export-mode requires a value (combined, segments, raw)." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      if (!parse_meshlet_export_mode(argv[arg_idx + 1], mesh_export_mode))
      {
        std::cerr << "Error: Unknown export mode: " << argv[arg_idx + 1]
                  << " (expected combined, segments, or raw)." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      arg_idx += 2;
    }
    else if (arg == "--export-path")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --export-path requires a path." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      mesh_export_path = argv[arg_idx + 1];
      arg_idx += 2;
    }
    else if (arg == "--untransformed")
    {
      mesh_export_profile_space = true;
      mesh_transform_at_construction = false;
      ++arg_idx;
    }
    else if (arg == "--transform-at-construction")
    {
      mesh_transform_at_construction = true;
      ++arg_idx;
    }
    else if (arg == "--transform-at-export")
    {
      mesh_transform_at_construction = false;
      ++arg_idx;
    }
    else if (arg == "--validate")
    {
      mesh_validate_vertex_sources = true;
      ++arg_idx;
    }
    else if (arg == "--validate-log")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --validate-log requires a file path." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      mesh_validate_log_path = argv[arg_idx + 1];
      ++arg_idx;
      ++arg_idx;
    }
    else if (arg == "--section-shading")
    {
      mesh_alternate_section_shading = true;
      ++arg_idx;
    }
    else if (arg == "--start")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --start requires a section index." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      try
      {
        const long long parsed = std::stoll(argv[arg_idx + 1]);
        if (parsed < 0)
        {
          throw std::out_of_range("negative");
        }
        mesh_start_section = static_cast<size_t>(parsed);
      }
      catch (const std::exception&)
      {
        std::cerr << "Error: --start expects a non-negative integer section index." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      arg_idx += 2;
    }
    else if (arg == "--end")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --end requires a section index." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      try
      {
        const long long parsed = std::stoll(argv[arg_idx + 1]);
        if (parsed < 0)
        {
          throw std::out_of_range("negative");
        }
        mesh_end_section = static_cast<size_t>(parsed);
      }
      catch (const std::exception&)
      {
        std::cerr << "Error: --end expects a non-negative integer section index." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      arg_idx += 2;
    }
    else if (arg == "--cutoff")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --cutoff requires a numeric value." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      try
      {
        const double parsed = std::stod(argv[arg_idx + 1]);
        if (!(parsed >= 0.0) || !std::isfinite(parsed))
        {
          throw std::out_of_range("non-finite or negative");
        }
        mesh_alpha_cutoff = parsed;
      }
      catch (const std::exception&)
      {
        std::cerr << "Error: --cutoff expects a non-negative finite number." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      arg_idx += 2;
    }
    else if (arg == "--debug-files")
    {
      mesh_visual_debug = true;
      // Optional directory: consume next arg only if it is present and not another option/command.
      if (arg_idx + 1 < argc)
      {
        const std::string next = argv[arg_idx + 1];
        if (!next.empty() && next[0] != '-')
        {
          mesh_visual_debug_output_root = std::filesystem::path(next);
          arg_idx += 2;
          continue;
        }
      }
      arg_idx += 1;
    }
    else if (arg == "--svg-separate-pending-splits")
    {
      mesh_visual_debug_separate_pending_splits = true;
      arg_idx += 1;
    }
    else if (arg == "--error-files")
    {
      mesh_error_files = true;
      if (arg_idx + 1 < argc)
      {
        const std::string next = argv[arg_idx + 1];
        if (!next.empty() && next[0] != '-')
        {
          mesh_visual_debug_output_root = std::filesystem::path(next);
          arg_idx += 2;
          continue;
        }
      }
      arg_idx += 1;
    }
    else if (arg == "--help" || arg == "-h")
    {
      if (!command.empty() && command != "help")
      {
        std::cerr << "Error: Multiple commands specified (existing: " << command << ", new: help)." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      command = "help";
      ++arg_idx;
    }
    else if (arg == "--demo")
    {
      if (!command.empty() && command != "demo")
      {
        std::cerr << "Error: Multiple commands specified (existing: " << command << ", new: demo)." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      command = "demo";
      ++arg_idx;
    }
    else if (arg == "--mesh")
    {
      if (arg_idx + 1 >= argc)
      {
        std::cerr << "Error: --mesh requires a filename argument." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      if (!command.empty() && command != "mesh")
      {
        std::cerr << "Error: Multiple commands specified (existing: " << command << ", new: mesh)." << std::endl;
        print_usage(argv[0]);
        return 1;
      }
      command = "mesh";
      mesh_file = argv[arg_idx + 1];
      arg_idx += 2;
    }
    else
    {
      report_unknown_cli_argument(arg);
      print_usage(argv[0]);
      return 1;
    }
  }

  if (command.empty())
  {
    std::cerr << "Error: No command specified." << std::endl;
    print_usage(argv[0]);
    return 1;
  }

  if (mesh_end_section.has_value() && mesh_end_section.value() < mesh_start_section)
  {
    std::cerr << "Error: --end (" << mesh_end_section.value() << ") is before --start (" << mesh_start_section
              << ")." << std::endl;
    return 1;
  }

  // Second pass: execute the chosen command (options already applied)
  if (command == "help")
  {
    print_usage(argv[0]);
    return 0;
  }
  else if (command == "demo")
  {
    std::cout << "Running demo (kinetic_delaunay_example)..." << std::endl;
    if (!kinetic_delaunay_example(mesh_validate_vertex_sources, mesh_validate_log_path, mesh_alternate_section_shading))
    {
      return 1;
    }
    return 0;
  }
  else if (command == "mesh")
  {
    std::cout << "Running TreeMesher on file: " << mesh_file << std::endl;
    if (!mesh_from_file(mesh_file, mesh_export_mode, mesh_export_path, mesh_export_profile_space,
          mesh_transform_at_construction, mesh_validate_vertex_sources, mesh_validate_log_path,
          mesh_alternate_section_shading, mesh_start_section, mesh_end_section, mesh_alpha_cutoff, mesh_visual_debug,
          mesh_error_files, mesh_visual_debug_separate_pending_splits, mesh_visual_debug_output_root))
    {
      return 1;
    }
    return 0;
  }

  // Should not reach here
  std::cerr << "Error: Unknown command state: " << command << std::endl;
  return 1;
}
