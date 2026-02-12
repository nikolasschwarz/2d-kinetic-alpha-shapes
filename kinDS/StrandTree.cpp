#include "StrandTree.hpp"

#include "Logger.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace kinDS;

static glm::dvec3 ProfileToModelCoordinatesBranch(
  const std::vector<std::vector<glm::mat4>>& profile_to_model_transforms, glm::dvec3 point, float t,
  const std::vector<size_t>& branch_indices, float w = 1.0f)
{
  size_t lower_section_index = static_cast<size_t>(std::max(0.0f, glm::floor(t)));

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
  glm::vec4 local_pos(point[0], (1.0f - w) * point[2], point[1], w);
  if (lower_section_index >= branch_indices.size())
  {
    std::cout << ("ProfileToModelCoordinates: lower bound of point z-coordinate out of range: " + coord_str).c_str()
              << ", index: " << lower_section_index << ", size: " << branch_indices.size() << std::endl;
  }
  size_t lower_branch_index = branch_indices[lower_section_index];
  glm::vec4 global_pos = profile_to_model_transforms[lower_section_index][lower_branch_index] * local_pos;

  if (upper_section_index != lower_section_index)
  {
    size_t upper_branch_index = branch_indices[upper_section_index];
    glm::vec4 upper_global_pos = profile_to_model_transforms[upper_section_index][upper_branch_index] * local_pos;
    double frac = t - static_cast<double>(lower_section_index);
    global_pos = glm::mix(global_pos, upper_global_pos, frac);
  }

  if (w == 0.0f)
  {
    global_pos = glm::normalize(global_pos);
  }

  return { global_pos.x, global_pos.y, global_pos.z };
}

StrandTree::StrandTree(const std::vector<std::vector<glm::dvec2>>& support_points,
  const std::vector<std::vector<double>>& subdivisions_by_strand,
  const std::vector<std::vector<int>>& physics_strand_to_segment_indices,
  const std::vector<std::vector<glm::mat4>>& transforms_by_height_and_branch,
  const std::vector<std::vector<size_t>>& branch_indices,
  const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id)
  : support_points(support_points)
  , subdivisions_by_strand(subdivisions_by_strand)
  , physics_strand_to_segment_indices(physics_strand_to_segment_indices)
  , transforms_by_height_and_branch(transforms_by_height_and_branch)
  , branch_indices(branch_indices)
  , strands_by_branch_id(strands_by_branch_id)
{
  // compute height
  for (const auto& pts : support_points)
  {
    // always one less than size because the parameter is in range from smallest to largest index
    height = std::max(height, pts.size() - 1);
  }
}

const std::vector<std::vector<glm::dvec2>>& StrandTree::getPoints() const { return support_points; }

size_t StrandTree::getHeight() const { return height; }

size_t StrandTree::addTrajectory(const std::vector<glm::dvec2>& traj)
{
  size_t index = support_points.size();
  support_points.push_back(traj);
  return index;
}

glm::dvec2 StrandTree::evaluate(size_t strand_id, double t) const
{
  if (t < 0)
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  size_t lower_index = std::floor(t);
  size_t upper_index = lower_index + 1;
  double frac = t - lower_index;

  if (upper_index >= support_points[strand_id].size() && lower_index < support_points[strand_id].size()
    && frac < std::numeric_limits<double>::epsilon())
  {
    return support_points[strand_id].back();
  }

  if (upper_index >= support_points[strand_id].size())
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  const glm::dvec2& lower = support_points[strand_id][lower_index];
  const glm::dvec2& upper = support_points[strand_id][upper_index];

  return lower * (1.0 - frac) + upper * frac;
}

glm::dvec2 StrandTree::getPointTransformed(size_t strand_id, size_t index, size_t reference_branch) const
{
  size_t actual_branch;
  // dummy strands might not be mapped to a branch:
  if (strand_id >= branch_indices.size())
  {
    actual_branch = reference_branch;
  }
  else
  {
    actual_branch = branch_indices[strand_id][index];
  }

  glm::dvec2 point = support_points[strand_id][index];
  if (actual_branch == reference_branch)
  {
    return point;
  }

  PlaneProjector plane_projector(
    transforms_by_height_and_branch[index][actual_branch], transforms_by_height_and_branch[index][reference_branch]);

  auto result = plane_projector.project(glm::vec2(point[0], point[1]));
  return glm::dvec2 { result.x, result.y };
}

glm::dvec2 StrandTree::evaluateTransformed(size_t strand_id, double t, size_t reference_branch) const
{
  if (t < 0)
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  size_t lower_index = std::floor(t);
  size_t upper_index = lower_index + 1;
  double frac = t - lower_index;

  if (upper_index >= support_points[strand_id].size() && lower_index < support_points[strand_id].size()
    && frac < std::numeric_limits<double>::epsilon())
  {
    return getPointTransformed(strand_id, support_points[strand_id].size() - 1, reference_branch);
  }

  if (upper_index >= support_points[strand_id].size())
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  const glm::dvec2& lower = getPointTransformed(strand_id, lower_index, reference_branch);
  const glm::dvec2& upper = getPointTransformed(strand_id, upper_index, reference_branch);

  return lower * (1.0 - frac) + upper * frac;
}

glm::dvec3 StrandTree::getPointInObjectSpace(size_t strand_id, double t) const
{
  glm::dvec2 v = evaluate(strand_id, t);

  glm::dvec3 v_3d { v[0], 0.0, v[1] };
  return ProfileToModelCoordinatesBranch(transforms_by_height_and_branch, v_3d, t, branch_indices[strand_id]);
}

glm::dvec3 StrandTree::transformToObjectSpace(glm::dvec3& v_3d, size_t strand_id, double t) const
{
  return ProfileToModelCoordinatesBranch(transforms_by_height_and_branch, v_3d, t, branch_indices[strand_id]);
}

// TODO: also adjust to different reference frame
std::array<Polynomial, 2> StrandTree::getPiecePolynomial(size_t strand_id, size_t index) const
{
  if (strand_id >= support_points.size())
  {
    throw std::out_of_range("Strand id " + std::to_string(strand_id) + " out of range.");
  }

  if (index >= support_points[strand_id].size() - 1)
  {
    throw std::out_of_range("Index " + std::to_string(index) + " out of range for piece polynomial.");
  }
  const auto& P0 = support_points[strand_id][index];
  const auto& P1 = support_points[strand_id][index + 1];

  std::array<Polynomial, 2> result;
  // Create linear polynomials for each dimension
  for (int i = 0; i < 2; ++i)
  {
    result[i] = POLYNOMIAL(P0[i] + (P1[i] - P0[i]) * x);
  }

  return result;
}

namespace
{
const char* FORMAT_MAGIC = "StrandTree 1";

void writeVec2(std::ostream& out, const glm::dvec2& v) { out << v.x << ' ' << v.y << '\n'; }

void readVec2(std::istream& in, glm::dvec2& v) { in >> v.x >> v.y; }

void writeMat4(std::ostream& out, const glm::mat4& m)
{
  for (int row = 0; row < 4; ++row)
    for (int col = 0; col < 4; ++col)
      out << (col ? " " : "") << m[col][row];
  out << '\n';
}

void readMat4(std::istream& in, glm::mat4& m)
{
  for (int row = 0; row < 4; ++row)
    for (int col = 0; col < 4; ++col)
      in >> m[col][row];
}
} // namespace

void StrandTree::saveToFile(const std::filesystem::path& path) const
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error(
      std::string("StrandTree::saveToFile: cannot open ") + path.string() + std::string(" for writing"));
  out << FORMAT_MAGIC << '\n';
  out << height << '\n';
  out << support_points.size() << '\n';

  out << "support_points\n";
  for (const auto& strand : support_points)
  {
    out << strand.size() << '\n';
    for (const auto& p : strand)
      writeVec2(out, p);
  }

  out << "subdivisions_by_strand\n";
  for (const auto& strand : subdivisions_by_strand)
  {
    out << strand.size();
    for (double d : strand)
      out << ' ' << d;
    out << '\n';
  }

  out << "physics_strand_to_segment_indices\n";
  for (const auto& strand : physics_strand_to_segment_indices)
  {
    out << strand.size();
    for (int i : strand)
      out << ' ' << i;
    out << '\n';
  }

  out << "transforms_by_height_and_branch\n";
  out << transforms_by_height_and_branch.size() << '\n';
  for (const auto& by_branch : transforms_by_height_and_branch)
  {
    out << by_branch.size() << '\n';
    for (const auto& m : by_branch)
      writeMat4(out, m);
  }

  out << "branch_indices\n";
  for (const auto& strand : branch_indices)
  {
    out << strand.size();
    for (size_t s : strand)
      out << ' ' << s;
    out << '\n';
  }

  out << "strands_by_branch_id\n";
  out << strands_by_branch_id.size() << '\n';
  for (const auto& by_height : strands_by_branch_id)
  {
    out << by_height.size() << '\n';
    for (const auto& branch : by_height)
    {
      out << branch.size();
      for (size_t s : branch)
        out << ' ' << s;
      out << '\n';
    }
  }
}

StrandTree StrandTree::loadFromFile(const std::filesystem::path& path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error(
      std::string("StrandTree::loadFromFile: cannot open ") + path.string() + std::string(" for reading"));
  std::string line;
  auto getLine = [&in, &line]() -> bool { return static_cast<bool>(std::getline(in, line)); };
  auto expect = [&line](const char* tag)
  {
    if (line != tag)
      throw std::runtime_error(std::string("StrandTree::loadFromFile: expected '") + tag + "', got '" + line + "'");
  };

  if (!getLine() || line != FORMAT_MAGIC)
    throw std::runtime_error("StrandTree::loadFromFile: invalid or unsupported format");
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  size_t read_height = 0;
  std::istringstream(line) >> read_height;
  size_t num_strands = 0;
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  std::istringstream(line) >> num_strands;

  std::vector<std::vector<glm::dvec2>> support_points;
  support_points.reserve(num_strands);
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("support_points");
  for (size_t s = 0; s < num_strands; ++s)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (support_points count)");
    size_t n = 0;
    std::istringstream(line) >> n;
    std::vector<glm::dvec2> strand;
    strand.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      if (!getLine())
        throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (support point)");
      glm::dvec2 p;
      std::istringstream iss(line);
      readVec2(iss, p);
      strand.push_back(p);
    }
    support_points.push_back(std::move(strand));
  }

  std::vector<std::vector<double>> subdivisions_by_strand;
  subdivisions_by_strand.reserve(num_strands);
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("subdivisions_by_strand");
  for (size_t s = 0; s < num_strands; ++s)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (subdivisions)");
    std::istringstream iss(line);
    size_t n = 0;
    iss >> n;
    std::vector<double> strand;
    double d;
    while (iss >> d)
      strand.push_back(d);
    if (strand.size() != n)
      throw std::runtime_error("StrandTree::loadFromFile: subdivisions count mismatch");
    subdivisions_by_strand.push_back(std::move(strand));
  }

  std::vector<std::vector<int>> physics_strand_to_segment_indices;
  physics_strand_to_segment_indices.reserve(num_strands);
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("physics_strand_to_segment_indices");
  for (size_t s = 0; s < num_strands; ++s)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (physics indices)");
    std::istringstream iss(line);
    size_t n = 0;
    iss >> n;
    std::vector<int> strand;
    int i;
    while (iss >> i)
      strand.push_back(i);
    if (strand.size() != n)
      throw std::runtime_error("StrandTree::loadFromFile: physics_strand_to_segment_indices count mismatch");
    physics_strand_to_segment_indices.push_back(std::move(strand));
  }

  std::vector<std::vector<glm::mat4>> transforms_by_height_and_branch;
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("transforms_by_height_and_branch");
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  size_t num_heights = 0;
  std::istringstream(line) >> num_heights;
  transforms_by_height_and_branch.reserve(num_heights);
  for (size_t h = 0; h < num_heights; ++h)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (transforms branch count)");
    size_t num_branches = 0;
    std::istringstream(line) >> num_branches;
    std::vector<glm::mat4> by_branch;
    by_branch.reserve(num_branches);
    for (size_t b = 0; b < num_branches; ++b)
    {
      if (!getLine())
        throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (transform matrix)");
      glm::mat4 m;
      std::istringstream iss(line);
      readMat4(iss, m);
      by_branch.push_back(m);
    }
    transforms_by_height_and_branch.push_back(std::move(by_branch));
  }

  std::vector<std::vector<size_t>> branch_indices;
  branch_indices.reserve(num_strands);
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("branch_indices");
  for (size_t s = 0; s < num_strands; ++s)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (branch indices)");
    std::istringstream iss(line);
    size_t n = 0;
    iss >> n;
    std::vector<size_t> strand;
    size_t val;
    while (iss >> val)
      strand.push_back(val);
    if (strand.size() != n)
      throw std::runtime_error("StrandTree::loadFromFile: branch_indices count mismatch");
    branch_indices.push_back(std::move(strand));
  }

  std::vector<std::vector<std::vector<size_t>>> strands_by_branch_id;
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  expect("strands_by_branch_id");
  if (!getLine())
    throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF");
  num_heights = 0;
  std::istringstream(line) >> num_heights;
  strands_by_branch_id.reserve(num_heights);
  for (size_t h = 0; h < num_heights; ++h)
  {
    if (!getLine())
      throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (strands_by_branch branch count)");
    size_t num_branches = 0;
    std::istringstream(line) >> num_branches;
    std::vector<std::vector<size_t>> by_branch;
    by_branch.reserve(num_branches);
    for (size_t b = 0; b < num_branches; ++b)
    {
      if (!getLine())
        throw std::runtime_error("StrandTree::loadFromFile: unexpected EOF (strands list)");
      std::istringstream iss(line);
      size_t n = 0;
      iss >> n;
      std::vector<size_t> branch;
      size_t val;
      while (iss >> val)
        branch.push_back(val);
      if (branch.size() != n)
        throw std::runtime_error("StrandTree::loadFromFile: strands_by_branch_id count mismatch");
      by_branch.push_back(std::move(branch));
    }
    strands_by_branch_id.push_back(std::move(by_branch));
  }

  StrandTree tree(support_points, subdivisions_by_strand, physics_strand_to_segment_indices,
    transforms_by_height_and_branch, branch_indices, strands_by_branch_id);

  // Validate height: compare stored height with derived height
  size_t derived_height = tree.getHeight();
  if (read_height != derived_height)
  {
    std::ostringstream ss;
    ss << "StrandTree::loadFromFile: Height mismatch - file contains height " << read_height
       << " but derived height from support_points is " << derived_height;
    KINDS_WARNING(ss.str());
  }

  return tree;
}