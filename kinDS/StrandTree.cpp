#include "StrandTree.hpp"

#include "Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace kinDS;

namespace
{
bool profileAxesAreUsable(const glm::dmat4& transform)
{
  return glm::length(glm::dvec3(transform[0])) > 0.0 && glm::length(glm::dvec3(transform[2])) > 0.0;
}

const glm::dmat4& requireProfileTransform(
  const std::vector<std::vector<glm::dmat4>>& transforms, size_t height, size_t branch_id)
{
  if (height >= transforms.size())
  {
    throw std::out_of_range(
      "StrandTree: no profile transform at height " + std::to_string(height) + " (count=" + std::to_string(transforms.size()) + ")");
  }
  if (branch_id >= transforms[height].size())
  {
    throw std::out_of_range("StrandTree: no profile transform at height " + std::to_string(height) + " branch "
      + std::to_string(branch_id) + " (count=" + std::to_string(transforms[height].size()) + ")");
  }
  const glm::dmat4& transform = transforms[height][branch_id];
  if (!profileAxesAreUsable(transform))
  {
    throw std::runtime_error("StrandTree: profile transform at height " + std::to_string(height) + " branch "
      + std::to_string(branch_id) + " has zero profile axes");
  }
  return transform;
}
} // namespace

static std::string profileInputBranchSuffix(const std::vector<size_t>& branch_indices, size_t requested_height)
{
  if (branch_indices.empty())
  {
    return " input_branch=none";
  }
  return " input_branch=" + std::to_string(branch_indices[std::min(requested_height, branch_indices.size() - 1)]);
}

static glm::dvec3 ProfileToModelCoordinatesBranch(
  const std::vector<std::vector<glm::dmat4>>& profile_to_model_transforms, glm::dvec3 point, float t,
  const std::vector<size_t>& branch_indices, double w = 1.0)
{
  if (profile_to_model_transforms.empty())
  {
    throw std::out_of_range("ProfileToModelCoordinates: no profile transforms");
  }

  const size_t lower_section_index = static_cast<size_t>(std::max(0.0f, glm::floor(t)));
  const size_t upper_section_index = static_cast<size_t>(glm::ceil(t));
  if (lower_section_index >= profile_to_model_transforms.size() || lower_section_index >= branch_indices.size())
  {
    throw std::out_of_range("ProfileToModelCoordinates: t=" + std::to_string(t) + " height " + std::to_string(lower_section_index)
      + " is out of range (transforms=" + std::to_string(profile_to_model_transforms.size())
      + " branch_indices=" + std::to_string(branch_indices.size()) + ")"
      + profileInputBranchSuffix(branch_indices, lower_section_index));
  }

  glm::dvec4 local_pos(point[0], (1.0f - w) * point[2], point[1], w);
  const size_t lower_branch_index = branch_indices[lower_section_index];
  glm::dvec4 global_pos
    = requireProfileTransform(profile_to_model_transforms, lower_section_index, lower_branch_index) * local_pos;

  if (upper_section_index != lower_section_index)
  {
    if (upper_section_index >= profile_to_model_transforms.size() || upper_section_index >= branch_indices.size())
    {
      throw std::out_of_range("ProfileToModelCoordinates: t=" + std::to_string(t) + " height " + std::to_string(upper_section_index)
        + " is out of range (transforms=" + std::to_string(profile_to_model_transforms.size())
        + " branch_indices=" + std::to_string(branch_indices.size()) + ")"
        + profileInputBranchSuffix(branch_indices, upper_section_index));
    }
    const size_t upper_branch_index = branch_indices[upper_section_index];
    const glm::dvec4 upper_global_pos
      = requireProfileTransform(profile_to_model_transforms, upper_section_index, upper_branch_index) * local_pos;
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
  const std::vector<std::vector<glm::dmat4>>& transforms_by_height_and_branch,
  const std::vector<std::vector<size_t>>& branch_indices,
  const std::vector<std::vector<std::vector<size_t>>>& strands_by_branch_id)
  : support_points(support_points)
  , subdivisions_by_strand(subdivisions_by_strand)
  , physics_strand_to_segment_indices(physics_strand_to_segment_indices)
  , transforms_by_height_and_branch(transforms_by_height_and_branch)
  , branch_indices(branch_indices)
  , strands_by_branch_id(strands_by_branch_id)
{
  // Snap all profile transforms to similarities before anything derived from them (normals, projections) is computed.
  conditionProfileTransforms();

  // compute height
  for (const auto& pts : support_points)
  {
    // always one less than size because the parameter is in range from smallest to largest index
    tree_height = std::max(tree_height, pts.size() - 1);
    computeNormalTransforms();
  }
}

const std::vector<std::vector<glm::dvec2>>& StrandTree::getPoints() const { return support_points; }

size_t StrandTree::getHeight() const { return tree_height; }

size_t StrandTree::getBranchCount(size_t height) const
{
  if (height >= strands_by_branch_id.size())
  {
    return 0;
  }
  return strands_by_branch_id[height].size();
}

size_t StrandTree::resolveBranchAtHeight(size_t height, size_t branch_id, size_t branch_lookup_height) const
{
  if (height >= transforms_by_height_and_branch.size())
  {
    throw std::out_of_range("StrandTree::resolveBranchAtHeight: height " + std::to_string(height)
      + " is out of range (count=" + std::to_string(transforms_by_height_and_branch.size()) + ")");
  }

  if (branch_id < getBranchCount(height))
  {
    return branch_id;
  }

  if (strands_by_branch_id.empty())
  {
    throw std::runtime_error("StrandTree::resolveBranchAtHeight: no strand-to-branch tables");
  }

  size_t lookup = branch_lookup_height;
  while (lookup < strands_by_branch_id.size() && branch_id >= getBranchCount(lookup))
  {
    ++lookup;
  }

  if (lookup >= strands_by_branch_id.size() || branch_id >= getBranchCount(lookup))
  {
    throw std::runtime_error("StrandTree::resolveBranchAtHeight: branch " + std::to_string(branch_id)
      + " is not present at height " + std::to_string(height) + " or lookup " + std::to_string(branch_lookup_height));
  }

  const auto& strands = strands_by_branch_id[lookup][branch_id];
  if (strands.empty())
  {
    throw std::runtime_error("StrandTree::resolveBranchAtHeight: branch " + std::to_string(branch_id)
      + " has no strands at lookup height " + std::to_string(lookup));
  }

  const size_t representative_strand = strands.front();
  if (representative_strand >= branch_indices.size() || height >= branch_indices[representative_strand].size())
  {
    throw std::out_of_range("StrandTree::resolveBranchAtHeight: strand " + std::to_string(representative_strand)
      + " has no branch id at height " + std::to_string(height));
  }

  return branch_indices[representative_strand][height];
}

size_t StrandTree::getBranchIndex(size_t strand_id, size_t height) const
{
  if (strand_id >= branch_indices.size())
  {
    throw std::out_of_range("StrandTree::getBranchIndex: strand " + std::to_string(strand_id)
      + " is out of range (count=" + std::to_string(branch_indices.size()) + ")");
  }

  const std::vector<size_t>& per_strand = branch_indices[strand_id];
  if (per_strand.empty())
  {
    throw std::runtime_error("StrandTree::getBranchIndex: strand " + std::to_string(strand_id) + " has no branch ids");
  }

  const size_t clamped_height = std::min(height, per_strand.size() - 1);
  return per_strand[clamped_height];
}

const glm::dmat4& StrandTree::getTransformByHeightAndBranch(size_t height, size_t branch_id) const
{
  return requireProfileTransform(transforms_by_height_and_branch, height, branch_id);
}

const glm::dmat4& StrandTree::getResolvedTransformByHeightAndBranch(
  size_t height, size_t branch_id, size_t branch_lookup_height) const
{
  return getTransformByHeightAndBranch(height, resolveBranchAtHeight(height, branch_id, branch_lookup_height));
}

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

  if (lower_index >= support_points[strand_id].size())
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  if (frac < std::numeric_limits<double>::epsilon())
  {
    return support_points[strand_id][lower_index];
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
  return getPointTransformedAtSection(strand_id, index, reference_branch, index);
}

glm::dvec2 StrandTree::getPointTransformedAtSection(
  size_t strand_id, size_t index, size_t reference_branch, size_t branch_lookup_height) const
{
  size_t actual_branch;
  // dummy strands might not be mapped to a branch:
  if (strand_id >= branch_indices.size())
  {
    actual_branch = reference_branch;
  }
  else
  {
    actual_branch = getBranchIndex(strand_id, index);
  }

  actual_branch = resolveBranchAtHeight(index, actual_branch, branch_lookup_height);
  reference_branch = resolveBranchAtHeight(index, reference_branch, branch_lookup_height);

  glm::dvec2 point = support_points[strand_id][index];
  if (actual_branch == reference_branch)
  {
    return point;
  }

  PlaneProjector plane_projector(
    getTransformByHeightAndBranch(index, reference_branch), getTransformByHeightAndBranch(index, actual_branch));

  glm::dvec2 result = plane_projector.project(glm::dvec2(point[0], point[1]));
  return glm::dvec2 { result.x, result.y };
}

glm::dvec2 StrandTree::evaluateTransformed(size_t strand_id, double t, size_t reference_branch) const
{
  if (t < 0)
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  size_t lower_index = std::floor(t);
  double frac = t - lower_index;

  if (lower_index >= support_points[strand_id].size())
  {
    throw std::runtime_error("Parameter t out of bounds");
  }

  if (frac < std::numeric_limits<double>::epsilon())
  {
    return getPointTransformedAtSection(strand_id, lower_index, reference_branch, lower_index);
  }

  const Trajectory<2> piece = getPiecePolynomial(strand_id, lower_index, reference_branch);
  return glm::dvec2(piece[0](frac), piece[1](frac));
}

glm::dvec3 StrandTree::getPointInObjectSpace(size_t strand_id, double t) const
{
  glm::dvec2 v = evaluate(strand_id, t);

  glm::dvec3 v_3d { v[0], 0.0, v[1] };
  return transformToObjectSpace(v_3d, strand_id, t);
}

glm::dvec3 StrandTree::transformToObjectSpace(glm::dvec3& v_3d, size_t strand_id, double t) const
{
  if (strand_id >= branch_indices.size() || branch_indices[strand_id].empty()
    || static_cast<size_t>(std::ceil(std::max(0.0, t))) >= branch_indices[strand_id].size())
  {
    const size_t samples = strand_id < branch_indices.size() ? branch_indices[strand_id].size() : 0;
    throw std::out_of_range("StrandTree::transformToObjectSpace: strand " + std::to_string(strand_id) + " has "
      + std::to_string(samples) + " profile samples; t=" + std::to_string(t) + " needs height "
      + std::to_string(static_cast<size_t>(std::ceil(std::max(0.0, t))))
      + (samples == 0 ? " input_branch=none" : " input_branch=" + std::to_string(getBranchIndex(strand_id,
          static_cast<size_t>(std::ceil(std::max(0.0, t)))))));
  }
  try
  {
    return transformToObjectSpace(v_3d, t, branch_indices[strand_id]);
  }
  catch (const std::out_of_range& ex)
  {
    throw std::out_of_range(
      std::string(ex.what()) + " strand=" + std::to_string(strand_id)
      + " input_branch=" + std::to_string(getBranchIndex(strand_id, static_cast<size_t>(std::ceil(std::max(0.0, t))))));
  }
}

glm::dvec3 StrandTree::transformToObjectSpace(
  glm::dvec3 v_3d, double t, const std::vector<size_t>& branch_indices_by_height) const
{
  return ProfileToModelCoordinatesBranch(
    transforms_by_height_and_branch, v_3d, static_cast<float>(t), branch_indices_by_height);
}

glm::dvec3 StrandTree::transformToObjectSpaceForBranch(glm::dvec3 v_3d, double t, size_t branch_id) const
{
  if (transforms_by_height_and_branch.empty())
  {
    throw std::runtime_error("StrandTree::transformToObjectSpaceForBranch: no profile transforms");
  }

  std::vector<size_t> constant_branch(transforms_by_height_and_branch.size(), branch_id);
  return transformToObjectSpace(std::move(v_3d), t, constant_branch);
}

Trajectory<2> StrandTree::getPiecePolynomial(size_t strand_id, size_t index, size_t reference_branch) const
{
  if (strand_id >= support_points.size())
  {
    throw std::out_of_range("Strand id " + std::to_string(strand_id) + " out of range.");
  }

  if (index >= support_points[strand_id].size() - 1)
  {
    throw std::out_of_range("Index " + std::to_string(index) + " out of range for piece polynomial.");
  }

  const glm::dvec2 P0 = getPointTransformedAtSection(strand_id, index, reference_branch, index);
  const glm::dvec2 P1 = getPointTransformedAtSection(strand_id, index + 1, reference_branch, index + 1);

  Trajectory<2> result;
  for (int i = 0; i < 2; ++i)
  {
    result[i] = POLYNOMIAL(P0[i] + (P1[i] - P0[i]) * x);
  }

  return result;
}

Trajectory<2> StrandTree::getPiecePolynomialBlendingReference(size_t strand_id, size_t index,
  size_t start_reference_branch, size_t end_reference_branch) const
{
  if (strand_id >= support_points.size())
  {
    throw std::out_of_range("Strand id " + std::to_string(strand_id) + " out of range.");
  }

  if (index >= support_points[strand_id].size() - 1)
  {
    throw std::out_of_range("Index " + std::to_string(index) + " out of range for piece polynomial.");
  }

  const glm::dvec2 P0
    = getPointTransformedAtSection(strand_id, index, start_reference_branch, index);
  const glm::dvec2 P1
    = getPointTransformedAtSection(strand_id, index + 1, end_reference_branch, index + 1);

  Trajectory<2> result;
  for (int i = 0; i < 2; ++i)
  {
    result[i] = POLYNOMIAL(P0[i] + (P1[i] - P0[i]) * x);
  }

  return result;
}

Trajectory<2> StrandTree::getLocalPiecePolynomial(size_t strand_id, size_t index) const
{
  if (strand_id >= support_points.size())
  {
    throw std::out_of_range("Strand id " + std::to_string(strand_id) + " out of range.");
  }

  if (index >= support_points[strand_id].size() - 1)
  {
    throw std::out_of_range("Index " + std::to_string(index) + " out of range for piece polynomial.");
  }

  const glm::dvec2& p0 = support_points[strand_id][index];
  const glm::dvec2& p1 = support_points[strand_id][index + 1];

  Trajectory<2> result;
  for (int i = 0; i < 2; ++i)
  {
    result[i] = POLYNOMIAL(p0[i] + (p1[i] - p0[i]) * x);
  }

  return result;
}

namespace
{
const char* FORMAT_MAGIC = "StrandTree 1";
const char* SECTION_PREFIX = "[Section] ";

void writeComment(std::ostream& out, const char* text) { out << "# " << text << '\n'; }

void writeSectionHeader(std::ostream& out, const char* name) { out << SECTION_PREFIX << name << '\n'; }

std::string trim(std::string s)
{
  const auto not_space = [](unsigned char c) { return !std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
  s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
  return s;
}

std::string stripInlineComment(std::string line)
{
  const size_t hash = line.find('#');
  if (hash != std::string::npos)
  {
    line.resize(hash);
  }
  return trim(std::move(line));
}

bool isCommentOrBlankLine(const std::string& line) { return stripInlineComment(line).empty(); }

bool parseSectionHeader(const std::string& line, const char* expected_name)
{
  std::string content = stripInlineComment(line);
  const std::string prefix = SECTION_PREFIX;
  if (content.rfind(prefix, 0) == 0)
  {
    content = trim(content.substr(prefix.size()));
  }
  return content == expected_name;
}

void writeVec2(std::ostream& out, const glm::dvec2& v) { out << v.x << ' ' << v.y << '\n'; }

void readVec2(std::istream& in, glm::dvec2& v) { in >> v.x >> v.y; }

void writeDMat4(std::ostream& out, const glm::dmat4& m)
{
  // One space between every coefficient. Do not use (col ? " " : "") — that omits the space
  // between the last column of one row and the first column of the next, merging tokens on read.
  const auto flags = out.flags();
  const auto precision = out.precision();
  out << std::scientific << std::setprecision(17);
  for (int row = 0; row < 4; ++row)
  {
    for (int col = 0; col < 4; ++col)
    {
      if (row != 0 || col != 0)
      {
        out << ' ';
      }
      out << m[col][row];
    }
  }
  out.flags(flags);
  out.precision(precision);
  out << '\n';
}

void readDMat4(std::istream& in, glm::dmat4& m)
{
  for (int row = 0; row < 4; ++row)
  {
    for (int col = 0; col < 4; ++col)
    {
      if (!(in >> m[col][row]))
      {
        throw std::runtime_error("StrandTree::loadFromFile: transform matrix has too few values");
      }
    }
  }
  std::string trailing;
  if (in >> trailing)
  {
    throw std::runtime_error("StrandTree::loadFromFile: transform matrix has too many values");
  }
}
} // namespace

void StrandTree::saveToFile(const std::filesystem::path& path) const
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error(
      std::string("StrandTree::saveToFile: cannot open ") + path.string() + std::string(" for writing"));

  writeComment(out, "kinDS StrandTree text format. Lines starting with # are comments.");
  writeComment(out, "Section headers use [Section] <name>. Inline comments after # are ignored on load.");
  out << FORMAT_MAGIC << '\n';
  writeComment(out, "Stored tree height (max support polyline length - 1).");
  out << tree_height << '\n';
  writeComment(out, "Number of strands (moving sites).");
  out << support_points.size() << '\n';

  writeSectionHeader(out, "support_points");
  writeComment(out, "Per strand: point count, then one (x y) profile position per height index.");
  for (const auto& strand : support_points)
  {
    out << strand.size() << '\n';
    for (const auto& p : strand)
      writeVec2(out, p);
  }

  writeSectionHeader(out, "subdivisions_by_strand");
  writeComment(out, "Per strand: count, then subdivision parameters in (0, 1) along each segment.");
  for (const auto& strand : subdivisions_by_strand)
  {
    out << strand.size();
    for (double d : strand)
      out << ' ' << d;
    out << '\n';
  }

  writeSectionHeader(out, "physics_strand_to_segment_indices");
  writeComment(out, "Per strand: count, then physics segment index per height (or -1 if unused).");
  for (const auto& strand : physics_strand_to_segment_indices)
  {
    out << strand.size();
    for (int i : strand)
      out << ' ' << i;
    out << '\n';
  }

  writeSectionHeader(out, "transforms_by_height_and_branch");
  writeComment(out, "Profile-to-world transforms: height count, then per height branch count and 4x4 matrices.");
  writeComment(out, "Matrix layout (row-major file order): m[col][row] for row,col in 0..3; local (u,0,v) axes.");
  out << transforms_by_height_and_branch.size() << '\n';
  for (const auto& by_branch : transforms_by_height_and_branch)
  {
    out << by_branch.size() << '\n';
    for (const auto& m : by_branch)
      writeDMat4(out, m);
  }

  writeSectionHeader(out, "branch_indices");
  writeComment(out, "Per strand: count, then data-structure branch id at each height index.");
  for (const auto& strand : branch_indices)
  {
    out << strand.size();
    for (size_t s : strand)
      out << ' ' << s;
    out << '\n';
  }

  writeSectionHeader(out, "strands_by_branch_id");
  writeComment(out, "Inverse branch membership: height count, branch count per height, strand ids per branch.");
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
  auto getContentLine = [&]() -> bool
  {
    while (getLine())
    {
      if (!isCommentOrBlankLine(line))
      {
        return true;
      }
    }
    return false;
  };
  auto expectSection = [&](const char* tag)
  {
    if (!getContentLine())
    {
      throw std::runtime_error(
        std::string("StrandTree::loadFromFile: unexpected EOF, expected section '") + tag + "'");
    }
    if (!parseSectionHeader(line, tag))
    {
      throw std::runtime_error(std::string("StrandTree::loadFromFile: expected section '") + tag + "', got '" + line
        + "'");
    }
  };
  auto parseContentLine = [&](const char* context) -> std::istringstream
  {
    if (!getContentLine())
    {
      throw std::runtime_error(
        std::string("StrandTree::loadFromFile: unexpected EOF (") + context + ")");
    }
    return std::istringstream(stripInlineComment(line));
  };

  if (!getContentLine() || stripInlineComment(line) != FORMAT_MAGIC)
  {
    throw std::runtime_error("StrandTree::loadFromFile: invalid or unsupported format");
  }

  size_t read_height = 0;
  parseContentLine("tree height") >> read_height;
  size_t num_strands = 0;
  parseContentLine("strand count") >> num_strands;

  std::vector<std::vector<glm::dvec2>> support_points;
  support_points.reserve(num_strands);
  expectSection("support_points");
  for (size_t s = 0; s < num_strands; ++s)
  {
    size_t n = 0;
    parseContentLine("support_points count") >> n;
    std::vector<glm::dvec2> strand;
    strand.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      glm::dvec2 p;
      std::istringstream iss = parseContentLine("support point");
      readVec2(iss, p);
      strand.push_back(p);
    }
    support_points.push_back(std::move(strand));
  }

  std::vector<std::vector<double>> subdivisions_by_strand;
  subdivisions_by_strand.reserve(num_strands);
  expectSection("subdivisions_by_strand");
  for (size_t s = 0; s < num_strands; ++s)
  {
    std::istringstream iss = parseContentLine("subdivisions");
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
  expectSection("physics_strand_to_segment_indices");
  for (size_t s = 0; s < num_strands; ++s)
  {
    std::istringstream iss = parseContentLine("physics indices");
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

  std::vector<std::vector<glm::dmat4>> transforms_by_height_and_branch;
  expectSection("transforms_by_height_and_branch");
  size_t num_heights = 0;
  parseContentLine("transform height count") >> num_heights;
  transforms_by_height_and_branch.reserve(num_heights);
  for (size_t h = 0; h < num_heights; ++h)
  {
    size_t num_branches = 0;
    parseContentLine("transforms branch count") >> num_branches;
    std::vector<glm::dmat4> by_branch;
    by_branch.reserve(num_branches);
    for (size_t b = 0; b < num_branches; ++b)
    {
      glm::dmat4 m;
      std::istringstream iss = parseContentLine("transform matrix");
      readDMat4(iss, m);
      by_branch.push_back(m);
    }
    transforms_by_height_and_branch.push_back(std::move(by_branch));
  }

  std::vector<std::vector<size_t>> branch_indices;
  branch_indices.reserve(num_strands);
  expectSection("branch_indices");
  for (size_t s = 0; s < num_strands; ++s)
  {
    std::istringstream iss = parseContentLine("branch indices");
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
  expectSection("strands_by_branch_id");
  num_heights = 0;
  parseContentLine("strands_by_branch_id height count") >> num_heights;
  strands_by_branch_id.reserve(num_heights);
  for (size_t h = 0; h < num_heights; ++h)
  {
    size_t num_branches = 0;
    parseContentLine("strands_by_branch branch count") >> num_branches;
    std::vector<std::vector<size_t>> by_branch;
    by_branch.reserve(num_branches);
    for (size_t b = 0; b < num_branches; ++b)
    {
      std::istringstream iss = parseContentLine("strands list");
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

void StrandTree::conditionProfileTransforms()
{
  for (size_t height = 0; height < transforms_by_height_and_branch.size(); ++height)
  {
    auto& transforms_at_height = transforms_by_height_and_branch[height];
    for (size_t branch_id = 0; branch_id < transforms_at_height.size(); ++branch_id)
    {
      glm::dmat4& transform = transforms_at_height[branch_id];
      const glm::dvec3 u(transform[0]);
      const glm::dvec3 normal(transform[1]);
      const glm::dvec3 v(transform[2]);

      const double lu = glm::length(u);
      const double lv = glm::length(v);
      if (lu <= 0.0 || lv <= 0.0)
      {
        const bool live = height < strands_by_branch_id.size() && branch_id < strands_by_branch_id[height].size()
          && !strands_by_branch_id[height][branch_id].empty();
        if (live)
        {
          throw std::runtime_error("StrandTree::conditionProfileTransforms: zero profile axes at height "
            + std::to_string(height) + " branch " + std::to_string(branch_id));
        }
        continue;
      }

      // Symmetric (minimal-rotation) orthonormalization of the in-plane axes: rotate u and v by equal and opposite
      // amounts about their bisector until they are orthogonal, then give both the geometric-mean length so the
      // in-plane block becomes scalar * orthonormal.
      const glm::dvec3 a = u / lu;
      const glm::dvec3 b = v / lv;
      const glm::dvec3 sum = a + b; // bisector direction
      const glm::dvec3 diff = a - b; // in-plane, perpendicular to the bisector (|a| == |b|)
      const double sum_len = glm::length(sum);
      const double diff_len = glm::length(diff);
      if (sum_len <= 1e-12 || diff_len <= 1e-12)
      {
        continue;
      }

      const glm::dvec3 bisector = sum / sum_len;
      const glm::dvec3 perp = diff / diff_len;
      const double scale = std::sqrt(lu * lv); // representative in-plane scale, preserved per branch

      // e_u = (bisector + perp)/sqrt(2) stays near a, e_v = (bisector - perp)/sqrt(2) stays near b; both orthonormal.
      const glm::dvec3 new_u = scale * glm::normalize(bisector + perp);
      const glm::dvec3 new_v = scale * glm::normalize(bisector - perp);

      // Give the normal the same length as the in-plane axes so the whole 3x3 is scalar * orthogonal. The normal does
      // not affect point positions (it is zeroed for points), but this keeps normal transforms consistent too.
      const double ln = glm::length(normal);
      const glm::dvec3 unit_normal = (ln > 0.0) ? normal / ln : glm::normalize(glm::cross(new_u, new_v));
      const glm::dvec3 new_normal = scale * unit_normal;

      transform[0] = glm::dvec4(new_u, transform[0].w);
      transform[1] = glm::dvec4(new_normal, transform[1].w);
      transform[2] = glm::dvec4(new_v, transform[2].w);
    }
  }
}

void StrandTree::computeNormalTransforms()
{
  normal_transforms_by_height_and_branch.resize(getTransformsByHeightAndBranch().size());

  for (size_t i = 0; i < getTransformsByHeightAndBranch().size(); i++)
  {
    normal_transforms_by_height_and_branch[i].resize(getTransformsByHeightAndBranch()[i].size());
    for (size_t j = 0; j < normal_transforms_by_height_and_branch[i].size(); j++)
    {
      normal_transforms_by_height_and_branch[i][j]
        = glm::transpose(glm::inverse(getTransformsByHeightAndBranch()[i][j]));
      normal_transforms_by_height_and_branch[i][j][3] = glm::dvec4(0.0f, 0.0f, 0.0f, 1.0f);
    }
  }
}
