#pragma once
#include "Polynomial.hpp"
#include "RationalFunction.hpp"

namespace kinDS
{

class VoronoiSiteTrajectory : public std::array<RationalFunction, 2>
{
 public:
  // double start_time; // Start time of the trajectory
  // double end_time; // End time of the trajectory
  static VoronoiSiteTrajectory circumcenterTrajectory(const Trajectory<2>& a, const Trajectory<2>& b, const Trajectory<2>& c)
  {
    // Calculate the circumcenter of the triangle formed by points a, b, c
    Polynomial D = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    RationalFunction Ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / D;
    RationalFunction Uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / D;
    return VoronoiSiteTrajectory { Ux, Uy };
  }
};

class RuledSurface
{
 private:
  std::vector<VoronoiSiteTrajectory> left_trajectory; // piecewise trajectory for the left side of the ruled surface
  std::vector<double> left_bounds; // bounds for the left trajectory
  std::vector<VoronoiSiteTrajectory> right_trajectory; // piecewise trajectory for the right side of the ruled surface
  std::vector<double> right_bounds; // bounds for the right trajectory
  bool is_initialized = false; // Flag to check if the ruled surface is initialized
  // use faces to establish an absolute order at vertex extraction
  size_t left_face_index; // Index of the left face in the delaunay triangulation
  size_t right_face_index; // Index of the right face in the delaunay triangulation
  std::vector<std::pair<double, Point<2>>> left_event_points; // Left event points by time
  std::vector<std::pair<double, Point<2>>> right_event_points; // Right event points by time
  std::vector<double> strand_subdivision; // Subdivision points along the strand, used for mesh extraction

 public:
  RuledSurface(size_t left_face_index, size_t right_face_index)
    : is_initialized(false)
    , left_face_index(left_face_index)
    , right_face_index(right_face_index)
  {
  }

  double lowerBound() const
  {
    assert(is_initialized && "Ruled surface must be initialized before accessing bounds.");
    return left_bounds[0];
  }
  double upperBound() const
  {
    assert(is_initialized && "Ruled surface must be initialized before accessing bounds.");
    if (left_bounds.back() != right_bounds.back())
    {
      logger.log(WARNING, "Left and right bounds of the ruled surface do not match. Returning the maximum of both.");
    }

    return std::max(left_bounds.back(), right_bounds.back());
  }
  std::pair<Point<2>, Point<2>> operator()(double t)
  {
    if (t < lowerBound() || t > left_bounds.back() || t > right_bounds.back())
    {
      throw std::out_of_range("Time t is out of bounds of the ruled surface.");
    }
    // Evaluate the left and right trajectories at time t
    int left_index = std::lower_bound(left_bounds.begin(), left_bounds.end(), t) - left_bounds.begin() - 1;
    int right_index = std::lower_bound(right_bounds.begin(), right_bounds.end(), t) - right_bounds.begin() - 1;

    return {
      Point<2> { left_trajectory[left_index][0](t), left_trajectory[left_index][1](t) },
      Point<2> { right_trajectory[right_index][0](t), right_trajectory[right_index][1](t) }
    };
  }

  void insertLeft(const VoronoiSiteTrajectory& traj, double lb)
  {
    // only allow whole numbers, else the event function should be used
    assert(std::floor(lb) == lb && "Lower bound must be a whole number.");
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > left_bounds.back() && "Lower bound must be greater than the last left bound.");
    left_trajectory.push_back(traj);
    left_bounds.push_back(lb);
  }

  void insertLeftEvent(const Point<2>& p, const VoronoiSiteTrajectory& traj, double lb)
  {
    assert(is_initialized && "Ruled surface must be initialized before inserting event points.");
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > left_bounds.back() && "Lower bound must be greater than the last left bound.");
    left_trajectory.push_back(traj);
    left_bounds.push_back(lb);
    left_event_points.emplace_back(lb, p);
  }

  void insertRight(const VoronoiSiteTrajectory& traj, double lb)
  {
    // only allow whole numbers, else the event function should be used
    assert(std::floor(lb) == lb && "Lower bound must be a whole number.");
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > right_bounds.back() && "Lower bound must be greater than the last right bound.");
    right_trajectory.push_back(traj);
    right_bounds.push_back(lb);
  }

  void insertRightEvent(const Point<2>& p, const VoronoiSiteTrajectory& traj, double lb)
  {
    assert(is_initialized && "Ruled surface must be initialized before inserting event points.");
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > left_bounds.back() && "Lower bound must be greater than the last left bound.");
    right_trajectory.push_back(traj);
    right_bounds.push_back(lb);
    right_event_points.emplace_back(lb, p);
  }

  void init(const VoronoiSiteTrajectory& left_traj, VoronoiSiteTrajectory& right_traj, double lb)
  {
    left_trajectory.push_back(left_traj);
    right_trajectory.push_back(right_traj);
    left_bounds.push_back(lb);
    right_bounds.push_back(lb);
    is_initialized = true;
  }

  void initEvent(const Point<2>& p, const VoronoiSiteTrajectory& left_traj, VoronoiSiteTrajectory& right_traj, double lb)
  {
    left_trajectory.push_back(left_traj);
    right_trajectory.push_back(right_traj);
    left_bounds.push_back(lb);
    right_bounds.push_back(lb);
    left_event_points.emplace_back(lb, p);
    right_event_points.emplace_back(lb, p);
    is_initialized = true;
  }

  void finalize(double upper_bound)
  {
    assert(!isFinalized() && "Ruled surface is already finalized.");
    assert(is_initialized && "Ruled surface must be initialized before finalizing.");
    left_bounds.push_back(upper_bound);
    right_bounds.push_back(upper_bound);
  }

  void finalizeEvent(const Point<2>& p, double upper_bound)
  {
    assert(!isFinalized() && "Ruled surface is already finalized.");
    assert(is_initialized && "Ruled surface must be initialized before finalizing.");
    left_bounds.push_back(upper_bound);
    right_bounds.push_back(upper_bound);
    left_event_points.emplace_back(upper_bound, p);
    right_event_points.emplace_back(upper_bound, p);
  }

  bool isFinalized() const
  {
    // Note: The other conditions (right side and initialization) are implied
    return left_bounds.size() == left_trajectory.size() + 1;
  }

  std::vector<Mesh> extractTriangles(bool invert, std::vector<std::vector<double>>& strand_subdivisions, size_t strand_id) const
  {
    // TODO: take inversion into account
    assert(isFinalized() && "Ruled surface must be finalized before extracting triangles.");

    std::vector<Mesh> meshes;
    std::vector<Point<3>> vertices;
    std::vector<size_t> indices;

    double t = left_bounds[0]; // Start time for the triangles
    double upper_bound = upperBound();

    size_t left_event_index = 0;
    size_t right_event_index = 0;

    // Initialize the first vertices
    size_t left_vertex_index = 0;
    // distinguish between event points and regular points
    if (!left_event_points.empty() && left_event_points[0].first == t)
    {
      vertices.emplace_back(Point<3> { left_event_points[0].second[0], left_event_points[0].second[1], t });
      left_event_index++;
    }
    else
    {
      vertices.emplace_back(Point<3> { left_trajectory[0][0](t - std::floor(t)), left_trajectory[0][1](t - std::floor(t)), t });
    }

    size_t right_vertex_index = 0;

    if (!right_event_points.empty() && right_event_points[0].first == t)
    {
      vertices.emplace_back(Point<3> { right_event_points[0].second[0], right_event_points[0].second[1], t });
      right_event_index++;
    }
    else
    {
      vertices.emplace_back(Point<3> { right_trajectory[0][0](t - std::floor(t)), right_trajectory[0][1](t - std::floor(t)), t });
    }

    size_t left_index = 0;
    size_t right_index = 0;
    size_t subdivision_index = 0;

    // TODO: whenever t exceeds a strand subdivision, create a new mesh

    while (t <= upper_bound)
    {
      if (left_index >= left_trajectory.size() && right_index >= right_trajectory.size())
      {
        break; // No more left trajectories to process
      }

      // Decide whether we have to subdivide here
      std::function<bool()> subdivide = [&]()
      {
        if (strand_subdivision.empty())
          return false;

        if (subdivision_index >= strand_subdivision.size())
          return false;

        if (left_index + 1 >= left_bounds.size() || right_index + 1 >= right_bounds.size())
          return false; // No more bounds to compare

        return (strand_subdivision[subdivision_index] < left_bounds[left_index + 1]) && (strand_subdivision[subdivision_index] < right_bounds[right_index + 1]);
      };

      if (subdivide())
      {
        // Finish the current mesh and start a new one
        t = strand_subdivision[subdivision_index];
        subdivision_index++;

        // Add the two top points at time t
        Point<3> left_point { left_trajectory[left_index][0](t - std::floor(t)), left_trajectory[left_index][1](t - std::floor(t)), t };
        vertices.emplace_back(left_point);
        size_t prev_left_vertex_index = left_vertex_index;
        left_vertex_index = vertices.size() - 1;

        Point<3> right_point { right_trajectory[right_index][0](t - std::floor(t)), right_trajectory[right_index][1](t - std::floor(t)), t };
        vertices.emplace_back(right_point);
        size_t prev_right_vertex_index = right_vertex_index;
        right_vertex_index = vertices.size() - 1;

        // Put triangles into index buffer
        indices.push_back(prev_left_vertex_index);
        indices.push_back(prev_right_vertex_index);
        indices.push_back(left_vertex_index);

        indices.push_back(prev_right_vertex_index);
        indices.push_back(right_vertex_index);
        indices.push_back(left_vertex_index);

        meshes.push_back(Mesh(vertices, indices));

        vertices.clear();
        indices.clear();

        left_vertex_index = 0;
        right_vertex_index = 0;
        // Re-insert the current vertices at time t
        vertices.emplace_back(left_point);
        left_vertex_index = 0;
        vertices.emplace_back(right_point);
        right_vertex_index = 1;
        continue; // Re-evaluate subdivision after resetting
      }

      std::function<bool()> advance_left = [&]()
      {
        if (left_bounds[left_index] < right_bounds[right_index])
          return true; // Left trajectory has a smaller bound

        if (left_bounds[left_index] > right_bounds[right_index])
        {
          return false;
        }

        // if both are equal, first check if one index is always at maximum and then return the other side
        if (left_index == left_bounds.size() - 1)
        {
          return false; // Left trajectory is at its maximum, prefer right
        }
        else if (right_index == right_bounds.size() - 1)
        {
          return true; // Right trajectory is at its maximum, prefer left
        }

        // If bounds are equal, prefer the side that is closer
        return left_bounds[left_index + 1] < right_bounds[right_index + 1] || (left_bounds[left_index + 1] == right_bounds[right_index + 1] && left_face_index < right_face_index);
      };

      if (advance_left())
      {

        left_index++;

        // Insert new vertices for the left trajectory
        t = left_bounds[left_index];

        // test if t also exist in the event points
        while (left_event_index < left_event_points.size() && left_event_points[left_event_index].first < t)
        {
          left_event_index++;
        }
        if (left_event_index < left_event_points.size() && left_event_points[left_event_index].first == t)
        {
          vertices.emplace_back(Point<3> { left_event_points[left_event_index].second[0], left_event_points[left_event_index].second[1], t });
          left_event_index++;
        }
        else
        {
          // assert that t is a whole number
          assert(std::floor(t) == t && "Left bound must be a whole number if it is not an event point.");

          // Last index must be treated differently because it can be at frac = 1.0 in the last trajectory piece
          if (left_index == left_trajectory.size())
          {
            double frac = t - std::floor(t);
            if (frac == 0.0)
              frac = 1.0;
            vertices.emplace_back(Point<3> { left_trajectory.back()[0](frac), left_trajectory.back()[1](frac), t });
          }
          else
          {
            vertices.emplace_back(Point<3> { left_trajectory[left_index][0](t - std::floor(t)), left_trajectory[left_index][1](t - std::floor(t)), t });
          }
        }

        // Create triangles with the previous right vertex
        indices.push_back(left_vertex_index);
        indices.push_back(right_vertex_index);

        left_vertex_index = vertices.size() - 1;
        indices.push_back(left_vertex_index); // New left vertex index
      }
      else
      {
        right_index++;

        // Insert new vertices for the right trajectory
        t = right_bounds[right_index];

        // test if t also exist in the event points
        while (right_event_index < right_event_points.size() && right_event_points[right_event_index].first < t)
        {
          right_event_index++;
        }
        if (right_event_index < right_event_points.size() && right_event_points[right_event_index].first == t)
        {
          vertices.emplace_back(Point<3> { right_event_points[right_event_index].second[0], right_event_points[right_event_index].second[1], t });
          right_event_index++;
        }
        else
        {
          // assert that t is a whole number
          assert(std::floor(t) == t && "Right bound must be a whole number if it is not an event point.");
          // Last index must be treated differently because it can be at frac = 1.0 in the last trajectory piece
          if (right_index == right_trajectory.size())
          {
            double frac = t - std::floor(t);
            if (frac == 0.0)
              frac = 1.0;
            vertices.emplace_back(Point<3> { right_trajectory.back()[0](frac), right_trajectory.back()[1](frac), t });
          }
          else
          {
            vertices.emplace_back(Point<3> { right_trajectory[right_index][0](t - std::floor(t)), right_trajectory[right_index][1](t - std::floor(t)), t });
          }
        }

        // Create triangles with the previous left vertex
        indices.push_back(left_vertex_index);
        indices.push_back(right_vertex_index);

        right_vertex_index = vertices.size() - 1;
        indices.push_back(right_vertex_index); // New right vertex index
      }
    }

    meshes.push_back(Mesh(vertices, indices));
    return meshes;
  }

  void printDebugInfo(std::string prefix = "") const
  {
    logger.log(INFO, "%sRuledSurface Debug Info:", prefix.c_str());
    logger.log(INFO, "%sLeft Trajectory Count: %zu", prefix.c_str(), left_trajectory.size());
    logger.log(INFO, "%sRight Trajectory Count: %zu", prefix.c_str(), right_trajectory.size());
    logger.log(INFO, "%sLeft Bounds Count: %zu", prefix.c_str(), left_bounds.size());

    // print all left bounds
    for (const auto& lb : left_bounds)
    {
      logger.log(INFO, "%sLeft Bound: %f", prefix.c_str(), lb);
    }

    logger.log(INFO, "%sRight Bounds Count: %zu", prefix.c_str(), right_bounds.size());

    // print all right bounds
    for (const auto& rb : right_bounds)
    {
      logger.log(INFO, "%sRight Bound: %f", prefix.c_str(), rb);
    }

    logger.log(INFO, "%sIs Initialized: %s", prefix.c_str(), is_initialized ? "true" : "false");
    logger.log(INFO, "%sIs Finalized: %s", prefix.c_str(), isFinalized() ? "true" : "false");
  }

  void debugCompareSamplesAtKnots()
  {
    logger.log(INFO, "Debugging RuledSurface samples at knots:");
    for (size_t i = 1; i < left_bounds.size() - 1; ++i)
    {
      double t_lower = left_bounds[i] - std::floor(left_bounds[i]);
      double t_upper = t_lower;
      if (t_upper == 0.0)
      {
        t_upper = 1.0; // Handle the case where the upper bound is exactly an integer
      }

      auto left_upper_sample = Point<2> {
        left_trajectory[i - 1][0](t_upper), left_trajectory[i - 1][1](t_upper)
      };
      auto left_lower_sample = Point<2> {
        left_trajectory[i][0](t_lower), left_trajectory[i][1](t_lower)
      };

      // now compare and print the results
      if (left_upper_sample.dist(left_lower_sample) > std::numeric_limits<double>::epsilon())
      {
        logger.log(INFO, "Upper Left Trajectory Sample at t = %f: %s", left_bounds[i], left_upper_sample.toString().c_str());
        logger.log(INFO, "Lower Left Trajectory Sample at t = %f: %s", left_bounds[i], left_lower_sample.toString().c_str());
      }
    }

    // same for right trajectory
    for (size_t i = 1; i < right_bounds.size() - 1; ++i)
    {
      double t_lower = right_bounds[i] - std::floor(right_bounds[i]);
      double t_upper = t_lower;
      if (t_upper == 0.0)
      {
        t_upper = 1.0; // Handle the case where the upper bound is exactly an integer
      }

      auto right_upper_sample = Point<2> {
        right_trajectory[i - 1][0](t_upper), right_trajectory[i - 1][1](t_upper)
      };
      auto right_lower_sample = Point<2> {
        right_trajectory[i][0](t_lower), right_trajectory[i][1](t_lower)
      };

      // now compare and print the results
      if (right_upper_sample.dist(right_lower_sample) > std::numeric_limits<double>::epsilon())
      {
        logger.log(INFO, "Upper Right Trajectory Sample at t = %f: %s", right_bounds[i], right_upper_sample.toString().c_str());
        logger.log(INFO, "Lower Right Trajectory Sample at t = %f: %s", right_bounds[i], right_lower_sample.toString().c_str());
      }
    }
  }

  void applySubdivision(const std::vector<double>& subdivision)
  {
    // if (subdivision.empty())
    //   return; // No subdivision to apply

    // Assert that the subdivision points are sorted
    assert(std::is_sorted(subdivision.begin(), subdivision.end()) && "Subdivision points must be sorted.");

    // Insert subdivision points into the ruled surface
    // We acheive this by merging the subdivision points with the existing bounds into a new vector and then copying it back
    // Anything outside the bounds of the ruled surface is ignored
    /* std::vector<double> new_left_bounds;
    std::vector<double> new_right_bounds;
    new_left_bounds.reserve(left_bounds.size() + subdivision.size());
    new_right_bounds.reserve(right_bounds.size() + subdivision.size());
    size_t left_index = 0;
    size_t right_index = 0;
    size_t subdiv_index = 0;
    double lower_bound = lowerBound();
    double upper_bound = upperBound();

    while (left_index < left_bounds.size() || right_index < right_bounds.size() || subdiv_index < subdivision.size())
    {
      double left_val = left_index < left_bounds.size() ? left_bounds[left_index] : std::numeric_limits<double>::infinity();
      double right_val = right_index < right_bounds.size() ? right_bounds[right_index] : std::numeric_limits<double>::infinity();
      double subdiv_val = subdiv_index < subdivision.size() ? subdivision[subdiv_index] : std::numeric_limits<double>::infinity();
      double next_val = std::min({ left_val, right_val, subdiv_val });
      if (next_val == std::numeric_limits<double>::infinity())
        break; // All vectors are exhausted

      // if subdiv_val is the smallest, insert it into both
      if (subdiv_val == next_val)
      {
        if (subdiv_val >= lower_bound && subdiv_val <= upper_bound)
        {
          new_left_bounds.push_back(subdiv_val);
          new_right_bounds.push_back(subdiv_val);
        }
      }
      else
      {
        if (left_val == next_val)
        {
          new_left_bounds.push_back(left_val);
        }
        if (right_val == next_val)
        {
          new_right_bounds.push_back(right_val);
        }
      }

      if (next_val == left_val)
        left_index++;
      if (next_val == right_val)
        right_index++;
      if (next_val == subdiv_val)
        subdiv_index++;
    }

    left_bounds = std::move(new_left_bounds);
    right_bounds = std::move(new_right_bounds);*/

    // This above does not work, instead we simply store the subdivision, if needed merge it with the previous one, and take care of it during mesh extraction
    std::vector<double> new_subdivision;
    new_subdivision.reserve(strand_subdivision.size() + subdivision.size());
    size_t i = 0, j = 0;
    while (i < strand_subdivision.size() && j < subdivision.size())
    {
      if (strand_subdivision[i] < subdivision[j])
      {
        new_subdivision.push_back(strand_subdivision[i]);
        i++;
      }
      else if (strand_subdivision[i] > subdivision[j])
      {
        new_subdivision.push_back(subdivision[j]);
        j++;
      }
      else
      {
        new_subdivision.push_back(strand_subdivision[i]); // Both are equal, add one
        i++;
        j++;
      }
    }
    // copy remaining elements
    if (i < strand_subdivision.size())
    {
      new_subdivision.insert(new_subdivision.end(), strand_subdivision.begin() + i, strand_subdivision.end());
    }
    if (j < subdivision.size())
    {
      new_subdivision.insert(new_subdivision.end(), subdivision.begin() + j, subdivision.end());
    }

    strand_subdivision = std::move(new_subdivision);

    assert(std::is_sorted(strand_subdivision.begin(), strand_subdivision.end()) && "Merged subdivisions should be sorted.");
  }
};
} // namespace kinDS
