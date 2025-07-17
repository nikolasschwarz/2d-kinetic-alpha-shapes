#pragma once
#include "Polynomial.hpp"

namespace kinDS
{

class RationalFunction
{
 private:
  Polynomial numerator; // The numerator polynomial
  Polynomial denominator; // The denominator polynomial

 public:
  RationalFunction(const Polynomial& num, const Polynomial& den)
    : numerator(num)
    , denominator(den)
  {
    if (denominator.getCoefficients().size() == 0 || (denominator.getCoefficients().size() == 1 && denominator.getCoefficients()[0] == 0))
    {
      throw std::invalid_argument("Denominator cannot be zero.");
    }
  }
  // Evaluate the rational function at a given point
  double operator()(double x) const
  {
    return numerator(x) / denominator(x);
  }

  void reduce()
  {
    Polynomial::reduce(numerator, denominator);
  }

  // Now implement operators
  RationalFunction operator+(const RationalFunction& other) const
  {
    Polynomial num = numerator * other.denominator + other.numerator * denominator;
    Polynomial den = denominator * other.denominator;
    RationalFunction result(num, den);
    result.reduce();
    return result;
  }

  RationalFunction operator-() const
  {
    return RationalFunction(-numerator, denominator);
  }

  RationalFunction operator-(const RationalFunction& other) const
  {
    return *this + (-other); // Use the addition operator for subtraction
  }

  RationalFunction operator*(const RationalFunction& other) const
  {
    Polynomial num = numerator * other.numerator;
    Polynomial den = denominator * other.denominator;
    RationalFunction result(num, den);
    result.reduce();
    return result;
  }

  RationalFunction inv() const
  {
    // Invert the rational function by swapping numerator and denominator
    if (numerator.degree() == -1)
    {
      throw std::invalid_argument("Cannot invert a rational function with zero numerator.");
    }

    RationalFunction result(denominator, numerator);
    return result;
  }

  RationalFunction operator/(const RationalFunction& other) const
  {
    return *this * other.inv(); // Use the multiplication operator for division
  }

  // Get the numerator polynomial
  const Polynomial& getNumerator() const { return numerator; }
  // Get the denominator polynomial
  const Polynomial& getDenominator() const { return denominator; }
};

RationalFunction operator/(const Polynomial& num, const Polynomial& den)
{
  return RationalFunction { num, den };
}

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
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > left_bounds.back() && "Lower bound must be greater than the last left bound.");
    left_trajectory.push_back(traj);
    left_bounds.push_back(lb);
  }

  void insertRight(const VoronoiSiteTrajectory& traj, double lb)
  {
    assert(is_initialized && "Ruled surface must be initialized before inserting trajectories.");
    assert(lb > right_bounds.back() && "Lower bound must be greater than the last right bound.");
    right_trajectory.push_back(traj);
    right_bounds.push_back(lb);
  }

  void init(const VoronoiSiteTrajectory& left_traj, VoronoiSiteTrajectory& right_traj, double lb)
  {
    left_trajectory.push_back(left_traj);
    right_trajectory.push_back(right_traj);
    left_bounds.push_back(lb);
    right_bounds.push_back(lb);
    is_initialized = true;
  }

  void finalize(double upper_bound)
  {
    assert(!isFinalized() && "Ruled surface is already finalized.");
    assert(is_initialized && "Ruled surface must be initialized before finalizing.");
    left_bounds.push_back(upper_bound);
    right_bounds.push_back(upper_bound);
  }

  bool isFinalized() const
  {
    // Note: The other conditions (right side and initialization) are implied
    return left_bounds.size() == left_trajectory.size() + 1;
  }

  void extractTriangles(std::vector<Point<3>>& vertices, std::vector<size_t>& indices) const
  {
    assert(isFinalized() && "Ruled surface must be finalized before extracting triangles.");

    double t = left_bounds[0]; // Start time for the triangles
    double upper_bound = upperBound();

    // Initialize the first vertices
    size_t left_vertex_index = vertices.size();
    vertices.emplace_back(Point<3> { left_trajectory[0][0](t - std::floor(t)), left_trajectory[0][1](t - std::floor(t)), t });
    size_t right_vertex_index = vertices.size();
    vertices.emplace_back(Point<3> { right_trajectory[0][0](t - std::floor(t)), right_trajectory[0][1](t - std::floor(t)), t });

    size_t left_index = 0;
    size_t right_index = 0;

    while (t <= upper_bound)
    {
      if (left_index >= left_trajectory.size() && right_index >= right_trajectory.size())
      {
        break; // No more left trajectories to process
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

        // Create triangles with the previous left vertex
        indices.push_back(left_vertex_index);
        indices.push_back(right_vertex_index);

        right_vertex_index = vertices.size() - 1;
        indices.push_back(right_vertex_index); // New right vertex index
      }
    }
  }

  void printDebugInfo() const
  {
    logger.log(INFO, "RuledSurface Debug Info:");
    logger.log(INFO, "Left Trajectory Count: %zu", left_trajectory.size());
    logger.log(INFO, "Right Trajectory Count: %zu", right_trajectory.size());
    logger.log(INFO, "Left Bounds Count: %zu", left_bounds.size());

    // print all left bounds
    for (const auto& lb : left_bounds)
    {
      logger.log(INFO, "Left Bound: %f", lb);
    }

    logger.log(INFO, "Right Bounds Count: %zu", right_bounds.size());

    // print all right bounds
    for (const auto& rb : right_bounds)
    {
      logger.log(INFO, "Right Bound: %f", rb);
    }

    logger.log(INFO, "Is Initialized: %s", is_initialized ? "true" : "false");
    logger.log(INFO, "Is Finalized: %s", isFinalized() ? "true" : "false");
  }
};
} // namespace kinDS
