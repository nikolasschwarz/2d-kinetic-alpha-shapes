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

 public:
  RuledSurface()
    : is_initialized(false)
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
};
} // namespace kinDS
