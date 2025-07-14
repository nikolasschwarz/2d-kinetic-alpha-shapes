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

class VoronoiSiteTrajectory : public std::vector<std::array<RationalFunction, 2>>
{
 public:
  double start_time; // Start time of the trajectory
  double end_time; // End time of the trajectory
};

class RuledSurface
{
 private:
  std::vector<VoronoiSiteTrajectory> left_trajectory; // piecewise trajectory for the left side of the ruled surface
  std::vector<double> left_bounds; // bounds for the left trajectory
  VoronoiSiteTrajectory right_trajectory; // piecewise trajectory for the right side of the ruled surface
  VoronoiSiteTrajectory right_bounds; // bounds for the right trajectory
  double lower_bound; // Lower bound of the ruled surface
  double upper_bound; // Upper bound of the ruled surface

 public:
  RuledSurface(const std::vector<VoronoiSiteTrajectory>& left, const std::vector<double>& left_bounds,
    const VoronoiSiteTrajectory& right, const VoronoiSiteTrajectory& right_bounds,
    double lower, double upper)
    : left_trajectory(left)
    , left_bounds(left_bounds)
    , right_trajectory(right)
    , right_bounds(right_bounds)
    , lower_bound(lower)
    , upper_bound(upper)
  {
    if (left_trajectory.size() != left_bounds.size() - 1)
    {
      throw std::invalid_argument("Left trajectory and bounds size mismatch.");
    }

    if (right_trajectory.size() != right_bounds.size() - 1)
    {
      throw std::invalid_argument("Right trajectory and bounds size mismatch.");
    }
  }

  double getLowerBound() const { return lower_bound; }
  double getUpperBound() const { return upper_bound; }
  std::pair<Point<2>, Point<2>> operator()(double t)
  {
    if (t < lower_bound || t > upper_bound)
    {
      throw std::out_of_range("Time t is out of bounds of the ruled surface.");
    }
    // Evaluate the left and right trajectories at time t
    int left_index = std::lower_bound(left_bounds.begin(), left_bounds.end(), t) - left_bounds.begin() - 1;
    int right_index = std::lower_bound(right_bounds.begin(), right_bounds.end(), t) - right_bounds.begin() - 1;

    /*return {
      { left_trajectory[left_index](t), left_trajectory[left_index].getSiteTrajectory()[1](t) },
      { right_trajectory.getSiteTrajectory()[0](t), right_trajectory.getSiteTrajectory()[1](t) }
    };*/
  }
};
} // namespace kinDS
