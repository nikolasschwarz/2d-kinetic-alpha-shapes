#pragma once
#include "Polynomial.hpp"

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
        if (denominator.coeffs.size() == 0 || (denominator.coeffs.size() == 1 && denominator.coeffs[0] == 0))
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

class VoronoiSiteTrajectory
{
private:
    std::vector<std::array<RationalFunction, 2>> site_trajectories; // Each site trajectory is a pair of rational functions (x(t), y(t))
    double start_time; // Start time of the trajectory
    double end_time; // End time of the trajectory
};

class RuledSurface
{
};
