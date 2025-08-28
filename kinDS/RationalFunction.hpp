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
  RationalFunction(const Polynomial& num, const Polynomial& den);

  // Evaluate the rational function at a given point
  double operator()(double x) const;

  void reduce();

  // Now implement operators
  RationalFunction operator+(const RationalFunction& other) const;

  RationalFunction operator-() const;

  RationalFunction operator-(const RationalFunction& other) const;

  RationalFunction operator*(const RationalFunction& other) const;

  RationalFunction inv() const;

  RationalFunction operator/(const RationalFunction& other) const;

  // Get the numerator polynomial
  const Polynomial& getNumerator() const;
  // Get the denominator polynomial
  const Polynomial& getDenominator() const;
};

RationalFunction operator/(const Polynomial& num, const Polynomial& den);
}
