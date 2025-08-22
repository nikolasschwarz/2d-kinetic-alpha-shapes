#pragma once

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
}
