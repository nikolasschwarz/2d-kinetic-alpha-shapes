#pragma once

#include "Polynomial.hpp"

#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>

namespace kinDS
{

/**
 * \brief A pointwise polynomial trajectory in \R^dim.
 *
 * This is essentially an array of `Polynomial` of length `dim`, but wrapped
 * in a type with convenience operations defined componentwise.
 */
template<std::size_t dim> class Trajectory : public std::array<Polynomial, dim>
{
  using Base = std::array<Polynomial, dim>;

 public:
  using Base::Base; // inherit std::array constructors

  // Default constructor: leaves component polynomials default-constructed.
  Trajectory() = default;

  explicit Trajectory(const Base& arr)
    : Base(arr)
  {
  }

  // Allow brace-initialization from a list of Polynomials, e.g.
  // Trajectory<3> t{p0, p1, p2};
  Trajectory(std::initializer_list<Polynomial> init)
  {
    std::size_t i = 0;
    for (const auto& p : init)
    {
      if (i >= dim)
        break;
      (*this)[i++] = p;
    }
    // Any remaining components (if init.size() < dim) stay default-constructed.
  }

  // Componentwise addition
  Trajectory operator+(const Trajectory& other) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i] + other[i];
    }
    return result;
  }

  // Componentwise subtraction
  Trajectory operator-(const Trajectory& other) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i] - other[i];
    }
    return result;
  }

  // Componentwise multiplication
  Trajectory operator*(const Trajectory& other) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i] * other[i];
    }
    return result;
  }

  // Scalar multiplication (on the right)
  Trajectory operator*(double scalar) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i] * scalar;
    }
    return result;
  }

  // Scalar division
  Trajectory operator/(double scalar) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i] / scalar;
    }
    return result;
  }

  // Compound assignment operators
  Trajectory& operator+=(const Trajectory& other)
  {
    for (std::size_t i = 0; i < dim; ++i)
    {
      (*this)[i] = (*this)[i] + other[i];
    }
    return *this;
  }

  Trajectory& operator-=(const Trajectory& other)
  {
    for (std::size_t i = 0; i < dim; ++i)
    {
      (*this)[i] = (*this)[i] - other[i];
    }
    return *this;
  }

  Trajectory& operator*=(double scalar)
  {
    for (std::size_t i = 0; i < dim; ++i)
    {
      (*this)[i] = (*this)[i] * scalar;
    }
    return *this;
  }

  Trajectory& operator/=(double scalar)
  {
    for (std::size_t i = 0; i < dim; ++i)
    {
      (*this)[i] = (*this)[i] / scalar;
    }
    return *this;
  }

  // Pointwise derivative: returns a new trajectory whose components are the derivatives.
  Trajectory derivative(std::size_t order = 1) const
  {
    Trajectory result;
    for (std::size_t i = 0; i < dim; ++i)
    {
      result[i] = (*this)[i].derivative(order);
    }
    return result;
  }

  // Pointwise dot product: returns a Polynomial sum_i a_i * b_i.
  static Polynomial dot(const Trajectory& a, const Trajectory& b)
  {
    Polynomial result(0.0);
    for (std::size_t i = 0; i < dim; ++i)
    {
      result = result + a[i] * b[i];
    }
    return result;
  }

  Polynomial dot(const Trajectory& other) const { return dot(*this, other); }

  // Pointwise squared norm: dot product with itself.
  Polynomial squaredNorm() const { return dot(*this, *this); }

  // Cross product: only defined for dim == 3.
  template<std::size_t D = dim>
  static std::enable_if_t<D == 3, Trajectory> cross(const Trajectory& a, const Trajectory& b)
  {
    Trajectory result;
    // a x b in homogeneous coordinates (componentwise polynomials)
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
  }
};

// Scalar multiplication (scalar on the left)
template<std::size_t dim> Trajectory<dim> operator*(double scalar, const Trajectory<dim>& traj)
{
  return traj * scalar;
}

// Polynomial * Trajectory: multiply each component polynomial by `poly`.
template<std::size_t dim> Trajectory<dim> operator*(const Polynomial& poly, const Trajectory<dim>& traj)
{
  Trajectory<dim> result;
  for (std::size_t i = 0; i < dim; ++i)
  {
    result[i] = poly * traj[i];
  }
  return result;
}

// Trajectory * Polynomial: multiply each component polynomial by `poly`.
template<std::size_t dim> Trajectory<dim> operator*(const Trajectory<dim>& traj, const Polynomial& poly)
{
  return poly * traj;
}

} // namespace kinDS
