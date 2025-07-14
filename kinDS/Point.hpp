#pragma once
#include <array>

namespace kinDS
{
template<size_t dim>
using Point = std::array<double, dim>;

// operators for Point
template<size_t dim>
Point<dim> operator+(const Point<dim>& a, const Point<dim>& b)
{
  Point<dim> result {};
  for (size_t i = 0; i < dim; ++i)
  {
    result[i] = a[i] + b[i];
  }
  return result;
}

template<size_t dim>
Point<dim> operator-(const Point<dim>& a, const Point<dim>& b)
{
  Point<dim> result {};
  for (size_t i = 0; i < dim; ++i)
  {
    result[i] = a[i] - b[i];
  }
  return result;
}

template<size_t dim>
Point<dim> operator*(const Point<dim>& a, double scalar)
{
  Point<dim> result {};
  for (size_t i = 0; i < dim; ++i)
  {
    result[i] = a[i] * scalar;
  }
  return result;
}

// allow multiplication with a scalar before a point

template<size_t dim>
Point<dim> operator*(double scalar, const Point<dim>& a)
{
  Point<dim> result {};
  for (size_t i = 0; i < dim; ++i)
  {
    result[i] = a[i] * scalar;
  }
  return result;
}

// also provide Vector as alias
template<size_t dim>
using Vector = Point<dim>;
}
