#pragma once
#include "../eigen/Eigen/Dense"
#include "../eigen/unsupported/Eigen/Polynomials"
#include <iostream>

namespace kinDS
{
  struct CoefficientTerm {
    double coefficient; // Coefficient of the term
    int exponent;       // Exponent of the term
    CoefficientTerm(double coeff, int exp = 0) : coefficient(coeff), exponent(exp) {}

    // unary minus for CoefficientTerm
    CoefficientTerm operator-() const {
      return { -coefficient, exponent };
    }
  };

  // Allow polyonomials to be expressed symbolically
  struct Var {
    // overload ^ operator for powers
    CoefficientTerm operator^(int exp) const {
      return { 1, exp }; // coefficient 1, exponent `exp`
    }

    // allow X by itself: X == X^1
    operator CoefficientTerm() const {
      return { 1, 1 };
    }
  };

  // allow multiplication: 3 * (X^2)
  template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  CoefficientTerm operator*(T coeff, const CoefficientTerm& term) {
    return { coeff * term.coefficient, term.exponent };
  }

  // allow multiplication: (X^2) * 3
  template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  CoefficientTerm operator*(const CoefficientTerm& term, T coeff) {
    return { coeff * term.coefficient, term.exponent };
  }

  /**
   * \brief Class representing a polynomial with real coefficients.
   */
  class Polynomial {
    Eigen::VectorXd coeffs;

  private:
    static Eigen::VectorXd add_poly(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
      int n = std::max(a.size(), b.size());
      Eigen::VectorXd result = Eigen::VectorXd::Zero(n);
      result.head(a.size()) += a;
      result.head(b.size()) += b;
      return result;
    }

    static Eigen::VectorXd multiply_poly(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
      Eigen::VectorXd result = Eigen::VectorXd::Zero(a.size() + b.size() - 1);
      for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < b.size(); ++j) {
          result[i + j] += a[i] * b[j];
        }
      }
      return result;
    }

  public:
    Polynomial(const Eigen::VectorXd& c) : coeffs(c) {}

    Polynomial(std::initializer_list<CoefficientTerm> init) {
      // Determine the maximum exponent
      int max_exp = 0;
      for (auto& [coeff, exp] : init) {
        if (exp > max_exp) {
          max_exp = exp;
        }
      }

      coeffs = Eigen::VectorXd::Zero(max_exp + 1);

      for (auto& [coeff, exp] : init) {
        if (exp < 0) {
          throw std::invalid_argument("Negative exponent in polynomial initialization.");
        }
        coeffs[exp] += coeff; // Add coefficients for the same exponent
      }
    }

    Polynomial(const Var& x) : coeffs(Eigen::VectorXd::Zero(2)) {
      coeffs[1] = 1; // x^1
    }

    Polynomial(double constant) : coeffs(Eigen::VectorXd::Zero(1)) {
      coeffs[0] = constant; // Constant polynomial
    }

    Polynomial(int constant) : coeffs(Eigen::VectorXd::Zero(1)) {
      coeffs[0] = static_cast<double>(constant); // Constant polynomial
    }

    Polynomial(std::function<Polynomial(const Var&)> func) {
      Var X;
      Polynomial poly = func(X);
      coeffs = poly.coeffs;
    }

    double eval(double x) const {
      double result = 0;
      for (int i = coeffs.size() - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
      }
      return result;
    }

    Eigen::VectorXcd roots() const {
      Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
      solver.compute(coeffs);
      Eigen::VectorXcd result = solver.roots();
      return result;
    }

    Eigen::VectorXd realRoots() const {
      Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
      solver.compute(coeffs);
      const Eigen::VectorXcd& complex_roots = solver.roots();

      Eigen::VectorXd real_roots;
      for (auto& root : complex_roots) {
        if (root.imag() == 0)
        {
          real_roots << root.real();
        }
      }

      return real_roots;
    }

    Polynomial derivative(size_t order = 1) const {
      if (order == 0) {
        return *this; // The zeroth derivative is the polynomial itself
      }
      Eigen::VectorXd derived_coeffs = coeffs;

      // do computation in place
      for (size_t i = 0; i < order; ++i) {
        for (int j = 1; j < derived_coeffs.size(); ++j) {
          derived_coeffs[j - 1] = derived_coeffs[j] * j;
        }
        derived_coeffs.tail(1).setZero(); // Set the last coefficient to zero after differentiation
        derived_coeffs.resize(derived_coeffs.size() - 1);
      }
      return Polynomial(derived_coeffs);
    }

    Polynomial operator+(const Polynomial& other) const {
      return Polynomial(add_poly(coeffs, other.coeffs));
    }
    Polynomial operator-(const Polynomial& other) const {
      return Polynomial(add_poly(coeffs, -other.coeffs));
    }

    // unary minus
    Polynomial operator-() const {
      return Polynomial(-coeffs);
    }

    Polynomial operator*(const Polynomial& other) const {
      return Polynomial(multiply_poly(coeffs, other.coeffs));
    }

    // allow multiplication with a scalar
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    Polynomial operator*(T scalar) const {
      return Polynomial(coeffs * scalar);
    }

    // allow multiplication with a scalar after a term
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    friend Polynomial operator*(const Polynomial& poly, T scalar) {
      return Polynomial(poly.coeffs * scalar);
    }

    // Operators for adding symbolic terms
    /*Polynomial operator+(const CoefficientTerm& term) const {
      Polynomial result = *this;
      // Ensure the polynomial has enough coefficients
      if (result.coeffs.size() <= term.exponent) {
        size_t old_size = result.coeffs.size();
        result.coeffs.conservativeResize(term.exponent + 1);
        result.coeffs.tail(term.exponent + 1 - old_size).setZero(); // Ensure new coefficients are zero
      }

      result.coeffs[term.exponent] += term.coefficient;
      return result;
    }*/
    

    void print() const {
      for (int i = coeffs.size() - 1; i >= 0; --i)
        if (coeffs[i] != 0)
        {
          std::cout << (i != coeffs.size() - 1 && coeffs[i] > 0 ? "+" : "")
            << coeffs[i];

          if (i > 1)
            std::cout << "x^" << i << " ";
          else if (i == 1)
            std::cout << "x ";
        }
      std::cout << "\n";
    }
  };

  // combine coeffient terms to a polynomial
  Polynomial operator+(const CoefficientTerm& term1, const CoefficientTerm& term2);

  Polynomial operator-(const CoefficientTerm& term1, const CoefficientTerm& term2);

#define POLYNOMIAL(expr) kinDS::Polynomial([&](kinDS::Var x) { return (kinDS::Polynomial{(expr)}); })
} // namespace kinDS