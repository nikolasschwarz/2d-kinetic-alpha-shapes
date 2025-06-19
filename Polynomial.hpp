#pragma once
#include "eigen/Eigen/Dense"

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

  double eval(double x) const {
    double result = 0;
    for (int i = coeffs.size() - 1; i >= 0; --i) {
      result = result * x + coeffs[i];
    }
    return result;
  }

  Polynomial operator+(const Polynomial& other) const {
    return Polynomial(add_poly(coeffs, other.coeffs));
  }

  Polynomial operator*(const Polynomial& other) const {
    return Polynomial(multiply_poly(coeffs, other.coeffs));
  }

  void print() const {
    for (int i = coeffs.size() - 1; i >= 0; --i)
      if (coeffs[i] != 0)
        std::cout << (i != coeffs.size() - 1 && coeffs[i] > 0 ? "+" : "")
        << coeffs[i] << "*x^" << i << " ";
    std::cout << "\n";
  }
};
