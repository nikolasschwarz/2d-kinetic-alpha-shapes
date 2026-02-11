#include "kinDS/Polynomial.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using namespace kinDS;
using Catch::Approx;

TEST_CASE("Polynomial creation from Eigen vectors", "[Polynomial]")
{
  Eigen::VectorXd a(3);
  a << 1, 2, 3; // 1 + 2x + 3x^2
  Eigen::VectorXd b(2);
  b << 4, 5; // 4 + 5x

  Polynomial p1(a);
  Polynomial p2(b);

  // Test evaluation at x=0
  REQUIRE(p1(0.0) == Approx(1.0));
  REQUIRE(p2(0.0) == Approx(4.0));

  // Test evaluation at x=1
  REQUIRE(p1(1.0) == Approx(6.0)); // 1 + 2 + 3 = 6
  REQUIRE(p2(1.0) == Approx(9.0)); // 4 + 5 = 9
}

TEST_CASE("Polynomial addition", "[Polynomial]")
{
  Eigen::VectorXd a(3);
  a << 1, 2, 3; // 1 + 2x + 3x^2
  Eigen::VectorXd b(2);
  b << 4, 5; // 4 + 5x

  Polynomial p1(a);
  Polynomial p2(b);
  Polynomial sum = p1 + p2;

  // Sum should be: (1+4) + (2+5)x + 3x^2 = 5 + 7x + 3x^2
  REQUIRE(sum(0.0) == Approx(5.0));
  REQUIRE(sum(1.0) == Approx(15.0)); // 5 + 7 + 3 = 15
  REQUIRE(sum(2.0) == Approx(31.0)); // 5 + 14 + 12 = 31
}

TEST_CASE("Polynomial multiplication", "[Polynomial]")
{
  Eigen::VectorXd a(3);
  a << 1, 2, 3; // 1 + 2x + 3x^2
  Eigen::VectorXd b(2);
  b << 4, 5; // 4 + 5x

  Polynomial p1(a);
  Polynomial p2(b);
  Polynomial prod = p1 * p2;

  // Product: (1 + 2x + 3x^2) * (4 + 5x)
  // = 4 + 5x + 8x + 10x^2 + 12x^2 + 15x^3
  // = 4 + 13x + 22x^2 + 15x^3
  REQUIRE(prod(0.0) == Approx(4.0));
  REQUIRE(prod(1.0) == Approx(54.0)); // 4 + 13 + 22 + 15 = 54
  REQUIRE(prod(2.0) == Approx(238.0));
}

TEST_CASE("Polynomial evaluation", "[Polynomial]")
{
  Eigen::VectorXd a(3);
  a << 1, 2, 3; // 1 + 2x + 3x^2
  Eigen::VectorXd b(2);
  b << 4, 5; // 4 + 5x

  Polynomial p1(a);
  Polynomial p2(b);
  Polynomial sum = p1 + p2;
  Polynomial prod = p1 * p2;

  double x = 2.0;
  double result_sum = sum(x);
  double result_prod = prod(x);

  // Sum at x=2: 5 + 7*2 + 3*4 = 5 + 14 + 12 = 31
  REQUIRE(result_sum == Approx(31.0));

  // Product at x=2: 4 + 13*2 + 22*4 + 15*8 = 4 + 26 + 88 + 120 = 238
  REQUIRE(result_prod == Approx(238.0));
}

TEST_CASE("Polynomial creation with POLYNOMIAL macro", "[Polynomial]")
{
  // Test POLYNOMIAL macro: POLYNOMIAL(x^2) should create x^2
  auto test = POLYNOMIAL(x ^ 2); // X^2
  REQUIRE(test(0.0) == Approx(0.0));
  REQUIRE(test(1.0) == Approx(1.0));
  REQUIRE(test(2.0) == Approx(4.0));
  REQUIRE(test(3.0) == Approx(9.0));
}

TEST_CASE("Polynomial creation with lambda", "[Polynomial]")
{
  // Test creating polynomial from lambda expression
  Polynomial p3 = Polynomial([&](Var x) { return (Polynomial { (10 * (x ^ 4) + 4 * (x ^ 2) + 2) }); });

  // p3 = 10x^4 + 4x^2 + 2
  REQUIRE(p3(0.0) == Approx(2.0));
  REQUIRE(p3(1.0) == Approx(16.0)); // 10 + 4 + 2 = 16
  REQUIRE(p3(2.0) == Approx(178.0)); // 10*16 + 4*4 + 2 = 160 + 16 + 2 = 178
}

TEST_CASE("Polynomial with negative coefficients", "[Polynomial]")
{
  Eigen::VectorXd c(3);
  c << 4, 5, -2; // 4 + 5x - 2x^2

  Polynomial p3(c);

  REQUIRE(p3(0.0) == Approx(4.0));
  REQUIRE(p3(1.0) == Approx(7.0)); // 4 + 5 - 2 = 7
  REQUIRE(p3(2.0) == Approx(6.0));
}
