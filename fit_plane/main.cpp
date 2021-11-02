#include "fit_plane.h"

#include <vector>
#include <fmt/printf.h>
#include <Eigen/Dense>
#include <iostream>

int main() {
  std::vector<std::vector<double>> points{
    {-1.0, 1.5, 2.0},
    {-6.2, -6.5, -4.3},
    {2.3, 1.9, 0.5},
    {-1.0, 1.0, -1.4},
    {7.6, 1.8, -6.9}
  };
  fit_plane_calculus f1;
  f1.calc_value(points);
  const auto eq = f1.factors();
  fmt::print("Solution using calculus:\n");
  fmt::print("k0 = {:12.7f}\nk1 = {:12.7f}\nk2 = {:12.7f}\n",
             std::get<0>(eq), std::get<1>(eq), std::get<2>(eq));
  // use linear algebra
  Eigen::MatrixXd A(points.size(), 3);
  Eigen::VectorXd b(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    A(i, 0) = points[i][0];
    A(i, 1) = points[i][1];
    A(i, 2) = 1.0;
    b(i) = points[i][2];
  }
//   std::cout << A << std::endl;
//   std::cout << b << std::endl;
  // normal equations
  const auto L = (A.transpose() * A).ldlt();
  const auto solution_linalg = L.solve(A.transpose() * b);
  std::cout << "Solution using linear algebra:\n";
  fmt::print("k0 = {:12.7f}\nk1 = {:12.7f}\nk2 = {:12.7f}\n",
             solution_linalg(0), solution_linalg(1), solution_linalg(2));
  return 0;
}
