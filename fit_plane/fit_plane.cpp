#include "fit_plane.h"

#include <stdexcept>
#include <cmath>
#include <fmt/printf.h>

fit_plane_calculus::fit_plane_calculus():
  sum_xi(0.0), sum_yi(0.0), sum_zi(0.0),
  sum_xi_yi(0.0), sum_xi_zi(0.0), sum_xi_square(0.0),
  sum_yi_square(0.0), g(0.0), norm(0.0), k0(0.0),
  k1(0.0), k2(0.0) {}

void fit_plane_calculus::calc_value(const std::vector<std::vector<double>>& points) {
  sum_xi = 0;
  sum_yi = 0;
  sum_zi = 0;
  sum_xi_yi = 0;
  sum_xi_zi = 0;
  sum_yi_zi = 0;
  sum_xi_square = 0;
  sum_yi_square = 0;
  g = 0;
  m_points = points;
  const double n = m_points.size();
  for (const auto& p: m_points) {
    if (p.size() == 3) {
      sum_xi += p[0];
      sum_yi += p[1];
      sum_zi += p[2];
      sum_xi_yi += p[0] * p[1];
      sum_xi_zi += p[0] * p[2];
      sum_yi_zi += p[1] * p[2];
      sum_xi_square += p[0] * p[0];
      sum_yi_square += p[1] * p[1];
    } else {
      throw std::runtime_error("invalid point size.\n");
    }
  }
  g = sum_xi_square * (n * sum_yi_square - sum_yi * sum_yi) - sum_xi_yi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_xi * (sum_xi_yi * sum_yi - sum_xi * sum_yi_square);
  if (g == 0) {
    throw std::runtime_error("The plane may go through the Z axis! Please consider applying a restraint to circumvent this.\n");
  }
  k0 = (sum_xi_zi * (sum_yi * sum_yi - n * sum_yi_square) + sum_yi_zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_zi * (sum_xi * sum_yi_square - sum_yi * sum_xi_yi)) / g;
  k1 = (sum_xi_zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_yi_zi * (sum_xi * sum_xi - n * sum_xi_square) + sum_zi * (sum_xi_square * sum_yi - sum_xi * sum_xi_yi)) / g;
  k2 = 1.0;
  k2_alt = (sum_xi_zi * (sum_yi_square * sum_xi - sum_yi * sum_xi_yi) +
            sum_yi_zi * (sum_xi_square * sum_yi - sum_xi * sum_xi_yi) +
            sum_zi * (sum_xi_yi * sum_xi_yi - sum_xi_square * sum_yi_square)) / g;
  norm = 1.0 / std::sqrt(k0 * k0 + k1 * k1 + 1);
//   fmt::print("k0 = {:12.7f} ; k1 = {:12.7f} ; norm = {:12.7f}\n",
//              k0, k1, norm);
}

std::tuple<double, double, double> fit_plane_calculus::normal_vector() const {
  return std::make_tuple(k0 * norm, k1 * norm, k2 * norm);
}

std::tuple<double, double, double> fit_plane_calculus::factors() const {
  return std::make_tuple(-k0, -k1, -k2_alt);
}

void fit_plane_calculus::calc_gradients(std::vector<std::vector<double>>& gradients) {
  // gradients
  const double n = m_points.size();
  const double f1 = sum_xi_zi * (sum_yi * sum_yi - n * sum_yi_square) + sum_yi_zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_zi * (sum_xi * sum_yi_square - sum_yi * sum_xi_yi);
  const double f2 = sum_xi_zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_yi_zi * (sum_xi * sum_xi - n * sum_xi_square) + sum_zi * (sum_xi_square * sum_yi - sum_xi * sum_xi_yi);
  gradients.resize(n, std::vector<double>(3, 0));
  for (size_t i = 0; i < m_points.size(); ++i) {
    const double& xi = m_points[i][0];
    const double& yi = m_points[i][1];
    const double& zi = m_points[i][2];
    const double dg_dxi = 2 * xi * (n * sum_yi_square - sum_yi * sum_yi) - yi * (n * sum_xi_yi - sum_xi * sum_yi) - sum_xi_yi * (n * yi - sum_yi) + sum_xi_yi * sum_yi - sum_xi * sum_yi_square + yi * sum_xi * sum_yi - sum_xi * sum_yi_square;
    const double dg_dyi = 2 * n * yi * sum_xi_square - 2 * sum_xi_square * sum_yi - n * xi * sum_xi_yi + xi * sum_xi * sum_yi - n * xi * sum_xi_yi + sum_xi * sum_xi_yi + xi * sum_xi * sum_yi + sum_xi * sum_xi_yi - 2 * yi * sum_xi * sum_xi;
    const double df1_dxi = zi * (sum_yi * sum_yi - n * sum_yi_square) + sum_yi_zi * (n * yi - sum_yi) + sum_zi * (sum_yi_square - yi * sum_yi);
    const double df1_dyi = sum_xi_zi * (2 * sum_yi - 2 * n * yi) + zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_yi_zi * (n * xi - sum_xi) + sum_zi * (2 * yi * sum_xi - (sum_xi_yi + xi * sum_yi));
    const double df1_dzi = xi * (sum_yi * sum_yi - n * sum_yi_square) + yi * (n * sum_xi_yi - sum_xi * sum_yi) + (sum_xi * sum_yi_square - sum_yi * sum_xi_yi);
    const double dk0_dxi = (df1_dxi * g - dg_dxi * f1) / (g * g);
    const double dk0_dyi = (df1_dyi * g - dg_dyi * f1) / (g * g);
    const double dk0_dzi = df1_dzi / g;
    const double df2_dxi = zi * (n * sum_xi_yi - sum_xi * sum_yi) + sum_xi_zi * (n * yi - sum_yi) + sum_yi_zi * (2 * sum_xi - 2 * n *xi) + sum_zi * (2 * xi * sum_yi - (sum_xi_yi + yi * sum_xi));
    const double df2_dyi = sum_xi_zi * (n * xi - sum_xi) + zi * (sum_xi * sum_xi - n * sum_xi_square) + sum_zi * (sum_xi_square - xi * sum_xi);
    const double df2_dzi = xi * (n * sum_xi_yi - sum_xi * sum_yi) + yi * (sum_xi * sum_xi - n * sum_xi_square) + (sum_xi_square * sum_yi - sum_xi * sum_xi_yi);
    const double dk1_dxi = (df2_dxi * g - dg_dxi * f2) / (g * g);
    const double dk1_dyi = (df2_dyi * g - dg_dyi * f2) / (g * g);
    const double dk1_dzi = df2_dzi / g;
    const double dnorm_dxi = 0.5 * norm * (2 * k0 * dk0_dxi + 2 * k1 * dk1_dxi) * (-1.0 / (k0 * k0 + k1 * k1 + 1));
    const double dnorm_dyi = 0.5 * norm * (2 * k0 * dk0_dyi + 2 * k1 * dk1_dyi) * (-1.0 / (k0 * k0 + k1 * k1 + 1));
    const double dnorm_dzi = 0.5 * norm * (2 * k0 * dk0_dzi + 2 * k1 * dk1_dzi) * (-1.0 / (k0 * k0 + k1 * k1 + 1));
    gradients[i][0] = dk0_dxi * norm + k0 * dnorm_dxi;
    gradients[i][1] = dk0_dyi * norm + k0 * dnorm_dyi;
    gradients[i][2] = dk0_dzi * norm + k0 * dnorm_dzi;
    gradients[i][0] = dk1_dxi * norm + k1 * dnorm_dxi;
    gradients[i][1] = dk1_dyi * norm + k1 * dnorm_dyi;
    gradients[i][2] = dk1_dzi * norm + k1 * dnorm_dzi;
    gradients[i][0] = 1.0 * dnorm_dxi;
    gradients[i][1] = 1.0 * dnorm_dyi;
    gradients[i][2] = 1.0 * dnorm_dzi;
  }
}
