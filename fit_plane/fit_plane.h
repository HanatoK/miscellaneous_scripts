#ifndef FIT_PLANE_H
#define FIT_PLANE_H

#include <vector>
#include <tuple>

class fit_plane_calculus {
public:
  fit_plane_calculus();
  void calc_value(const std::vector<std::vector<double>>& points);
  std::tuple<double, double, double> normal_vector() const;
  std::tuple<double, double, double> factors() const;
  void calc_gradients(std::vector<std::vector<double>>& gradients);
private:
  double sum_xi;
  double sum_yi;
  double sum_zi;
  double sum_xi_yi;
  double sum_xi_zi;
  double sum_yi_zi;
  double sum_xi_square;
  double sum_yi_square;
  double g;
  double norm;
  double k0;
  double k1;
  double k2;
  double k2_alt;
  std::vector<std::vector<double>> m_points;
};

#endif // FIT_PLANE_H
