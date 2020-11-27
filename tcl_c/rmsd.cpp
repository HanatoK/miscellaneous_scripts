#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "rmsd.h"

// equation (5)

std::array<std::array<double, 4>, 4> getMatrixF(const AtomSet &a,
                                                const AtomSet &b) {
  std::array<std::array<double, 4>, 4> matrix_F;
  double R11 = 0;
  double R22 = 0;
  double R33 = 0;
  double R12 = 0;
  double R13 = 0;
  double R23 = 0;
  double R21 = 0;
  double R31 = 0;
  double R32 = 0;
  // equation (5)
  for (size_t k = 0; k < a.size(); ++k) {
    R11 += a(1, k) * b(1, k);
    R22 += a(2, k) * b(2, k);
    R33 += a(3, k) * b(3, k);
    R12 += a(1, k) * b(2, k);
    R13 += a(1, k) * b(3, k);
    R23 += a(2, k) * b(3, k);
    R21 += a(2, k) * b(1, k);
    R31 += a(3, k) * b(1, k);
    R32 += a(3, k) * b(2, k);
  }
  // equation (10)
  matrix_F[0][0] = R11 + R22 + R33;
  matrix_F[1][0] = R23 - R32;
  matrix_F[0][1] = matrix_F[1][0];
  matrix_F[2][0] = R31 - R13;
  matrix_F[0][2] = matrix_F[2][0];
  matrix_F[3][0] = R12 - R21;
  matrix_F[0][3] = matrix_F[3][0];
  matrix_F[1][1] = R11 - R22 - R33;
  matrix_F[1][2] = R12 + R21;
  matrix_F[2][1] = matrix_F[1][2];
  matrix_F[1][3] = R13 + R31;
  matrix_F[3][1] = matrix_F[1][3];
  matrix_F[2][2] = -R11 + R22 - R33;
  matrix_F[2][3] = R23 + R32;
  matrix_F[3][2] = matrix_F[2][3];
  matrix_F[3][3] = -R11 - R22 + R33;
  return matrix_F;
}

// only this function depends the Eigen library
void computeEigen(const std::array<std::array<double, 4>, 4> &matrix_F,
                  std::array<double, 4> &eigenvalues,
                  std::array<std::array<double, 4>, 4> &eigenvectors) {
  Eigen::Matrix4d mat;
  mat << matrix_F[0][0], matrix_F[1][0], matrix_F[2][0], matrix_F[3][0],
      matrix_F[0][1], matrix_F[1][1], matrix_F[2][1], matrix_F[3][1],
      matrix_F[0][2], matrix_F[1][2], matrix_F[2][2], matrix_F[3][2],
      matrix_F[0][3], matrix_F[1][3], matrix_F[2][3], matrix_F[3][3];
  // compute eigenvalues and eigenvectors
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver(mat);
  if (eigensolver.info() != Eigen::Success)
    abort();
  eigenvalues[0] = eigensolver.eigenvalues()[0];
  eigenvalues[1] = eigensolver.eigenvalues()[1];
  eigenvalues[2] = eigensolver.eigenvalues()[2];
  eigenvalues[3] = eigensolver.eigenvalues()[3];
  eigenvectors[0][0] = eigensolver.eigenvectors().col(0)[0];
  eigenvectors[0][1] = eigensolver.eigenvectors().col(1)[0];
  eigenvectors[0][2] = eigensolver.eigenvectors().col(2)[0];
  eigenvectors[0][3] = eigensolver.eigenvectors().col(3)[0];
  eigenvectors[1][0] = eigensolver.eigenvectors().col(0)[1];
  eigenvectors[1][1] = eigensolver.eigenvectors().col(1)[1];
  eigenvectors[1][2] = eigensolver.eigenvectors().col(2)[1];
  eigenvectors[1][3] = eigensolver.eigenvectors().col(3)[1];
  eigenvectors[2][0] = eigensolver.eigenvectors().col(0)[2];
  eigenvectors[2][1] = eigensolver.eigenvectors().col(1)[2];
  eigenvectors[2][2] = eigensolver.eigenvectors().col(2)[2];
  eigenvectors[2][3] = eigensolver.eigenvectors().col(3)[2];
  eigenvectors[3][0] = eigensolver.eigenvectors().col(0)[3];
  eigenvectors[3][1] = eigensolver.eigenvectors().col(1)[3];
  eigenvectors[3][2] = eigensolver.eigenvectors().col(2)[3];
  eigenvectors[3][3] = eigensolver.eigenvectors().col(3)[3];
  // normalization
  double norm[4];
  norm[0] = std::sqrt(eigenvectors[0][0] * eigenvectors[0][0] +
                      eigenvectors[0][1] * eigenvectors[0][1] +
                      eigenvectors[0][2] * eigenvectors[0][2] +
                      eigenvectors[0][3] * eigenvectors[0][3]);
  norm[1] = std::sqrt(eigenvectors[1][0] * eigenvectors[1][0] +
                      eigenvectors[1][1] * eigenvectors[1][1] +
                      eigenvectors[1][2] * eigenvectors[1][2] +
                      eigenvectors[1][3] * eigenvectors[1][3]);
  norm[2] = std::sqrt(eigenvectors[2][0] * eigenvectors[2][0] +
                      eigenvectors[2][1] * eigenvectors[2][1] +
                      eigenvectors[2][2] * eigenvectors[2][2] +
                      eigenvectors[2][3] * eigenvectors[2][3]);
  norm[3] = std::sqrt(eigenvectors[3][0] * eigenvectors[3][0] +
                      eigenvectors[3][1] * eigenvectors[3][1] +
                      eigenvectors[3][2] * eigenvectors[3][2] +
                      eigenvectors[3][3] * eigenvectors[3][3]);
  eigenvectors[0][0] /= norm[0];
  eigenvectors[0][1] /= norm[0];
  eigenvectors[0][2] /= norm[0];
  eigenvectors[0][3] /= norm[0];
  eigenvectors[1][0] /= norm[1];
  eigenvectors[1][1] /= norm[1];
  eigenvectors[1][2] /= norm[1];
  eigenvectors[1][3] /= norm[1];
  eigenvectors[2][0] /= norm[2];
  eigenvectors[2][1] /= norm[2];
  eigenvectors[2][2] /= norm[2];
  eigenvectors[2][3] /= norm[2];
  eigenvectors[3][0] /= norm[3];
  eigenvectors[3][1] /= norm[3];
  eigenvectors[3][2] /= norm[3];
  eigenvectors[3][3] /= norm[3];
}

std::array<std::array<double, 3>, 3> rotationMatrixAtoB(AtomSet a, AtomSet b) {
  std::array<double, 4> eigenvalues;
  std::array<std::array<double, 4>, 4> eigenvectors;
  a.shiftToCenter();
  b.shiftToCenter();
  computeEigen(getMatrixF(a, b), eigenvalues, eigenvectors);
#ifdef DEBUG
  std::cout << "Eigenvalues: " << '\n';
  for (size_t i = 0; i < 4; ++i) {
    std::cout << eigenvalues[i] << ' ';
  }
  std::cout << '\n';
  std::cout << "Eigenvectors (in collumn vectors): " << '\n';
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      std::cout << eigenvectors[i][j] << ' ';
    }
    std::cout << '\n';
  }
#endif
  // find the maximum eigenvalue
  const size_t max_index =
      std::distance(eigenvalues.begin(),
                    std::max_element(eigenvalues.begin(), eigenvalues.end()));
  std::array<double, 4> q;
  q[0] = eigenvectors[0][max_index];
  q[1] = eigenvectors[1][max_index];
  q[2] = eigenvectors[2][max_index];
  q[3] = eigenvectors[3][max_index];
  std::array<std::array<double, 3>, 3> mat;
  mat[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  mat[0][1] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  mat[0][2] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
  mat[1][0] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  mat[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  mat[1][2] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  mat[2][0] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  mat[2][1] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  mat[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
  return mat;
}

double rmsd(const AtomSet &a, const AtomSet &b) {
  double rmsd = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    rmsd +=
        (a.m_atoms[i].x - b.m_atoms[i].x) * (a.m_atoms[i].x - b.m_atoms[i].x);
    rmsd +=
        (a.m_atoms[i].y - b.m_atoms[i].y) * (a.m_atoms[i].y - b.m_atoms[i].y);
    rmsd +=
        (a.m_atoms[i].z - b.m_atoms[i].z) * (a.m_atoms[i].z - b.m_atoms[i].z);
  }
  rmsd /= a.size();
  rmsd = std::sqrt(rmsd);
  return rmsd;
}

// minimal rmsd using equation (9)
double minimalRMSD(AtomSet a, AtomSet b) {
  std::array<double, 4> eigenvalues;
  std::array<std::array<double, 4>, 4> eigenvectors;
  a.shiftToCenter();
  b.shiftToCenter();
  computeEigen(getMatrixF(a, b), eigenvalues, eigenvectors);
  const size_t max_index =
      std::distance(eigenvalues.begin(),
                    std::max_element(eigenvalues.begin(), eigenvalues.end()));
  double minimal_rmsd = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    minimal_rmsd +=
        a.m_atoms[i].x * a.m_atoms[i].x + b.m_atoms[i].x * b.m_atoms[i].x;
    minimal_rmsd +=
        a.m_atoms[i].y * a.m_atoms[i].y + b.m_atoms[i].y * b.m_atoms[i].y;
    minimal_rmsd +=
        a.m_atoms[i].z * a.m_atoms[i].z + b.m_atoms[i].z * b.m_atoms[i].z;
  }
  minimal_rmsd -= 2.0 * eigenvalues[max_index];
  minimal_rmsd = std::sqrt(minimal_rmsd / a.size());
  return minimal_rmsd;
}

double minimalRMSD(const std::vector<std::vector<double>> &a,
                   const std::vector<std::vector<double>> &b) {
  AtomSet atom_set_a(a);
  AtomSet atom_set_b(b);
  return minimalRMSD(atom_set_a, atom_set_b);
}

// C interface
double minimalRMSD(double **a, double **b, const int len) {
  std::vector<std::vector<double>> vec_a(len, std::vector<double>(3));
  std::vector<std::vector<double>> vec_b(len, std::vector<double>(3));
  for (int i = 0; i < len; ++i) {
    vec_a[i][0] = a[i][0];
    vec_a[i][1] = a[i][1];
    vec_a[i][2] = a[i][2];
    vec_b[i][0] = b[i][0];
    vec_b[i][1] = b[i][1];
    vec_b[i][2] = b[i][2];
  }
  return minimalRMSD(vec_a, vec_b);
}

// C interface
// rotate A to optimally align to B
void optimalRotateAtoB(double **a, double **b, double **rotated_a,
                       const int len) {
  std::vector<std::vector<double>> vec_a(len, std::vector<double>(3));
  std::vector<std::vector<double>> vec_b(len, std::vector<double>(3));
  for (int i = 0; i < len; ++i) {
    vec_a[i][0] = a[i][0];
    vec_a[i][1] = a[i][1];
    vec_a[i][2] = a[i][2];
    vec_b[i][0] = b[i][0];
    vec_b[i][1] = b[i][1];
    vec_b[i][2] = b[i][2];
  }
  AtomSet atom_set_a(vec_a);
  AtomSet atom_set_b(vec_b);
  AtomSet atom_set_a_r(vec_a);
  auto mat = rotationMatrixAtoB(atom_set_a, atom_set_b);
  atom_set_a_r.rotate(mat);
  for (int i = 0; i < len; ++i) {
    rotated_a[i][0] = atom_set_a_r.m_atoms[i].x;
    rotated_a[i][1] = atom_set_a_r.m_atoms[i].y;
    rotated_a[i][2] = atom_set_a_r.m_atoms[i].z;
  }
}

// C interface
void rotationMatrixAtoB(double **a, double **b, double (*res)[3][3],
                        const int len) {
  std::vector<std::vector<double>> vec_a(len, std::vector<double>(3));
  std::vector<std::vector<double>> vec_b(len, std::vector<double>(3));
  for (int i = 0; i < len; ++i) {
    vec_a[i][0] = a[i][0];
    vec_a[i][1] = a[i][1];
    vec_a[i][2] = a[i][2];
    vec_b[i][0] = b[i][0];
    vec_b[i][1] = b[i][1];
    vec_b[i][2] = b[i][2];
  }
  AtomSet atom_set_a(vec_a);
  AtomSet atom_set_b(vec_b);
  auto mat = rotationMatrixAtoB(atom_set_a, atom_set_b);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      (*res)[i][j] = mat[i][j];
    }
  }
}

// C interface
void bringToCenter(double **atom_positions, const int num_atoms) {
  // compute center-of-geometry
  std::array<double, 3> center_of_geometry{0, 0, 0};
  for (int i = 0; i < num_atoms; ++i) {
    for (int j = 0; j < 3; ++j) {
      center_of_geometry[j] += atom_positions[i][j];
    }
  }
  for (int j = 0; j < 3; ++j) {
    center_of_geometry[j] /= static_cast<double>(num_atoms);
  }
  // move atoms
  for (int i = 0; i < num_atoms; ++i) {
    for (int j = 0; j < 3; ++j) {
      atom_positions[i][j] -= center_of_geometry[j];
    }
  }
} 
