#ifndef RMSD_H
#define RMSD_H
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

struct AtomPosition {
  double x;
  double y;
  double z;
};

class AtomSet {
public:
  AtomSet() {}
  AtomSet(const std::vector<AtomPosition> &atoms) : m_atoms(atoms) {}
  AtomSet(const std::vector<std::vector<double>> &pos_set) {
    m_atoms.resize(pos_set.size());
    for (size_t i = 0; i < pos_set.size(); ++i) {
      m_atoms[i].x = pos_set[i][0];
      m_atoms[i].y = pos_set[i][1];
      m_atoms[i].z = pos_set[i][2];
    }
  }
  void addAtom(const AtomPosition &a) { m_atoms.push_back(a); }
  void addDummyAtom(double x, double y, double z) {
    AtomPosition a;
    a.x = x;
    a.y = y;
    a.z = z;
    m_atoms.push_back(a);
  }
  size_t size() const { return m_atoms.size(); }
  double operator()(int i, size_t k) const {
    switch (i) {
    case 1:
      return m_atoms[k].x;
      break;
    case 2:
      return m_atoms[k].y;
      break;
    case 3:
      return m_atoms[k].z;
      break;
    default:
      return 0.0;
    }
  }
  std::array<double, 3> shiftToCenter() {
    std::array<double, 3> com;
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;
    for (const auto &a : m_atoms) {
      com[0] += a.x;
      com[1] += a.y;
      com[2] += a.z;
    }
    com[0] /= size();
    com[1] /= size();
    com[2] /= size();
    for (auto &a : m_atoms) {
      a.x -= com[0];
      a.y -= com[1];
      a.z -= com[2];
    }
    return com;
  }
  void rotate(const std::array<std::array<double, 3>, 3> &rotate_matrix) {
    for (size_t i = 0; i < size(); ++i) {
      const AtomPosition tmp_atom = m_atoms[i];
      m_atoms[i].x = rotate_matrix[0][0] * tmp_atom.x +
                     rotate_matrix[0][1] * tmp_atom.y +
                     rotate_matrix[0][2] * tmp_atom.z;
      m_atoms[i].y = rotate_matrix[1][0] * tmp_atom.x +
                     rotate_matrix[1][1] * tmp_atom.y +
                     rotate_matrix[1][2] * tmp_atom.z;
      m_atoms[i].z = rotate_matrix[2][0] * tmp_atom.x +
                     rotate_matrix[2][1] * tmp_atom.y +
                     rotate_matrix[2][2] * tmp_atom.z;
    }
  }
  void dumpAtoms(std::ostream &os) const {
    for (const auto &atom : m_atoms) {
      os << "{ " << atom.x << " " << atom.y << " " << atom.z << " }"
         << " ";
    }
    os << std::endl;
  }
  std::vector<AtomPosition> m_atoms;
};

// debug printing
template <typename T, size_t N>
void printMatrix2D(const std::array<std::array<T, N>, N> &mat) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      std::cout << mat[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

// equation (5)

std::array<std::array<double, 4>, 4> getMatrixF(const AtomSet &a,
                                                const AtomSet &b);

// only this function depends the Eigen library
void computeEigen(const std::array<std::array<double, 4>, 4> &matrix_F,
                  std::array<double, 4> &eigenvalues,
                  std::array<std::array<double, 4>, 4> &eigenvectors);

std::array<std::array<double, 3>, 3> rotationMatrixAtoB(AtomSet a, AtomSet b);

double rmsd(const AtomSet &a, const AtomSet &b);

// minimal rmsd using equation (9)
double minimalRMSD(AtomSet a, AtomSet b);

double minimalRMSD(const std::vector<std::vector<double>> &a,
                   const std::vector<std::vector<double>> &b);

// C interface
double minimalRMSD(double **a, double **b, const int len);

// C interface
// rotate A to optimally align to B
void optimalRotateAtoB(double **a, double **b, double **rotated_a,
                       const int len);

// C interface
void rotationMatrixAtoB(double **a, double **b, double (*res)[3][3],
                        const int len);

// C interface
void bringToCenter(double **atom_positions, const int num_atoms);

#endif // RMSD_H
