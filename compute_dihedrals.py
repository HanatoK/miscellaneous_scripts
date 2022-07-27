#!/usr/bin/env python3
import numpy as np


def compute_dihedral_angles(X):
    # X is a (4,3) array
    vec_ij = X[1] - X[0]
    vec_jk = X[2] - X[1]
    vec_kl = X[3] - X[2]
    normal_kji = np.cross(vec_ij, vec_jk)
    normal_jkl = np.cross(vec_jk, vec_kl)
    cos_alpha = np.sum(normal_kji * normal_jkl)
    sin_alpha = np.sum(normal_kji * vec_kl) * \
        np.sqrt(np.sum(np.square(vec_jk)))
    return np.arctan2(sin_alpha, cos_alpha)


def test():
    X = {1: [-0.197758333931069, -0.250678894164538, -0.02368355028707],
         5: [-0.0957022165356829, -0.151713094094331, 0.0136232367728677],
         7: [-0.137798040706579, -0.0269396992302593, 0.0210804173350176],
         9: [-0.0499756353480867, 0.080045033821632, 0.0541627679311519],
         11: [-0.00173298734984859, 0.0602820233654584, 0.197418755186087],
         15: [0.049248461867456, 0.128750167408009, -0.050832145642861],
         17: [0.157422330773366, 0.0606799517818354, -0.0660578969720305],
         19: [0.276296421230445, 0.0995745111121945, -0.145711584323163]}
    atom_coordinates = np.array([X[5], X[7], X[9], X[15]])
    result = compute_dihedral_angles(atom_coordinates)
    print(np.degrees(result))
    atom_coordinates = np.array([X[7], X[9], X[15], X[17]])
    result = compute_dihedral_angles(atom_coordinates)
    print(np.degrees(result))


if __name__ == '__main__':
    test()
