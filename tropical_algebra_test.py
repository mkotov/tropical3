"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import sys
import unittest
import tropical_algebra
import matrix_tools
import upper_matrix_tools
import lower_matrix_tools
import random


class TestTropicalAlgebra(unittest.TestCase):
    def test_min_plus_algebra(self):
        R = tropical_algebra.R_min_plus()

        self.assertEqual(5, R.sum(10, 5))
        self.assertEqual(15, R.mul(10, 5))
        self.assertEqual(0, R.one())
        self.assertEqual(tropical_algebra.INFTY, R.zero())
        self.assertEqual(5, R.mul(R.one(), 5))
        self.assertEqual(5, R.sum(R.zero(), 5))

    def test_max_plus_algebra(self):
        R = tropical_algebra.R_max_plus()

        self.assertEqual(10, R.sum(10, 5))
        self.assertEqual(15, R.mul(10, 5))
        self.assertEqual(0, R.one())
        self.assertEqual(tropical_algebra.MINFTY, R.zero())
        self.assertEqual(5, R.mul(R.one(), 5))
        self.assertEqual(5, R.sum(R.zero(), 5))

    def test_matrix_min_plus_algebra(self):
        M = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 2)

        S = [[0, 2], [2, -1]]
        P = [[1, 4], [1, 7]]
        A = [[1, 2], [5, -1]]
        B = [[0, 3], [2, 8]]
        Z = [[tropical_algebra.INFTY, tropical_algebra.INFTY],
             [tropical_algebra.INFTY, tropical_algebra.INFTY]]
        U = [[0, tropical_algebra.INFTY], [tropical_algebra.INFTY, 0]]
        C = [[3, 4], [7, 1]]
        self.assertEqual(S, M.sum(A, B))
        self.assertEqual(P, M.mul(A, B))
        self.assertEqual(Z, M.zero())
        self.assertEqual(U, M.one())
        self.assertEqual(C, M.mul_by_coef(2, A))

    def test_matrix_min_plus_algebra_pwr_1(self):
        M = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 3)

        A = [[-42, 13, -96], [-28, 16, 65], [-85, 31, -75]]

        P = [[-256, -168, -277], [-209, -93, -199], [-266, -150, -256]]

        self.assertEqual(P, M.pwr(A, 3))

    def test_matrix_min_plus_algebra_poly(self):
        M = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 2)

        p = [2, 1]
        A = [[1, 2], [5, -1]]
        P = [[2, 3], [6, 0]]

        self.assertEqual(P, M.calc_poly(p, A))

    def test_matrix_min_plus_algebra_pwr_2(self):
        M = tropical_algebra.MatrixSemiring(tropical_algebra.R_max_plus(), 3)

        A = [[0, tropical_algebra.MINFTY, 2], [2, 0, 4], [1, 2, 3]]

        P = [[6, 7, 8], [8, 9, 10], [7, 8, 9]]

        self.assertEqual(P, M.pwr(A, 3))

    def test_mul_basis_upper_t_circulant_matrix_and_matrix(self):
        for i in range(10):
            R = tropical_algebra.MatrixSemiring(
                tropical_algebra.R_min_plus(), 10)
            s = random.randint(-100, 100)
            k = random.randint(0, 9)
            Y = matrix_tools.generate_random_matrix(R, -100, 100)
            B = upper_matrix_tools.generate_basis_upper_t_circulant_matrix(R, s, k)
            self.assertEqual(R.mul(
                B, Y), upper_matrix_tools.mul_basis_upper_t_circulant_matrix_and_matrix(R, s, k, Y))
            self.assertEqual(R.mul(
                Y, B), upper_matrix_tools.mul_matrix_and_basis_upper_t_circulant_matrix(R, s, k, Y))

    def test_mul_basis_lower_t_circulant_matrix_and_matrix(self):
        for i in range(10):
            R = tropical_algebra.MatrixSemiring(
                tropical_algebra.R_min_plus(), 10)
            s = random.randint(-100, 100)
            k = random.randint(0, 9)
            Y = matrix_tools.generate_random_matrix(R, -100, 100)
            B = lower_matrix_tools.generate_basis_lower_t_circulant_matrix(R, s, k)
            self.assertEqual(R.mul(
                B, Y), lower_matrix_tools.mul_basis_lower_t_circulant_matrix_and_matrix(R, s, k, Y))
            self.assertEqual(R.mul(
                Y, B), lower_matrix_tools.mul_matrix_and_basis_lower_t_circulant_matrix(R, s, k, Y))


if __name__ == "__main__":
    unittest.main()
