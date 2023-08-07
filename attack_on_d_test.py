"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
"""

import unittest
import attack_on_d
import tropical_algebra
import matrix_tools


class TestTropicalAlgebra(unittest.TestCase):
    def test_attack(self):
        R = tropical_algebra.MatrixSemiring(
            tropical_algebra.R_min_plus(), 3)

        M = [
            [15, 63, 97],
            [57, 60, 83],
            [48, 100, 26]
        ]
        L = [
            [12, 62, 3],
            [49, 55, 77],
            [97, 98, 0]
        ]
        u = [
            [-83, -33, -92],
            [-41, 9, -50],
            [-50, 0, -73]
        ]
        v = [
            [-28, 22, -37],
            [14, 23, 5],
            [5, 55, -32]
        ]
        K = [
            [-123, -73, -132],
            [-81, -31, -90],
            [-90, -40, -105]
        ]

        attack_params = {
            "ring": R,
            "poly_deg_bound": 20
        }

        instance = {
            "M": M,
            "L": L,
            "u": u,
            "v": v
        }

        result = attack_on_d.run_attack(attack_params, instance)

        self.assertEqual(K,  R.mul(R.calc_poly(
            result[0], M), R.mul(v, R.calc_poly(result[1], M))))


if __name__ == "__main__":
    unittest.main()
