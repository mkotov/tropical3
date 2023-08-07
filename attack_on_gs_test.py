"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
"""

import unittest
import attack_on_gs
import tropical_algebra
import matrix_tools


class TestTropicalAlgebra(unittest.TestCase):
    def test_attack(self):
        R = tropical_algebra.MatrixSemiring(
            tropical_algebra.R_min_plus(), 3)

        A = [
            [765, 233, 999],
            [387, 458, 277],
            [115, 676, 189]
        ]
        B = [
            [658, 858, 279],
            [215, 369, 570],
            [834, 772, 468]
        ]

        u = [
            [-369, -215, -459],
            [-589, -435, -679],
            [-677, -523, -767]
        ]

        v = [
            [308, 99, 198],
            [80, 308, 143],
            [-19, 36, 55]
        ]

        K = [
            [-490, -336, -580],
            [-534, -380, -624],
            [-622, -468, -712]
        ]

        attack_params = {
            "ring": R,
            "max_poly_deg": 10,
            "min_poly_coef": -1000,
            "max_poly_coef": 1000
        }

        instance = {
            "A": A,
            "B": B,
            "u": u,
            "v": v
        }

        result = attack_on_gs.run_attack(attack_params, instance)

        self.assertEqual(K, R.mul(R.calc_poly(
            result[0], A), R.mul(v, R.calc_poly(result[1], B))))


if __name__ == "__main__":
    unittest.main()
