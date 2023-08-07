"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
"""

import unittest
import attack_on_ap_1
import tropical_algebra
import matrix_tools


class TestTropicalAlgebra(unittest.TestCase):
    def test_attack(self):
        R = tropical_algebra.MatrixSemiring(
            tropical_algebra.R_min_plus(), 3)
        attack_params = {
            "ring": R,
            "min_matrix_elem": -2**16,
            "max_matrix_elem": 2**16,
        }

        s = -154
        t = 1797

        Y = [
            [-601, -615, 54332],
            [-554, 98, 45],
            [4325, 65, 3232],
        ]

        Ka = [
            [-1032, -1436, -1450],
            [-985, -1389, -1261],
            [-1052, -1456, -1470]
        ]

        Kb = [
            [-1273, -621, -674],
            [-1336, -1350, -51],
            [-1474, -1488, -690]
        ]

        instance = {
            "Y": Y,
            "s": s,
            "t": t,
            "Ka": Ka,
            "Kb": Kb
        }

        K = [
            [-1704, -2108, -2051],
            [-1771, -2175, -2189],
            [-1905, -2309, -2323]
        ]

        result = attack_on_ap_1.run_attack(attack_params, instance)
        P = matrix_tools.generate_lower_t_circulant_matrix(R, result[0], s)
        Q = matrix_tools.generate_lower_t_circulant_matrix(R, result[1], t)

        self.assertEqual(Ka, R.mul(R.mul(P, Y), Q))
        self.assertEqual(K, R.mul(R.mul(P, Kb), Q))


if __name__ == "__main__":
    unittest.main()
