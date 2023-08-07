"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
"""

import unittest
import attack_on_ap_2
import tropical_algebra
import matrix_tools


class TestTropicalAlgebra(unittest.TestCase):
    def test_attack(self):
        R = tropical_algebra.MatrixSemiring(
            tropical_algebra.R_min_plus(), 3)
        attack_params = {
            "ring": R,
            "min_matrix_elem": -2**30,
            "max_matrix_elem": 2**30,
        }

        s = -71
        t = -98876
        p = 2

        Y = [
            [1090000, -33434, 32434543],
            [23, -2251, 955543],
            [-55432, 32455, 34442],
        ]

        Ka = [
            [-168593, -168597, -146668],
            [-168666, -168670, -146670],
            [-168662, -168666, -146601]
        ]

        Kb = [
            [-154799, -154803, -132874],
            [-154872, -154876, -132876],
            [-154868, -154872, -132807]
        ]

        instance = {
            "Y": Y,
            "s": s,
            "t": t,
            "p": p,
            "Ka": Ka,
            "Kb": Kb
        }

        K = [
            [-268112, -268110, -268114],
            [-268108, -268106, -268110],
            [-268110, -268108, -268112]
        ]

        result = attack_on_ap_2.run_attack(attack_params, instance)

        P = matrix_tools.generate_anti_t_p_circulant_matrix(
            R, result[0][0], s, p)
        Q = matrix_tools.generate_anti_t_p_circulant_matrix(
            R, result[1][0], t, p)

        self.assertEqual(K, R.mul(R.mul(P, Kb), Q))


if __name__ == "__main__":
    unittest.main()
