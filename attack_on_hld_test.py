"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
"""

import unittest
import attack_on_hld
import tropical_algebra
import matrix_tools


class TestTropicalAlgebra(unittest.TestCase):
    def test_attack(self):
        R = tropical_algebra.MatrixSemiring(
            tropical_algebra.R_min_plus(), 5)
        attack_params = {
            "ring": R,
            "min_matrix_elem": 0,
            "max_matrix_elem": 2**15,
        }

        s = t = 9361

        Y = [
            [8630, 29391, 21921, 18968, 25014],
            [15306, 5461, 18973, 800, 1786],
            [7986, 27430, 22510, 11233, 30900],
            [2398, 6071, 25269, 27186, 4328],
            [18306, 10527, 16873, 11565, 9569]
        ]

        Ka = [
            [26578, 19555, 38342, 32846, 29893],
            [3350, 25959, 16386, 21160, 11725],
            [24783, 18911, 30607, 33184, 22158],
            [5892, 13323, 16996, 23702, 26279],
            [11133, 29231, 21452, 27798, 21563]
        ]

        Kb = [
            [18245, 27756, 29434, 23095, 24081],
            [18102, 15076, 16754, 10415, 11401],
            [17601, 18918, 20596, 14257, 15243],
            [12013, 15686, 31029, 20282, 13943],
            [15855, 19528, 26488, 21180, 17785],
        ]

        instance = {
            "Y": Y,
            "s": s,
            "t": t,
            "Ka": Ka,
            "Kb": Kb
        }

        K = [
            [25645, 29170, 38681, 40359, 34020],
            [12965, 29027, 26001, 27679, 21340],
            [16807, 28526, 29843, 31521, 25182],
            [15507, 22938, 26611, 33317, 31207],
            [19349, 26780, 30453, 37159, 31178]
        ]

        result = attack_on_hld.run_attack(attack_params, instance)
        P = matrix_tools.generate_upper_t_circulant_matrix(R, result[0], s)
        Q = matrix_tools.generate_upper_t_circulant_matrix(R, result[1], t)

        self.assertEqual(Ka, R.mul(R.mul(P, Y), Q))
        self.assertEqual(K, R.mul(R.mul(P, Kb), Q))


if __name__ == "__main__":
    unittest.main()
