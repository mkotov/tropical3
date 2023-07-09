"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.

An attack on Protocol 1 from Huang, Li, and Deng, "Public-Key Cryptography Based on Tropical Circular Matrices", Applied Sciences, 12.15 (2022): 7401.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools
import random


def perform_one_experiment(params):
    def generate_instance(params):
        s = random.randint(params["a_min"], params["a_max"])
        t = random.randint(params["a_min"], params["a_max"])

        Y = matrix_tools.generate_random_matrix(
            params["R"], params["a_min"], params["a_max"])
        P1 = matrix_tools.generate_random_upper_t_circulant_matrix(
            params["R"], s, params["a_min"], params["a_max"])
        Q1 = matrix_tools.generate_random_upper_t_circulant_matrix(
            params["R"], t, params["a_min"], params["a_max"])
        P2 = matrix_tools.generate_random_upper_t_circulant_matrix(
            params["R"], s, params["a_min"], params["a_max"])
        Q2 = matrix_tools.generate_random_upper_t_circulant_matrix(
            params["R"], t, params["a_min"], params["a_max"])
        Ka = params["R"].mul(P1, params["R"].mul(Y, Q1))
        Kb = params["R"].mul(P2, params["R"].mul(Y, Q2))
        KA = params["R"].mul(P1, params["R"].mul(Kb, Q1))
        KB = params["R"].mul(P2, params["R"].mul(Ka, Q2))
        if KA != KB:
            return None

        return {
            "Y": Y,
            "s": s,
            "t": t,
            "Ka": Ka,
            "Kb": Kb,
            "K": KA
        }

    def run_attack(params, instance):
        def compute_base_element(i, j):
            return matrix_tools.minus_matrix_from_matrix(
                params["R"].mul(params["R"].mul(
                    matrix_tools.generate_basis_upper_t_circulant_matrix(
                        params["R"], instance["s"], i),
                    instance["Y"]),
                    matrix_tools.generate_basis_upper_t_circulant_matrix(params["R"], instance["t"], j)),
                params["R"].mul_by_coef(-2 * params["a_min"], instance["Ka"]))

        def extract_solution(r1, r2):
            return [[r + params["a_min"] for r in r1], [r + params["a_min"] for r in r2]]

        def heuristics_to_sort(a):
            s = matrix_tools.spectrum(a, 0)
            t = matrix_tools.spectrum(a, 1)
            return (-len(a), -matrix_tools.number_of_indexes(a, 0) * matrix_tools.number_of_indexes(a, 1), s[0] * t[0], s, t)

        return attack.apply_attack(params["R"].size() - 1, compute_base_element, extract_solution, heuristics_to_sort)

    def check_key(params, instance, result):
        KC = params["R"].mul(params["R"].mul(matrix_tools.generate_upper_t_circulant_matrix(params["R"], result[0], instance["s"]),
                         instance["Kb"]),
                   matrix_tools.generate_upper_t_circulant_matrix(params["R"], result[1], instance["t"]))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(params, generate_instance, run_attack, check_key)


test_tools.test_suite(perform_one_experiment, {
    "R":  tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 10),
    "a_min": 1,
    "a_max": 10000},
    25)
