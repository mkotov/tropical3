"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.

An attack on Protocol 1 from Huang, Li, and Deng, "Public-Key Cryptography Based on Tropical Circular Matrices", Applied Sciences, 12.15 (2022): 7401.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools
import random


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        s = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        t = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])

        Y = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P1 = matrix_tools.generate_random_upper_t_circulant_matrix(
            instance_params["ring"], s, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q1 = matrix_tools.generate_random_upper_t_circulant_matrix(
            instance_params["ring"], t, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P2 = matrix_tools.generate_random_upper_t_circulant_matrix(
            instance_params["ring"], s, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q2 = matrix_tools.generate_random_upper_t_circulant_matrix(
            instance_params["ring"], t, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Ka = instance_params["ring"].mul(
            P1, instance_params["ring"].mul(Y, Q1))
        Kb = instance_params["ring"].mul(
            P2, instance_params["ring"].mul(Y, Q2))
        KA = instance_params["ring"].mul(
            P1, instance_params["ring"].mul(Kb, Q1))
        KB = instance_params["ring"].mul(
            P2, instance_params["ring"].mul(Ka, Q2))
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

    def run_attack(attack_params, instance):
        def compute_base_element(i, j):
            return matrix_tools.minus_matrix_from_matrix(
                attack_params["ring"].mul(attack_params["ring"].mul(
                    matrix_tools.generate_basis_upper_t_circulant_matrix(
                        attack_params["ring"], instance["s"], i),
                    instance["Y"]),
                    matrix_tools.generate_basis_upper_t_circulant_matrix(attack_params["ring"], instance["t"], j)),
                instance["Ka"])

        def heuristics_to_sort(a):
            s = matrix_tools.spectrum(a, 0)
            t = matrix_tools.spectrum(a, 1)
            return (-len(a), -matrix_tools.number_of_indexes(a, 0) * matrix_tools.number_of_indexes(a, 1), s[0] * t[0], s, t)

        return attack.apply_attack(
            attack_params["ring"].size(),
            attack_params["ring"].size(),
            compute_base_element,
            heuristics_to_sort,
            bounds=(attack_params["min_matrix_elem"], attack_params["max_matrix_elem"]))

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(
            instance_params["ring"].mul(
                matrix_tools.generate_upper_t_circulant_matrix(
                    instance_params["ring"], result[0], instance["s"]),
                instance["Kb"]),
            matrix_tools.generate_upper_t_circulant_matrix(instance_params["ring"], result[1], instance["t"]))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(instance_params, attack_params, generate_instance, run_attack, check_key)


R = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 15)
test_tools.test_suite(perform_one_experiment,
                      {
                          "ring": R,
                          "min_matrix_elem": -10000,
                          "max_matrix_elem": 10000
                      },
                      {
                          "ring": R,
                          "min_matrix_elem": -10000,
                          "max_matrix_elem": 10000
                      },
                      10, 60)
