"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023

An attack on the protocol from D. Grigoriev, V. Shpilrain, "Tropical cryptography", Comm. Algebra 43 (2014), 2624–2632, Section 2.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        """Generates an instance of the protocol."""
        A = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        B = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        p1 = matrix_tools.generate_random_polynomial(
            instance_params["max_poly_deg"], instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        p2 = matrix_tools.generate_random_polynomial(
            instance_params["max_poly_deg"], instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        u = instance_params["ring"].mul(instance_params["ring"].calc_poly(
            p1, A), instance_params["ring"].calc_poly(p2, B))
        q1 = matrix_tools.generate_random_polynomial(
            instance_params["max_poly_deg"], instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        q2 = matrix_tools.generate_random_polynomial(
            instance_params["max_poly_deg"], instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        v = instance_params["ring"].mul(instance_params["ring"].calc_poly(
            q1, A), instance_params["ring"].calc_poly(q2, B))
        KA = instance_params["ring"].mul(instance_params["ring"].calc_poly(
            p1, A), instance_params["ring"].mul(v, instance_params["ring"].calc_poly(p2, B)))
        KB = instance_params["ring"].mul(instance_params["ring"].calc_poly(
            q1, A), instance_params["ring"].mul(u, instance_params["ring"].calc_poly(q2, B)))

        if KA != KB:
            return None
        return {
            "A": A,
            "B": B,
            "u": u,
            "v": v,
            "K": KA
        }

    def run_attack(attack_params, instance):
        def compute_base_element(i, j):
            return matrix_tools.minus_matrix_from_matrix(attack_params["ring"].mul(attack_params["ring"].pwr(instance["A"], i), attack_params["ring"].pwr(instance["B"], j)), instance["u"])

        def heuristics_to_sort(a):
            return (-len(a), -matrix_tools.number_of_indexes(a, 0) * matrix_tools.number_of_indexes(a, 1))

        return attack.apply_attack(attack_params["max_poly_deg"] + 1, attack_params["max_poly_deg"] + 1, compute_base_element, heuristics_to_sort, bounds=(attack_params["min_poly_coef"], None))

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(instance_params["ring"].calc_poly(result[0], instance["A"]), instance_params["ring"].mul(
            instance["v"], instance_params["ring"].calc_poly(result[1], instance["B"])))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(instance_params, attack_params, generate_instance, run_attack, check_key)


R = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 10)
test_tools.test_suite(perform_one_experiment,
                      {
                          "ring": R,
                          "min_matrix_elem": -100,
                          "max_matrix_elem": 100,
                          "max_poly_deg": 10,
                          "min_poly_coef": -100,
                          "max_poly_coef": 100},
                      {
                          "ring": R,
                          "max_poly_deg": 10,
                          "min_poly_coef": -100,
                          "max_poly_coef": 100
                      },
                      10, 60)
