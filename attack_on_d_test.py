"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023

An attack on the protocol from Durcheva, M. I. "TrES: Tropical Encryption Scheme Based on Double Key Exchange." European J. IT and CS, 2(4), 11–17.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools
import random


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        """Generates an instance of the protocol."""
        D1 = random.randint(
            instance_params["min_poly_deg"], instance_params["max_poly_deg"])
        D2 = random.randint(
            instance_params["min_poly_deg"], instance_params["max_poly_deg"])
        D3 = random.randint(
            instance_params["min_poly_deg"], instance_params["max_poly_deg"])
        D4 = random.randint(
            instance_params["min_poly_deg"], instance_params["max_poly_deg"])
        A = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        B = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        p1 = matrix_tools.generate_random_polynomial(
            D1, instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        p2 = matrix_tools.generate_random_polynomial(
            D2, instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        u = instance_params["ring"].mul(instance_params["ring"].calc_poly(
            p1, A), instance_params["ring"].calc_poly(p2, B))
        q1 = matrix_tools.generate_random_polynomial(
            D3, instance_params["min_poly_coef"], instance_params["max_poly_coef"])
        q2 = matrix_tools.generate_random_polynomial(
            D4, instance_params["min_poly_coef"], instance_params["max_poly_coef"])
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

        def get_degree_bound(A):
            d = matrix_tools.get_first_repeated(
                attack_params["ring"], A, attack_params["poly_deg_bound"])
            if d:
                return d
            return attack_params["poly_deg_bound"]

        d1 = get_degree_bound(instance["A"])
        d2 = get_degree_bound(instance["B"])
        return attack.apply_attack(d1 + 1, d2 + 1, compute_base_element, heuristics_to_sort)

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(instance_params["ring"].calc_poly(result[0], instance["A"]), instance_params["ring"].mul(
            instance["v"], instance_params["ring"].calc_poly(result[1], instance["B"])))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(instance_params, attack_params, generate_instance, run_attack, check_key)


R = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 10)
test_tools.test_suite(perform_one_experiment,
                      {
                          "ring": R,
                          "min_matrix_elem": -100000,
                          "max_matrix_elem": 100000,
                          "min_poly_deg": 5,
                          "max_poly_deg": 25,
                          "min_poly_coef": -100000,
                          "max_poly_coef": 100000},
                      {
                          "ring": R,
                          "poly_deg_bound": 25},
                      10)
