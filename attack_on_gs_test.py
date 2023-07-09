"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023

An attack on the protocol from D. Grigoriev, V. Shpilrain, "Tropical cryptography", Comm. Algebra 43 (2014), 2624–2632, Section 2.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools


def perform_one_experiment(params):
    def generate_instance(params):
        """Generates an instance of the protocol."""
        A = matrix_tools.generate_random_matrix(
            params["R"], params["mm"], params["mM"])
        B = matrix_tools.generate_random_matrix(
            params["R"], params["mm"], params["mM"])
        p1 = matrix_tools.generate_random_polynomial(
            params["D"], params["pm"], params["pM"])
        p2 = matrix_tools.generate_random_polynomial(
            params["D"], params["pm"], params["pM"])
        u = params["R"].mul(params["R"].calc_poly(p1, A),
                            params["R"].calc_poly(p2, B))
        q1 = matrix_tools.generate_random_polynomial(
            params["D"], params["pm"], params["pM"])
        q2 = matrix_tools.generate_random_polynomial(
            params["D"], params["pm"], params["pM"])
        v = params["R"].mul(params["R"].calc_poly(q1, A),
                            params["R"].calc_poly(q2, B))
        KA = params["R"].mul(params["R"].calc_poly(
            p1, A), params["R"].mul(v, params["R"].calc_poly(p2, B)))
        KB = params["R"].mul(params["R"].calc_poly(
            q1, A), params["R"].mul(u, params["R"].calc_poly(q2, B)))

        if KA != KB:
            return None

        return {
            "A": A,
            "B": B,
            "u": u,
            "v": v,
            "K": KA
        }

    def run_attack(params, instance):
        def compute_base_element(i, j):
            return matrix_tools.minus_matrix_from_matrix(
                params["R"].mul(params["R"].pwr(instance["A"], i),
                                params["R"].pwr(instance["B"], j)),
                params["R"].mul_by_coef(-2 * params["pm"], instance["u"]))

        def extract_solution(r1, r2):
            return [[r + params["pm"] for r in r1], [r + params["pm"] for r in r2]]

        def heuristics_to_sort(a):
            return (-len(a), -matrix_tools.number_of_indexes(a, 0) * matrix_tools.number_of_indexes(a, 1))

        return attack.apply_attack(params["D"], compute_base_element, extract_solution, heuristics_to_sort)

    def check_key(params, instance, result):
        KC = params["R"].mul(params["R"].calc_poly(result[0], instance["A"]), params["R"].mul(
            instance["v"], params["R"].calc_poly(result[1], instance["B"])))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(params, generate_instance, run_attack, check_key)


test_tools.test_suite(perform_one_experiment, {
    "R":  tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), 10),
    "mm": -100,
    "mM": 100,
    "D": 10,
    "pm": -100,
    "pM": 100},
    10)
