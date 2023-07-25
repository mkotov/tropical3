"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023

An attack on the protocol from D. Grigoriev, V. Shpilrain, "Tropical cryptography", Comm. Algebra 43 (2014), 2624–2632, Section 2.
"""

import tropical_algebra
import matrix_tools
import attack
import test_tools
import argparse
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

        return attack.apply_attack(attack_params["max_poly_deg"] + 1, attack_params["max_poly_deg"] + 1, compute_base_element, bounds=(attack_params["min_poly_coef"], None))

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(instance_params["ring"].calc_poly(result[0], instance["A"]), instance_params["ring"].mul(
            instance["v"], instance_params["ring"].calc_poly(result[1], instance["B"])))
        return instance["K"] == KC

    return test_tools.perform_one_experiment(instance_params, attack_params, generate_instance, run_attack, check_key)


def get_arguments_parser():
    parser = argparse.ArgumentParser(
        description="The script to check the attack.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--count",
        help="Number of tests",
        required=True,
        type=int
    )
    parser.add_argument(
        "--size",
        help="Size of matrices",
        required=True,
        type=int
    )
    parser.add_argument(
        "--timeout",
        help="Timeout for each experiment",
        required=True,
        type=int
    )
    parser.add_argument(
        "--min_matrix_elem",
        help="Lower bound to generate elements of matrices",
        required=True,
        type=int
    )
    parser.add_argument(
        "--max_matrix_elem",
        help="Upper bound to generate elements of matrices",
        required=True,
        type=int
    )
    parser.add_argument(
        "--min_poly_deg",
        help="Lower bound to generate polynomial degrees",
        required=True,
        type=int
    )
    parser.add_argument(
        "--max_poly_deg",
        help="Upper bound to generate polynomial degrees",
        required=True,
        type=int
    )
    parser.add_argument(
        "--min_poly_coef",
        help="Lower bound to generate polynomial coefficients",
        required=True,
        type=int
    )
    parser.add_argument(
        "--max_poly_coef",
        help="Upper bound to generate polynomial coefficients",
        required=True,
        type=int
    )
    return parser


if __name__ == "__main__":
    args = get_arguments_parser().parse_args()
    R = tropical_algebra.MatrixSemiring(
        tropical_algebra.R_min_plus(), args.size)
    test_tools.test_suite(perform_one_experiment,
                          {
                              "ring": R,
                              "min_matrix_elem": args.min_matrix_elem,
                              "max_matrix_elem": args.max_matrix_elem,
                              "min_poly_deg": args.min_poly_deg,
                              "max_poly_deg": args.max_poly_deg,
                              "min_poly_coef": args.min_poly_coef,
                              "max_poly_coef": args.max_poly_coef},
                          {
                              "ring": R,
                              "max_poly_deg": args.max_poly_deg,
                              "min_poly_coef": args.min_poly_coef,
                              "max_poly_coef": args.max_poly_coef
                          },
                          args.count, args.timeout)
