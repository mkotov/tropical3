"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023

An attack on the second step of the protocol from Durcheva, M. I. "TrES: Tropical Encryption Scheme Based on Double Key Exchange." European J. IT and CS, 2(4), 11–17.
"""

import tropical_algebra
from matrix_tools import generate_random_matrix
from matrix_tools import generate_random_polynomial
from matrix_tools import get_first_repeated
from matrix_tools import subtract_matrix_from_matrix
import attack
import test_tools
from random import randint
import argparse


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        """Generates an instance of the protocol."""
        R = instance_params["ring"]
        dm = instance_params["min_poly_deg"]
        dM = instance_params["max_poly_deg"]
        mm = instance_params["min_matrix_elem"]
        mM = instance_params["max_matrix_elem"]
        cm = instance_params["min_poly_coef"]
        cM = instance_params["max_poly_coef"]

        D1 = randint(dm, dM)
        D2 = randint(dm, dM)
        D3 = randint(dm, dM)
        D4 = randint(dm, dM)
        A = generate_random_matrix(R, mm, mM)
        B = generate_random_matrix(R, mm, mM)
        L = generate_random_matrix(R, mm, mM)
        p1 = generate_random_polynomial(D1, cm, cM)
        p2 = generate_random_polynomial(D2, cm, cM)
        q1 = generate_random_polynomial(D3, cm, cM)
        q2 = generate_random_polynomial(D4, cm, cM)
        u = R.mul(R.mul(R.calc_poly(p1, A), L), R.calc_poly(p2, B))
        v = R.mul(R.mul(R.calc_poly(q1, A), L), R.calc_poly(q2, B))
        KA = R.mul(R.calc_poly(p1, A), R.mul(v, R.calc_poly(p2, B)))
        KB = R.mul(R.calc_poly(q1, A), R.mul(u, R.calc_poly(q2, B)))

        if KA != KB:
            return None

        return {
            "A": A,
            "B": B,
            "L": L,
            "u": u,
            "v": v
        }, KA

    def run_attack(attack_params, instance):
        R = attack_params["ring"]
        db = attack_params["poly_deg_bound"]
        A = instance["A"]
        B = instance["B"]
        u = instance["u"]
        L = instance["L"]

        def get_degree_bound(A):
            d = get_first_repeated(R, A, db)
            if d:
                return d
            return db

        d1 = get_degree_bound(A)
        d2 = get_degree_bound(B)

        cache_Ai = dict()
        cache_Bj = dict()

        cache_Ai[0] = R.one()
        cache_Ai[1] = A
        for i in range(2, d1 + 1):
            cache_Ai[i] = R.mul(cache_Ai[i - 1], A)

        cache_AiL = dict()

        cache_Bj[0] = R.one()
        cache_Bj[1] = B
        for j in range(2, d2 + 1):
            cache_Bj[j] = R.mul(cache_Bj[j - 1], B)

        def compute_base_element(i, j):
            if i not in cache_AiL:
                cache_AiL[i] = R.mul(cache_Ai[i], L)
            return subtract_matrix_from_matrix(R.mul(cache_AiL[i], cache_Bj[j]), u)

        return attack.apply_attack(d1 + 1, d2 + 1, compute_base_element)

    def check_key(attack_params, instance, key, result):
        R = attack_params["ring"]
        A = instance["A"]
        B = instance["B"]
        v = instance["v"]

        KC = R.mul(R.calc_poly(result[0], A),
                   R.mul(v, R.calc_poly(result[1], B)))

        return key == KC

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
    parser.add_argument(
        "--poly_deg_bound",
        help="Upper bound for polynomial degrees, this is a parameter of the attack",
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
                              "poly_deg_bound": args.poly_deg_bound
                          },
                          args.count, args.timeout)
