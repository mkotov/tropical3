"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
An attack on Protocol 1 from B. Amutha and R. Perumal, Public key exchange protocols based on tropical lower circulant and anti circulant matrices, AIMS Mathematics, 8(7): 17307–17334.
"""

from tropical_algebra import MatrixSemiring
from tropical_algebra import R_min_plus
from matrix_tools import generate_random_lower_t_circulant_matrix
from matrix_tools import generate_random_matrix
from matrix_tools import mul_basis_lower_t_circulant_matrix_and_matrix
from matrix_tools import generate_lower_t_circulant_matrix
from matrix_tools import subtract_matrix_from_matrix
from matrix_tools import mul_matrix_and_basis_lower_t_circulant_matrix
import attack
import test_tools
from random import randint
import argparse

def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        R = instance_params["ring"]
        mm = instance_params["min_matrix_elem"]
        mM = instance_params["max_matrix_elem"]
        sm = instance_params["min_matrix_param"]
        sM = instance_params["max_matrix_param"]

        s = randint(sm, sM)
        t = randint(sm, sM)
        Y = generate_random_matrix(R, mm, mM)
        P1 = generate_random_lower_t_circulant_matrix(R, s, mm, mM)
        Q1 = generate_random_lower_t_circulant_matrix(R, t, mm, mM)
        P2 = generate_random_lower_t_circulant_matrix(R, s, mm, mM)
        Q2 = generate_random_lower_t_circulant_matrix(R, t, mm, mM)

        Ka = R.mul(P1, R.mul(Y, Q1))
        Kb = R.mul(P2, R.mul(Y, Q2))
        KA = R.mul(P1, R.mul(Kb, Q1))
        KB = R.mul(P2, R.mul(Ka, Q2))

        if KA != KB:
            return None

        return {
            "Y": Y,
            "s": s,
            "t": t,
            "Ka": Ka,
            "Kb": Kb,
        }, KA

    def run_attack(attack_params, instance):
        R = attack_params["ring"]
        mm = attack_params["min_matrix_elem"]
        mM = attack_params["max_matrix_elem"]
        n = R.size()

        s = instance["s"]
        t = instance["t"]
        Y = instance["Y"]
        Ka = instance["Ka"]

        cache_Bi_mul_Y = dict()

        def compute_base_element(i, j):
            if i not in cache_Bi_mul_Y:
                cache_Bi_mul_Y[i] = mul_basis_lower_t_circulant_matrix_and_matrix(
                    R, s, i, Y)
            return subtract_matrix_from_matrix(mul_matrix_and_basis_lower_t_circulant_matrix(R, t, j, cache_Bi_mul_Y[i]), Ka)

        return attack.apply_attack(n, n, compute_base_element, bounds=(mm, mM))

    def check_key(attack_params, instance, key, result):
        R = attack_params["ring"]

        s = instance["s"]
        t = instance["t"]
        Kb = instance["Kb"]

        P = generate_lower_t_circulant_matrix(R, result[0], s)
        Q = generate_lower_t_circulant_matrix(R, result[1], t)

        K = R.mul(R.mul(P, Kb), Q)

        return key == K

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
        help="Lower bound to generate elements of matrices",
        required=True,
        type=int
    )
    parser.add_argument(
        "--min_matrix_param",
        help="Lower bound to generate s and t",
        required=True,
        type=int
    )
    parser.add_argument(
        "--max_matrix_param",
        help="Lower bound to generate s and t",
        required=True,
        type=int
    )

    return parser


if __name__ == "__main__":
    args = get_arguments_parser().parse_args()
    R = MatrixSemiring(R_min_plus(), args.size)

    test_tools.test_suite(perform_one_experiment,
                          {
                              "ring": R,
                              "min_matrix_elem": args.min_matrix_elem,
                              "max_matrix_elem": args.max_matrix_elem,
                              "min_matrix_param": args.min_matrix_param,
                              "max_matrix_param": args.max_matrix_param,
                          },
                          {
                              "ring": R,
                              "min_matrix_elem": args.min_matrix_elem,
                              "max_matrix_elem": args.max_matrix_elem,
                          },
                          args.count,
                          args.timeout)
