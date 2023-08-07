"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.

An attack on Protocol 2 from B. Amutha and R. Perumal, Public key exchange protocols based on tropical lower circulant and anti circulant matrices, AIMS Mathematics, 8(7): 17307–17334.
"""

import tropical_algebra
import matrix_tools
from matrix_tools import generate_random_anti_t_p_circulant_matrix
from matrix_tools import generate_random_matrix
from matrix_tools import generate_basis_anti_t_p_circulant_matrix
from matrix_tools import generate_anti_t_p_circulant_matrix
from matrix_tools import subtract_matrix_from_matrix
import attack
import test_tools
from random import randint
import argparse


def generate_instance(instance_params):
    R = instance_params["ring"]
    mm = instance_params["min_matrix_elem"]
    mM = instance_params["max_matrix_elem"]
    sm = instance_params["min_matrix_param"]
    sM = instance_params["max_matrix_param"]
    pm = instance_params["min_matrix_step"]
    pM = instance_params["max_matrix_step"]

    s = randint(sm, sM)
    t = randint(sm, sM)
    p = randint(pm, pM)

    Y = generate_random_matrix(R, mm, mM)
    P1 = generate_random_anti_t_p_circulant_matrix(R, s, p, mm, mM)
    Q1 = generate_random_anti_t_p_circulant_matrix(R, t, p, mm, mM)
    P2 = generate_random_anti_t_p_circulant_matrix(R, s, p, mm, mM)
    Q2 = generate_random_anti_t_p_circulant_matrix(R, t, p, mm, mM)

    Ka = R.mul(P1, R.mul(Y, Q1))
    Kb = R.mul(P2, R.mul(Y, Q2))
    KA = R.mul(P1, R.mul(Kb, Q1))
    KB = R.mul(P2, R.mul(Ka, Q2))

    if KA != KB:
        return None

    return {
        "Y": Y,
        "s": s,
        "p": p,
        "t": t,
        "Ka": Ka,
        "Kb": Kb
    }, KA


def run_attack(attack_params, instance):
    R = attack_params["ring"]
    mm = attack_params["min_matrix_elem"]
    mM = attack_params["max_matrix_elem"]

    s = instance["s"]
    t = instance["t"]
    p = instance["p"]
    Y = instance["Y"]
    Ka = instance["Ka"]

    def compute_base_element(i, j):
        B1 = generate_basis_anti_t_p_circulant_matrix(R, s, p)
        B2 = generate_basis_anti_t_p_circulant_matrix(R, t, p)

        return subtract_matrix_from_matrix(R.mul(R.mul(B1, Y), B2), Ka)

    return attack.apply_attack(1, 1, compute_base_element, bounds=(mm, mM))


def check_key(attack_params, instance, key, result):
    R = attack_params["ring"]

    s = instance["s"]
    t = instance["t"]
    p = instance["p"]
    Kb = instance["Kb"]

    P = generate_anti_t_p_circulant_matrix(R, result[0][0], s, p)
    Q = generate_anti_t_p_circulant_matrix(R, result[1][0], t, p)
    KC = R.mul(R.mul(P, Kb), Q)

    return key == KC


def perform_one_experiment(instance_params, attack_params):
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
    parser.add_argument(
        "--min_matrix_step",
        help="Lower bound to generate p",
        required=True,
        type=int
    )
    parser.add_argument(
        "--max_matrix_step",
        help="Lower bound to generate p",
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
                              "min_matrix_param": args.min_matrix_param,
                              "max_matrix_param": args.max_matrix_param,
                              "min_matrix_step": args.min_matrix_step,
                              "max_matrix_step": args.max_matrix_step,
                          },
                          {
                              "ring": R,
                              "min_matrix_elem": args.min_matrix_elem,
                              "max_matrix_elem": args.max_matrix_elem,
                          },
                          args.count, args.timeout)
