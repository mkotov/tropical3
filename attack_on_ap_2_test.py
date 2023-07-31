"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
An attack on Protocol 2 from B. Amutha and R. Perumal, Public key exchange protocols based on tropical lower circulant and anti circulant matrices, AIMS Mathematics, 8(7): 17307–17334.
"""

import tropical_algebra
import matrix_tools
import anti_matrix_tools
import attack
import test_tools
import random
import argparse


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        R = instance_params["ring"]
        s = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        t = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        p = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])

        Y = matrix_tools.generate_random_matrix(
            R, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P1 = anti_matrix_tools.generate_random_anti_t_p_circulant_matrix(
            R, s, p, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q1 = anti_matrix_tools.generate_random_anti_t_p_circulant_matrix(
            R, t, p, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P2 = anti_matrix_tools.generate_random_anti_t_p_circulant_matrix(
            R, s, p, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q2 = anti_matrix_tools.generate_random_anti_t_p_circulant_matrix(
            R, t, p, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
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
            "Kb": Kb,
            "K": KA
        }

    def run_attack(attack_params, instance):
        def compute_base_element(i, j):
            return matrix_tools.minus_matrix_from_matrix(
                anti_matrix_tools.mul_matrix_and_basis_anti_t_p_circulant_matrix(
                    attack_params["ring"], instance["t"], instance["p"],
                    anti_matrix_tools.mul_basis_anti_t_p_circulant_matrix_and_matrix(
                        attack_params["ring"], instance["s"], instance["p"], instance["Y"])),
                instance["Ka"])

        return attack.apply_attack(
            1,
            1,
            compute_base_element,
            bounds=(attack_params["min_matrix_elem"], attack_params["max_matrix_elem"]))

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(
            instance_params["ring"].mul(
                anti_matrix_tools.generate_anti_t_p_circulant_matrix(
                    instance_params["ring"], result[0], instance["s"]),
                instance["Kb"]),
            anti_matrix_tools.generate_anti_t_p_circulant_matrix(instance_params["ring"], result[1], instance["t"]))
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
        help="Lower bound to generate elements of matrices",
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
                          },
                          {
                              "ring": R,
                              "min_matrix_elem": args.min_matrix_elem,
                              "max_matrix_elem": args.max_matrix_elem,
                          },
                          args.count, args.timeout)
