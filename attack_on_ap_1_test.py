"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023.
An attack on a protocol from B. Amutha and R. Perumal, Public key exchange protocols based on tropical lower circulant and anti circulant matrices, AIMS Mathematics, 8(7): 17307–17334.
"""

import tropical_algebra
import matrix_tools
import lower_matrix_tools
import attack
import test_tools
import random
import argparse


def perform_one_experiment(instance_params, attack_params):
    def generate_instance(instance_params):
        s = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        t = random.randint(
            instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])

        Y = matrix_tools.generate_random_matrix(
            instance_params["ring"], instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P1 = lower_matrix_tools.generate_random_lower_t_circulant_matrix(
            instance_params["ring"], s, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q1 = lower_matrix_tools.generate_random_lower_t_circulant_matrix(
            instance_params["ring"], t, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        P2 = lower_matrix_tools.generate_random_lower_t_circulant_matrix(
            instance_params["ring"], s, instance_params["min_matrix_elem"], instance_params["max_matrix_elem"])
        Q2 = lower_matrix_tools.generate_random_lower_t_circulant_matrix(
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
        cache_Bi_mul_Y = dict()

        def compute_base_element(i, j):
            if i not in cache_Bi_mul_Y:
                cache_Bi_mul_Y[i] = lower_matrix_tools.mul_basis_lower_t_circulant_matrix_and_matrix(
                    attack_params["ring"], instance["s"], i, instance["Y"])
            return matrix_tools.minus_matrix_from_matrix(
                lower_matrix_tools.mul_matrix_and_basis_lower_t_circulant_matrix(
                    attack_params["ring"], instance["t"], j, cache_Bi_mul_Y[i]),
                instance["Ka"])

        return attack.apply_attack(
            attack_params["ring"].size(),
            attack_params["ring"].size(),
            compute_base_element,
            bounds=(attack_params["min_matrix_elem"], attack_params["max_matrix_elem"]))

    def check_key(instance_params, instance, result):
        KC = instance_params["ring"].mul(
            instance_params["ring"].mul(
                lower_matrix_tools.generate_lower_t_circulant_matrix(
                    instance_params["ring"], result[0], instance["s"]),
                instance["Kb"]),
            lower_matrix_tools.generate_lower_t_circulant_matrix(instance_params["ring"], result[1], instance["t"]))
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
