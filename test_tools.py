"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import time

def perform_one_experiment(params, generate_instance, run_attack, check_key):
    """Runs one experiment: generates an instance, runs an attack, and check the obtained key."""

    instance = generate_instance(params)
    if not instance:
        return False

    result = run_attack(params, instance)
    if not result:
        return False

    return check_key(params, instance, result)


def test_suite(perform_one_experiment, params, number_of_tests):
    """Runs a set of tests."""
    st = time.time()
    ok = 0
    fl = 0
    for i in range(1, number_of_tests + 1):
        if perform_one_experiment(params):
            print("OK")
            ok += 1
        else:
            print("FAIL")
            fl += 1
    et = time.time()
    print(et - st)
    print("OK: ", ok)
    print("FAIL: ", fl)
