"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import time
import multiprocess


def perform_one_experiment(instance_params, attack_params, generate_instance, run_attack, check_key):
    """Runs one experiment: generates an instance, runs an attack, and check the obtained key."""

    instance, key = generate_instance(instance_params)
    if not instance or not key:
        return False

    result = run_attack(attack_params, instance)
    if not result:
        return False

    return check_key(attack_params, instance, key, result)


def test_suite(perform_one_experiment, instance_params, attack_params, number_of_tests, timeout):
    """Runs a set of tests."""
    st = time.time()
    ok = 0
    fl = 0
    fl_by_timeout = 0
    for i in range(1, number_of_tests + 1):
        q = multiprocess.Queue()
        p = multiprocess.Process(target=lambda q: q.put(
            perform_one_experiment(instance_params, attack_params)), args=(q,))
        p.start()
        p.join(timeout)
        if p.is_alive():
            print("TIMEOUT")
            p.terminate()
            p.join()
            fl_by_timeout += 1
        elif q.get():
            print("OK")
            ok += 1
        else:
            print("FAIL")
            fl += 1
    et = time.time()
    diff_time = et - st
    print("Total time: ", diff_time)
    print("Average time: ", diff_time / number_of_tests)
    print("OK: ", ok)
    print("FAIL: ", fl)
    print("FAIL (BY TIMEOUT): ", fl_by_timeout)
