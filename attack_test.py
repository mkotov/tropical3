"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import time
import tropical_algebra
import matrix_tools
import attack


def generate_instance(R, mm, mM, D, pm, pM):
    """Generates an instance of the protocol."""
    A = matrix_tools.generate_random_matrix(R.size(), mm, mM)
    B = matrix_tools.generate_random_matrix(R.size(), mm, mM)
    p1 = matrix_tools.generate_random_polynomial(D, pm, pM)
    p2 = matrix_tools.generate_random_polynomial(D, pm, pM)
    u = R.mul(R.calc_poly(p1, A), R.calc_poly(p2, B))
    q1 = matrix_tools.generate_random_polynomial(D, pm, pM)
    q2 = matrix_tools.generate_random_polynomial(D, pm, pM)
    v = R.mul(R.calc_poly(q1, A), R.calc_poly(q2, B))
    KA = R.mul(R.calc_poly(p1, A), R.mul(v, R.calc_poly(p2, B)))
    KB = R.mul(R.calc_poly(q1, A), R.mul(u, R.calc_poly(q2, B)))

    if KA != KB:
        return None

    return {
        "A": A,
        "B": B,
        "u": u,
        "v": v,
        "K": KA
    }


def test_attack(n, mm, mM, D, pm, pM):
    """Runs our attack."""
    R = tropical_algebra.MatrixSemiring(tropical_algebra.R_min_plus(), n)

    I = generate_instance(R, mm, mM, D, pm, pM)
    if not I:
        return False

    attack_result = attack.apply_attack(R, I["A"], I["B"], I["u"], D, pm)
    if not attack_result:
        return False

    KC = R.mul(R.calc_poly(attack_result[0], I["A"]), R.mul(
        I["v"], R.calc_poly(attack_result[1], I["B"])))

    if I["K"] != KC:
        return False

    return True


def test_suite(n, mm, mM, D, pm, pM, numberOfTests):
    """Runs a set of tests."""
    st = time.time()
    ok = 0
    fl = 0
    for i in range(1, numberOfTests + 1):
        if test_attack(n, mm, mM, D, pm, pM):
            print("OK")
            ok += 1
        else:
            print("FAIL")
            fl += 1
    et = time.time()
    print(et - st)
    print("OK: ", ok)
    print("FAIL: ", fl)


test_suite(5, -10**2, 10**2, 10, -10**2, 10**2, 10)
