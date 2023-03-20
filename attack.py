import random
import math
import time
import functools
from itertools import product


def generate_random_matrix(n, mm, mM):
    """Generates a random matrix A in Mat(ZZ, n), a_ij in [mm, mM]."""

    return [[random.randint(mm, mM) for j in range(n)] for i in range(n)]


def generate_random_polynomial(D, pm, pM):
    """Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM]."""

    return [[i, random.randint(pm, pM)] for i in range(random.randint(1, D) + 1)]


def product_of_two_scalar_min_plus(a, b):
    """Returns a \otimes b."""

    if math.isinf(a) or math.isinf(b):
        return float('inf')
    else:
        return a + b


def product_of_two_matrices_min_plus(A, B):
    """Retruns A \otimes b."""

    n = len(A)
    result = []
    for i in range(n):
        result.append([])
        for j in range(n):
            result[i].append(
                min([product_of_two_scalar_min_plus(A[i][k], B[k][j]) for k in range(n)]))
    return result


def zero_matrix_min_plus(n):
    """Returns zero matrix of size n x n over min-plus algebra."""

    return [[float('inf') for j in range(n)] for i in range(n)]


def identity_matrix_min_plus(n):
    """Returns ident matrix of size nxn over min-plus algebra."""

    return [[0 if i == j else float('inf') for j in range(n)] for i in range(n)]


def power_of_matrix_min_plus(A, n):
    """Returns A^n."""

    if n == 0:
        return identity_matrix_min_plus(len(A))
    elif n == 1:
        return A
    else:
        return product_of_two_matrices_min_plus(A, power_of_matrix_min_plus(A, n - 1))


def sum_of_two_matrices_max_plus(A, B):
    """Returns A + B."""

    return [[min(x, y) for x, y in zip(row_a, row_b)] for row_a, row_b in zip(A, B)]


def product_of_scalar_and_matrix_min_plus(s, A):
    """Returns s \otimes A."""

    return [[product_of_two_scalar_min_plus(s, x) for x in row] for row in A]


def apply_polynomial_min_plus(p, A):
    """Returns p(A)."""

    if len(p) == 0:
        return zero_matrix_min_plus(len(A))
    return functools.reduce(sum_of_two_matrices_max_plus, [product_of_scalar_and_matrix_min_plus(c[1], power_of_matrix_min_plus(A, c[0])) for c in p])


def minus_matrix_from_matrix(A, B):
    """Returns A - B. Some elements of the matrices can be infinity."""
    return [[x - y if x != float('inf') and y != float('inf') else float('inf') if x == float('inf') else None for x, y in zip(row_a, row_b)] for row_a, row_b in zip(A, B)]


def get_minimum_of_matrix(A):
    """Returns minimum of a matrix A and the set of corresponding indexes."""

    best = {'val': A[0][0], 'inds': set()}
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] < best['val']:
                best = {'val': A[i][j], 'inds': {(i, j)}}
            elif A[i][j] == best['val']:
                best['inds'].add((i, j))
    return best


def make_simplex_matrix(F, D, h):
    """Retruns the simpex-table constructed by a set of covers and a fixed cover."""

    A = [[0] * (2 * (D + 1) + (D + 1)**2 + 1) for _ in range((D + 1)**2 + 1)]

    for i in range(D + 1):
        for j in range(D + 1):
            A[i * (D + 1) + j][i] = 1
            A[i * (D + 1) + j][(D + 1) + j] = 1

    for S in F:
        A[S['ijminval'][0]['i'] * (D + 1) + S['ijminval'][0]['j']
          ][2 * (D + 1) + (D + 1)**2] = -S['ijminval'][0]['val']
    for i in range((D + 1)**2):
        A[i][2 * (D + 1) + i] = -1
    for S in h:
        A[S[0] * (D + 1) + S[1]][2 * (D + 1) + S[0] * (D + 1) + S[1]] = 0
    for i in range(len(A[0])):
        S = 0
        for j in range(len(A) - 1):
            S += A[j][i]
        A[(D + 1)**2][i] = -S

    return A


def find_pivot_column(M):
    """Returns index of pivot column."""

    best = 0
    best_i = None
    for i in range(len(M[0]) - 1):
        t = M[len(M) - 1][i]
        if t < best:
            best = t
            best_i = i
    return best_i


def find_pivot_row(M, l):
    """Returns index of pivot row."""

    best_i = None
    best = float("inf")
    for i in range(len(M) - 1):
        a = M[i][l]
        b = M[i][len(M[i]) - 1]
        if a > 0 and b / a < best:
            best_i = i
            best = b / a
    return best_i


def recalc(M, k, l, Bs, Ns):
    """Performs recalculation of a simplex-table. Changes M, Bs and Ns."""

    t = Ns[l]
    Ns[l] = Bs[k]
    Bs[k] = t

    for i in range(len(M)):
        for j in range(len(M[0])):
            if i != k and j != l:
                M[i][j] = M[i][j] - M[i][l] * M[k][j] / M[k][l]

    for i in range(len(M)):
        if i != k:
            M[i][l] = -M[i][l] / M[k][l]

    for i in range(len(M[0])):
        if i != l:
            M[k][i] = M[k][i] / M[k][l]

    M[k][l] = 1 / M[k][l]


def write_down_solution(M, Bs):
    """Returns the solution constructed by a simplex-table."""
    result = [0] * (len(M) + len(M[0]))
    for i in range(len(Bs)):
        result[Bs[i]] = M[i][len(M[0]) - 1]

    return result


def apply_simplex(M):
    """Applies the simplex-method to a matrix M and returns the corresponding solution. If there is no solution, then returns None."""

    Bs = list(range(len(M[0]) - 1, len(M[0]) + len(M) - 2))
    Ns = list(range(len(M[0]) - 1))
    while True:
        l = find_pivot_column(M)
        if l == None:
            break
        k = find_pivot_row(M, l)
        if k == None:
            return None
        recalc(M, k, l, Bs, Ns)
    if M[-1][-1] != 0:
        return None
    return write_down_solution(M, Bs)


def number_of_indexes(S, i):
    """Returns the number of i-th indexes."""

    return len(set([c[i-1] for c in S]))


def rar(G):
    """Compresses a set of pair [a1, L1], [a2, L2], ..., [an, Ln] to the set of pair, where all first components are unique."""

    def find(H, S):
        for i in range(len(H)):
            if S['inds'] == H[i]['inds']:
                return i
        return None

    H = []
    for S in G:
        i = find(H, S)
        if i is None:
            H.append(S)
        else:
            H[i]['ijminval'].extend(S['ijminval'])
    return H


def unite_sets(ss):
    if len(ss) == 0:
        return set()
    else:
        return set.union(*ss)


def get_compressed_covers(F):
    """Returns the compressed set of covers. This set contains all the minimal covers."""

    if len(F) == 0:
        return [[]]
    Z = rar(F)
    M = list(filter(lambda S: len(S['inds'].difference(unite_sets(
        [T['inds'] for T in filter(lambda T: S != T, Z)]))) != 0, Z))
    N = unite_sets([S['inds'] for S in M])
    P = rar(list(filter(lambda S: len(S['inds']) != 0, [
            {'ijminval': S['ijminval'], 'inds': S['inds'].difference(N)} for S in Z])))
    if len(P) > 0:
        P.sort(key=lambda S: len(S['inds']), reverse=True)
        return [M + S for S in [[P[0]] +
                                S for S in get_compressed_covers(list(filter(lambda S: len(S['inds']) != 0, [{'ijminval': S['ijminval'], 'inds': S['inds'].difference(P[0]['inds'])} for S in P])))] +
                get_compressed_covers(P[1:])]

    return [M]


def apply_attack(A, B, u, D, pm):
    """Applies our attack. Returns two polynomials p' and q'."""

    def repack(ij, m):
        return {'ijminval': [{'i': ij[0], 'j': ij[1], 'val': m['val']}], 'inds': m['inds']}

    F = list(filter(lambda S: S['ijminval'][0]['val'] <= 0,
                    map(lambda ij: repack(ij, get_minimum_of_matrix(
                        minus_matrix_from_matrix(product_of_two_matrices_min_plus(power_of_matrix_min_plus(A, ij[0]), power_of_matrix_min_plus(B, ij[1])),
                                                 product_of_scalar_and_matrix_min_plus(-2 * pm, u)))),
                        product(range(D + 1), range(D + 1)))))
    G = get_compressed_covers(F)
    H = unite_sets(
        [set(product(*[[(c['i'], c['j']) for c in T['ijminval']] for T in S])) for S in G])
    H = sorted(H, key=lambda a: (
        len(a), number_of_indexes(a, 0) * number_of_indexes(a, 1)))

    for S in H:
        T = apply_simplex(make_simplex_matrix(F, D, S))
        if T is not None:
            res1 = [[i, T[i] + pm] for i in range(D + 1)]
            res2 = [[i, T[D + 1 + i] + pm] for i in range(D + 1)]
            return [res1, res2]

    return None


def test_attack(n, mm, mM, D, pm, pM, ApplyAttack):
    """Runs our attack."""

    A = generate_random_matrix(n, mm, mM)
    B = generate_random_matrix(n, mm, mM)
    p1 = generate_random_polynomial(D, pm, pM)
    p2 = generate_random_polynomial(D, pm, pM)
    u = product_of_two_matrices_min_plus(
        apply_polynomial_min_plus(p1, A), apply_polynomial_min_plus(p2, B))
    q1 = generate_random_polynomial(D, pm, pM)
    q2 = generate_random_polynomial(D, pm, pM)
    v = product_of_two_matrices_min_plus(
        apply_polynomial_min_plus(q1, A), apply_polynomial_min_plus(q2, B))
    KA = product_of_two_matrices_min_plus(apply_polynomial_min_plus(
        p1, A), product_of_two_matrices_min_plus(v, apply_polynomial_min_plus(p2, B)))
    KB = product_of_two_matrices_min_plus(apply_polynomial_min_plus(
        q1, A), product_of_two_matrices_min_plus(u, apply_polynomial_min_plus(q2, B)))
    if KA != KB:
        return False
    attack_result = apply_attack(A, B, u, D, pm)
    if not attack_result:
        return False
    KC = product_of_two_matrices_min_plus(apply_polynomial_min_plus(attack_result[0], A),
                                          product_of_two_matrices_min_plus(v, apply_polynomial_min_plus(attack_result[1], B)))

    if KA != KC:
        return False

    return True


def test_suite(n, mm, mM, D, pm, pM, applyAttack, numberOfTests):
    """Runs a set of tests."""

    st = time.time()
    ok = 0
    fl = 0
    for i in range(1, numberOfTests + 1):
        if test_attack(n, mm, mM, D, pm, pM, applyAttack):
            print("OK")
            ok += 1
        else:
            print("FAIL")
            fl += 1
    et = time.time()
    print(et - st)
    print("OK: ", ok)
    print("FAIL: ", fl)


test_suite(5, -10**2, 10**2, 5, -10**2, 10**2, apply_attack, 100)
