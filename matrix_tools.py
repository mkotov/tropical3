"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import tropical_algebra


def generate_random_matrix(R, mm, mM):
    """Generates a random matrix A in Mat(ZZ, n), a_ij in [mm, mM]."""
    n = R.size()
    return [[random.randint(mm, mM) for j in range(n)] for i in range(n)]


def generate_random_upper_t_circulant_matrix(R, t, a_min, a_max):
    """Generates an upper-t-circulant matrix of size n x n."""
    n = R.size()
    array = [random.randint(a_min, a_max) for i in range(n)]

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            result[i][j] = array[(i - j + n) % n]
            if j > i:
                result[i][j] = R.semiring.mul(result[i][j], t)

    return result


def generate_upper_t_circulant_matrix(R, array, t):
    """Generates an upper-t-circulant matrix of size n x n by array."""
    n = R.size()
    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            result[i][j] = array[(i - j + n) % n]
            if j > i:
                result[i][j] = R.semiring.mul(result[i][j], t)

    return result


def generate_basis_upper_t_circulant_matrix(R, t, k):
    """Generates an upper-t-circulant matrix of size n x n."""
    n = R.size()

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            if (i - j + n) % n == k:
                if j <= i:
                    result[i][j] = R.semiring.one()
                else:
                    result[i][j] = t

    return result


def generate_random_polynomial(D, pm, pM):
    """Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM]."""
    return [random.randint(pm, pM) for i in range(D + 1)]


def minus_matrix_from_matrix(A, B):
    """Returns A - B. Some elements of the matrices can be infinity."""
    return [[x - y for x, y in zip(row_a, row_b)] for row_a, row_b in zip(A, B)]


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


def get_maximum_of_matrix(A):
    """Returns maximim of a matrix A and the set of corresponding indexes."""
    best = {'val': A[0][0], 'inds': set()}
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] > best['val']:
                best = {'val': A[i][j], 'inds': {(i, j)}}
            elif A[i][j] == best['val']:
                best['inds'].add((i, j))
    return best


def number_of_indexes(S, i):
    """Returns the number of i-th indexes."""
    return len(set([c[i] for c in S]))


def spectrum(S, j):
    """Returns the distribution of indexes."""
    a = [c[j] for c in S]
    return sorted([a.count(i) for i in set(a)], reverse=True)


def is_matrix_minus_matrix_const(R, A, B):
    """
    If A - B is a const, then returns this const. Returns None otherwise.
    """
    n = len(A)
    ra = None
    rb = None
    ok = False
    # Find a pair of non-zero finite elements of matrices, to compute a const.
    for i in range(n):
        for j in range(n):
            if A[i][j] != R.semiring.zero() and B[i][j] != R.semiring.zero():
                ra = A[i][j]
                rb = B[i][j]
                ok = True
                break
        if ok:
            break

    if not ra or not rb:
        return None

    for i in range(n):
        for j in range(n):
            if A[i][j] == R.semiring.zero() and B[i][j] == R.semiring.zero():
                continue
            if A[i][j] == R.semiring.zero() or B[i][j] == R.semiring.zero():
                return None

            if A[i][j] - B[i][j] != ra - rb:
                return None

    return ra - rb


def is_matrix_repeated(R, As):
    """
    Returns True iff the last matrix is const * As[i] for some i < len(As) - 1.
    """
    for i in range(len(As) - 1):
        if is_matrix_minus_matrix_const(R, As[-1], As[i]):
            return True

    return False


def get_first_repeated(R, A, bound):
    """
    For a matrix A, returns n s. t. A^n = A^m, m < n.
    """
    As = [A]
    for i in range(2, bound):
        As.append(R.mul(As[-1], A))
        if is_matrix_repeated(R, As):
            return i

    return None
