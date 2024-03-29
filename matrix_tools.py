"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import tropical_algebra


def generate_random_matrix(R, mm, mM):
    """Generates a random matrix A in Mat(ZZ, n), a_ij in [mm, mM]."""
    n = R.size()
    return [[random.randint(mm, mM) for j in range(n)] for i in range(n)]


def generate_random_polynomial(D, pm, pM):
    """Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM]."""
    return [random.randint(pm, pM) for i in range(D + 1)]


def subtract_matrix_from_matrix(A, B):
    """Returns A - B."""
    return [[x - y for x, y in zip(a, b)] for a, b in zip(A, B)]


def get_minimum_of_matrix(A):
    """Returns minimum of a matrix A and the set of corresponding indexes."""
    m = A[0][0]
    inds = {(0, 0)}

    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] < m:
                m = A[i][j]
                inds = {(i, j)}
            elif A[i][j] == m:
                inds.add((i, j))
    return m, inds


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


def generate_upper_t_circulant_matrix(R, array, t):
    """Generates an upper-t-circulant matrix of size n x n by array."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(array[(i - j + n) % n], t if j > i else one) for j in range(n)] for i in range(n)]


def generate_random_upper_t_circulant_matrix(R, t, a_min, a_max):
    """Generates an upper-t-circulant matrix of size n x n."""
    n = R.size()
    array = [random.randint(a_min, a_max) for i in range(n)]
    return generate_upper_t_circulant_matrix(R, array, t)


def generate_basis_upper_t_circulant_matrix(R, t, k):
    """Generates k-th basis upper-t-circulant matrix of size n x n."""
    n = R.size()
    zero = R.semiring.zero()
    one = R.semiring.one()
    array = [one if k == i else zero for i in range(n)]
    return generate_upper_t_circulant_matrix(R, array, t)


def mul_basis_upper_t_circulant_matrix_and_matrix(R, t, l, M):
    """Multiplies l-th basis upper-t-circulant matrix and matrix M."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(one if i - l >= 0 else t, M[(n - l + i) % n][j]) for j in range(n)] for i in range(n)]


def mul_matrix_and_basis_upper_t_circulant_matrix(R, t, l, M):
    """Multiplies matrix M and l-th basis upper-t-circulant matrix."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(M[i][(j + l) % n], one if j + l < n else t) for j in range(n)] for i in range(n)]


def generate_lower_t_circulant_matrix(R, array, t):
    """Generates a lower-t-circulant matrix of size n x n by array."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(array[(i - j + n) % n], t if j < i else one) for j in range(n)] for i in range(n)]


def generate_random_lower_t_circulant_matrix(R, t, a_min, a_max):
    """Generates a lower-t-circulant matrix of size n x n."""
    n = R.size()
    array = [random.randint(a_min, a_max) for i in range(n)]
    return generate_lower_t_circulant_matrix(R, array, t)


def generate_basis_lower_t_circulant_matrix(R, t, k):
    """Generates k-th basis lower-t-circulant matrix of size n x n."""
    n = R.size()
    zero = R.semiring.zero()
    one = R.semiring.one()
    array = [one if i == k else zero for i in range(n)]
    return generate_lower_t_circulant_matrix(R, array, t)


def mul_basis_lower_t_circulant_matrix_and_matrix(R, t, l, M):
    """Multiplies l-th basis lower-t-circulant matrix and matrix M."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(one if (n - l + i) % n >= i else t, M[(n - l + i) % n][j]) for j in range(n)] for i in range(n)]


def mul_matrix_and_basis_lower_t_circulant_matrix(R, t, l, M):
    """Multiplies matrix M and l-th basis lower-t-circulant matrix."""
    n = R.size()
    one = R.semiring.one()
    return [[R.semiring.mul(M[i][(j + l) % n], one if (j + l) % n <= j else t) for j in range(n)] for i in range(n)]


def generate_anti_t_p_circulant_matrix(R, c_1, t, p):
    """Generates an anti-t-p-circulant matrix of size n x n by its first element."""
    n = R.size()
    one = R.semiring.one()
    array = [c_1 - p * i for i in range(n)]
    return [[R.semiring.mul(array[(i - j + n) % n], t if i + j != n - 1 else one) for j in range(n)] for i in range(n)]


def generate_random_anti_t_p_circulant_matrix(R, t, p, a_min, a_max):
    """Generates an anti-t-p-circulant matrix of size n x n."""
    c_1 = random.randint(a_min, a_max)
    return generate_anti_t_p_circulant_matrix(R, c_1, t, p)


def generate_basis_anti_t_p_circulant_matrix(R, t, p):
    """Generates basis anti-t-p-circulant matrix of size n x n."""
    return generate_anti_t_p_circulant_matrix(R, 0, t, p)
