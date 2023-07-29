"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random

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
    """Generates k-th basis upper-t-circulant matrix of size n x n."""
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


def mul_basis_upper_t_circulant_matrix_and_matrix(R, t, l, M):
    """Multiplies l-th basis upper-t-circulant matrix and matrix M."""
    n = R.size()

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            k = (n - l + i) % n
            if i - l >= 0:
                b = R.semiring.one()
            else:
                b = t
            result[i][j] = R.semiring.mul(b, M[k][j])

    return result


def mul_matrix_and_basis_upper_t_circulant_matrix(R, t, l, M):
    """Multiplies matrix M and l-th basis upper-t-circulant matrix."""
    n = R.size()

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            k = (j + l) % n
            if j + l < n:
                b = R.semiring.one()
            else:
                b = t
            result[i][j] = R.semiring.mul(M[i][k], b)

    return result
