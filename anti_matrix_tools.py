"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import matrix_tools


def generate_random_anti_t_p_circulant_matrix(R, t, p, a_min, a_max):
    """Generates an anti-t-p-circulant matrix of size n x n."""
    n = R.size()
    c_1 = random.randint(a_min, a_max)
    array = [c_1] * n
    for i in range(1, n):
        array[i] = array[i - 1] - p

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            result[i][j] = array[(i - j + n) % n]
            if i + j != n - 1:
                result[i][j] = R.semiring.mul(result[i][j], t)

    return result


def generate_anti_t_p_circulant_matrix(R, array, t):
    """Generates an anti-t-p-circulant matrix of size n x n by array."""
    n = R.size()
    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            result[i][j] = array[(i - j + n) % n]
            if i + j != n - 1:
                result[i][j] = R.semiring.mul(result[i][j], t)

    return result


def generate_basis_anti_t_p_circulant_matrix(R, t, p):
    """Generates basis anti-t-p-circulant matrix of size n x n."""
    n = R.size()

    result = [[R.semiring.zero() for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                result[i][j] = 0
            elif i > j:
                result[i][j] = (i - j) * p
            else:
                result[i][j] = (n - j + i) * p
            if i + j != n - 1:
                result[i][j] = R.semiring.mul(result[i][j], t)

    return result


def mul_basis_anti_t_p_circulant_matrix_and_matrix(R, t, p, M):
    """Multiplies basis anti-t-p-circulant matrix and matrix M."""
    return R.mul(generate_basis_anti_t_p_circulant_matrix(R, t, p), M)


def mul_matrix_and_basis_anti_t_p_circulant_matrix(R, t, p, M):
    """Multiplies matrix M and basis anti-t-p-circulant matrix."""
    return R.mul(M, generate_basis_anti_t_p_circulant_matrix(R, t, p))
