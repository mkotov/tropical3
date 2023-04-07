"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import tropical_algebra


def generate_random_matrix(n, mm, mM):
    """Generates a random matrix A in Mat(ZZ, n), a_ij in [mm, mM]."""
    return [[random.randint(mm, mM) for j in range(n)] for i in range(n)]


def generate_random_polynomial(D, pm, pM):
    """Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM]."""
    return [random.randint(pm, pM) for i in range(random.randint(1, D) + 1)]


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

