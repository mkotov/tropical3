"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import math
import time
import functools
from itertools import product
import tropical_algebra
import matrix_tools
import simplex


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
    """Unites a collection of sets."""
    if len(ss) == 0:
        return set()
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


def apply_attack(R, A, B, u, D, pm):
    """Applies our attack. Returns two polynomials p' and q'."""
    def repack(ij, m):
        return {'ijminval': [{'i': ij[0], 'j': ij[1], 'val': m['val']}], 'inds': m['inds']}

    F = list(filter(lambda S: S['ijminval'][0]['val'] <= 0,
                    map(lambda ij: repack(ij, matrix_tools.get_minimum_of_matrix(
                        matrix_tools.minus_matrix_from_matrix(R.mul(R.pwr(A, ij[0]), R.pwr(B, ij[1])), R.mul_by_coef(-2 * pm, u)))),
                        product(range(D + 1), range(D + 1)))))
    G = get_compressed_covers(F)
    H = unite_sets(
        [set(product(*[[(c['i'], c['j']) for c in T['ijminval']] for T in S])) for S in G])
    H = sorted(H, key=lambda a: (
        len(a), number_of_indexes(a, 0) * number_of_indexes(a, 1)))

    for S in H:
        T = simplex.apply_simplex(make_simplex_matrix(F, D, S))
        if T is not None:
            res1 = [T[i] + pm for i in range(D + 1)]
            res2 = [T[D + 1 + i] + pm for i in range(D + 1)]
            return [res1, res2]

    return None
