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
import scipy.optimize


def compress(G):
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


def get_sets_with_unique_elements(Z):
    """Returns sets that has unique elements."""
    return list(filter(lambda S: len(S['inds'].difference(unite_sets([T['inds'] for T in filter(lambda T: S != T, Z)]))) != 0, Z))


def get_sets_without_elements(Z, N):
    return list(filter(lambda S: len(S['inds']) != 0, [{'ijminval': S['ijminval'], 'inds': S['inds'].difference(N)} for S in Z]))


def get_compressed_covers(F):
    """Returns the compressed set of covers."""
    if len(F) == 0:
        return [[]]
    Z = compress(F)
    M = get_sets_with_unique_elements(Z)
    P = get_sets_without_elements(Z, unite_sets([S['inds'] for S in M]))
    P.sort(key=lambda S: len(S['inds']), reverse=True)

    if len(P) == 0:
        return [M]

    X = [[P[0]] + S for S in get_compressed_covers(
        get_sets_without_elements(P[1:], P[0]['inds']))]
    Y = get_compressed_covers(P[1:])

    return [M + S for S in X + Y]


def make_matrices_for_simplex(M, S, d1, d2):
    """Returns matrices for simplex method."""
    c = [0 for _ in range(d1 + d2)]
    Aub = []
    bub = []
    Aeq = []
    beq = []

    def vecij(i, j):
        v = [0 for _ in range(d1 + d2)]
        v[i] = -1
        v[d1 + j] = -1
        return v
    for i in range(d1):
        for j in range(d2):
            if (i, j) in S:
                Aeq.append(vecij(i, j))
                beq.append(M[(i, j)]['val'])
            else:
                Aub.append(vecij(i, j))
                bub.append(M[(i, j)]['val'])
    return c, Aub, bub, Aeq, beq


def apply_attack(d1, d2, compute_base_element, heuristics_to_sort, bounds=(None, None)):
    """Applies our attack. Returns two polynomials p' and q'."""
    def repack(i, j, m):
        return {'ijminval': [{'i': i, 'j': j, 'val': m['val']}], 'inds': m['inds']}

    def check_bound(m):
        return not bounds[0] or M[m]['val'] <= -2 * bounds[0]

    M = {(i, j): matrix_tools.get_minimum_of_matrix(compute_base_element(i, j))
         for i in range(d1) for j in range(d2)}
    N = [repack(m[0], m[1], M[m]) for m in filter(check_bound, M)]
    G = get_compressed_covers(N)
    H = unite_sets([set(product(*[[(c['i'], c['j']) for c in T['ijminval']] for T in S])) for S in G])
    H = sorted(H, key=heuristics_to_sort, reverse=True)
    for S in H:
        c, Aub, bub, Aeq, beq = make_matrices_for_simplex(M, S, d1, d2)
        T = scipy.optimize.linprog(
            c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds)
        if T.success:
            return [[T.x[i] for i in range(d1)], [T.x[d1 + i] for i in range(d2)]]
    return None
