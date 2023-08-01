"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

import random
import math
import time
import functools
import tropical_algebra
from matrix_tools import get_minimum_of_matrix
import scipy.optimize
import itertools
import heapq


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
    return list(
        filter(lambda S: len(S['inds'].difference(unite_sets([T['inds'] for T in filter(lambda T: S != T, Z)]))) != 0,
               Z))


def get_sets_without_elements(Z, N):
    return list(filter(lambda S: len(S['inds']) != 0, [{
        'ijminval': S['ijminval'], 'inds': S['inds'].difference(N)} for S in Z]))


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
    Aub = None
    bub = None
    Aeq = None
    beq = None

    def vecij(i, j):
        v = [0 for _ in range(d1 + d2)]
        v[i] = -1
        v[d1 + j] = -1
        return v

    for i in range(d1):
        for j in range(d2):
            if (i, j) in S:
                if not Aeq:
                    Aeq = []
                Aeq.append(vecij(i, j))
                if not beq:
                    beq = []
                beq.append(M[(i, j)]['val'])
            else:
                if not Aub:
                    Aub = []
                Aub.append(vecij(i, j))
                if not bub:
                    bub = []
                bub.append(M[(i, j)]['val'])
    return c, Aub, bub, Aeq, beq


def get_weighted_sets(S):
    lins = {}
    for T in S:
        for p in T['ijminval']:
            if not p['i'] in lins:
                lins[p['i']] = 0
            lins[p['i']] += 1.0 / len(T['ijminval'])
    cols = {}
    for T in S:
        for p in T['ijminval']:
            if not p['j'] in cols:
                cols[p['j']] = 0
            cols[p['j']] += 1.0 / len(T['ijminval'])
    W = []
    for T in S:
        w = []
        for p in T['ijminval']:
            w.append((p['i'], p['j'], lins[p['i']] * cols[p['j']]))
        W.append([(p[0], p[1])
                  for p in sorted(w, reverse=True, key=lambda x: x[2])])
    return W


def enumerate_product_of_sets(W):
    def enumerate_product_of_sets_(W, s, i):
        if i == len(W) - 1:
            if s < len(W[i]):
                yield [W[i][s]]
            else:
                return
        else:
            for t in range(min(s + 1, len(W[i]))):
                for q in enumerate_product_of_sets_(W, s - t, i + 1):
                    yield [W[i][t]] + q

    l = sum(len(w) - 1 for w in W) + 1
    for s in range(l):
        yield from enumerate_product_of_sets_(W, s, 0)


def enumerate_covers(G):
    for S in sorted(G, key=len):
        W = get_weighted_sets(S)
        yield from enumerate_product_of_sets(W)


def apply_attack(d1, d2, compute_base_element, bounds=(None, None)):
    """Applies our attack. Returns two polynomials p' and q'."""

    def repack(i, j, m):
        return {'ijminval': [{'i': i, 'j': j, 'val': m['val']}], 'inds': m['inds']}

    def check_bound(m):
        if bounds[0] and M[m]['val'] > -2 * bounds[0]:
            return False
        if bounds[1] and M[m]['val'] < -2 * bounds[1]:
            return False
        return True

    M = {(i, j): get_minimum_of_matrix(compute_base_element(i, j))
         for i in range(d1) for j in range(d2)}
    N = [repack(m[0], m[1], M[m]) for m in filter(check_bound, M)]
    G = get_compressed_covers(N)
    for S in enumerate_covers(G):
        c, Aub, bub, Aeq, beq = make_matrices_for_simplex(M, S, d1, d2)

        T = scipy.optimize.linprog(
            c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds)
        if T.success:
            return [[T.x[i] for i in range(d1)], [T.x[d1 + i] for i in range(d2)]]
    return None
