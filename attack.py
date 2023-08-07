"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

from matrix_tools import get_minimum_of_matrix
import scipy.optimize
import heapq
import multiprocessing


def compress(G):
    """Compresses a set of pair [a1, L1], [a2, L2], ..., [an, Ln] to the set of pair, where all first components are unique."""

    def find(H, S):
        for i in range(len(H)):
            if S.inds == H[i].inds:
                return i
        return None

    H = []
    for S in G:
        i = find(H, S)
        if i is None:
            H.append(S)
        else:
            H[i].ijs.extend(S.ijs)
    return H


def unite_sets(ss):
    """Unites a collection of sets."""
    return set() if len(ss) == 0 else set.union(*ss)


def get_sets_with_unique_elements(Z):
    """Returns sets that has unique elements."""
    return list(filter(lambda S: len(S.inds.difference(unite_sets([T.inds for T in filter(lambda T: S != T, Z)]))) != 0, Z))


def get_sets_without_elements(Z, N):
    return list(filter(lambda S: len(S.inds) != 0, [Cover(S.inds.difference(N), S.ijs) for S in Z]))


def get_compressed_covers(F):
    """Returns the compressed set of covers."""
    if len(F) == 0:
        return [[]]
    Z = compress(F)
    M = get_sets_with_unique_elements(Z)
    P = get_sets_without_elements(Z, unite_sets([S.inds for S in M]))
    P.sort(key=lambda S: len(S.inds), reverse=True)

    if len(P) == 0:
        return [M]

    X = [[P[0]] + S for S in get_compressed_covers(
        get_sets_without_elements(P[1:], P[0].inds))]
    Y = get_compressed_covers(P[1:])

    return [M + S for S in X + Y]


def compute_preweights(S, R):
    lins = {}
    cols = {}
    for T in S:
        for p in T.ijs:
            i = p[0]
            j = p[1]
            if not i in lins:
                lins[i] = 0
            if not j in cols:
                cols[j] = 0
            if T != R:
                lins[i] += 1.0 / len(T.ijs)
                cols[j] += 1.0 / len(T.ijs)
    return lins, cols


def get_weighted_sets(S, solve_linprog):
    W = []
    mandatory = [(T.ijs[0][0], T.ijs[0][1]) for T in filter(lambda T: len(T.ijs) == 1, S)]

    for T in S:
        lins, cols = compute_preweights(S, T)
        if len(T.ijs) > 1:
            w = []
            for p in T.ijs:
                i = p[0]
                j = p[1]
                if solve_linprog(mandatory + [(i, j)]):
                    w.append((i, j, (lins[i] + 1) * (cols[j] + 1)))
            W.append([(p[0], p[1]) for p in sorted(w, reverse=True, key=lambda x: x[2])])
        else:
            W.append([(T.ijs[0][0], T.ijs[0][1])])

    return W


def enumerate_with_queue(E, chunk_size=10):
    """Enumerates elements using a priority queue with heuristics."""

    def heuristics_to_sort(S):
        def sum_of_cross(S, i, j):
            lins = sum(1 if s[0] == i else 0 for s in S)
            cols = sum(1 if s[1] == j else 0 for s in S)
            return lins * cols

        return -sum(sum_of_cross(S, s[0], s[1])**2 for s in S)

    q = []
    while True:
        k = 0
        for e in E:
            heapq.heappush(q, (heuristics_to_sort(e), e))
            k += 1
            if k == chunk_size:
                break

        if len(q) == 0:
            break

        r = heapq.heappop(q)
        yield r[1]


def enumerate_product_of_sets(W):
    """Enumerates the Cartesian product of sets trying to yield heavier tuples sooner."""
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


class Cover:
    def __init__(self, inds, ijs):
        self.inds = inds
        self.ijs = ijs

def apply_attack(d1, d2, compute_base_element, bounds=(None, None)):
    """Applies our attack. Returns two polynomials p' and q'."""

    M = dict()
    I = []
    for i in range(d1):
       for j in range(d2):
           m, inds = get_minimum_of_matrix(compute_base_element(i, j))
           M[(i, j)] = m
           if (not bounds[0] or m <= -2 * bounds[0]) and (not bounds[1] or m >= -2 * bounds[1]):
               I.append(Cover(inds, [(i, j)]))

    G = get_compressed_covers(I)

    def solve_linprog(S):
        """Solves the linear program corresponding to cover S."""

        def make_matrices_for_linprog(S):
            """Returns matrices for simplex method."""
            c = [0 for _ in range(d1 + d2)]
            Aub = []
            bub = []
            Aeq = []
            beq = []

            for i in range(d1):
                for j in range(d2):
                    v = [-1 if k == i or k == d1 +
                         j else 0 for k in range(d1 + d2)]
                    m = M[(i, j)]
                    if (i, j) in S:
                        Aeq.append(v)
                        beq.append(m)
                    else:
                        Aub.append(v)
                        bub.append(m)

            if Aub == []:
                Aub = None
                bub = None
            if Aeq == []:
                Aeq = None
                beq = None

            return c, Aub, bub, Aeq, beq

        c, Aub, bub, Aeq, beq = make_matrices_for_linprog(S)
        T = scipy.optimize.linprog(
            c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds)
        if T.success:
            return [[T.x[i] for i in range(d1)], [T.x[d1 + i] for i in range(d2)]]

    def enumerate_covers(G):
        """Enumerated covers generated from weighted sets."""

        for S in sorted(G, key=len):
            W = get_weighted_sets(S, solve_linprog)
            yield from enumerate_with_queue(enumerate_product_of_sets(W))

    for S in enumerate_covers(G):
        result = solve_linprog(S)
        if result:
            return result

    return None
