"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""


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
