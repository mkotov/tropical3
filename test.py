import random
import math

def mul(A, B):
    n = len(A)
    C = [[math.inf for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] = min(C[i][j], A[i][k] + B[k][j])

    return C

def pwr(A, n):
    if n == 1:
        return A
    else:
        return mul(A, pwr(A, n - 1))




A = [
[1259, 1239, 1264],
[1264, 1259, 1239],
[1239, 1264, 1259]
]

B = [
[1114, 1108, 1121],
[1134, 1120, 1119],
[1109, 1113, 1125]
]

print(mul(A, B))



