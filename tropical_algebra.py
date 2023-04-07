"""
(c) I. Buchinskiy, M. Kotov, A. Treier, 2023
"""

from abc import ABC, abstractmethod

INFTY = float('inf')
"""This constant represent +infinity."""

MINFTY = float('-inf')
"""This constant represent -infinity."""


class Semiring(ABC):
    """Semiring."""

    @abstractmethod
    def zero(self):
        """Returns the zero element of the semiring."""
        pass

    @abstractmethod
    def one(self):
        """Returns the identity element of the semiring."""
        pass

    @abstractmethod
    def sum(self, a, b):
        """Returns the sum of two elements of the semiring."""
        pass

    @abstractmethod
    def mul(self, a, b):
        """Returns the product of two elements of the semiring."""
        pass


class R_min_plus(Semiring):
    def zero(self):
        return INFTY

    def one(self):
        return 0

    def sum(self, a, b):
        return min(a, b)

    def mul(self, a, b):
        return a + b


class R_max_plus(Semiring):
    def zero(self):
        return MINFTY

    def one(self):
        return 0

    def sum(self, a, b):
        return max(a, b)

    def mul(self, a, b):
        return a + b


class MatrixSemiring:
    def __init__(self, semiring, n):
        self.semiring = semiring
        self.n = n

    def size(self):
        """Returns the size of matrices."""
        return self.n

    def zero(self):
        """Returns the zero matrix of size n over a semiring."""
        return [[self.semiring.zero() for row in range(self.n)] for col in range(self.n)]

    def one(self):
        """Returns the unit matrix of size n over a semiring."""
        return [[self.semiring.zero() if row != col else self.semiring.one() for col in range(self.n)] for row in range(self.n)]

    def sum(self, A, B):
        """Returns the sum of two matrices over a semiring."""
        return [[self.semiring.sum(A[i][j], B[i][j]) for j in range(self.n)] for i in range(self.n)]

    def mul(self, A, B):
        """Returns the product of two matrices over a semiring."""
        C = self.zero()

        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    C[i][j] = self.semiring.sum(
                        C[i][j], self.semiring.mul(A[i][k], B[k][j]))

        return C

    def mul_by_coef(self, coef, A):
        """Returns the product of an element of a semiring and a matrix over the semiring."""
        return [[self.semiring.mul(A[i][j], coef) for j in range(self.n)] for i in range(self.n)]

    def pwr(self, A, m):
        """Returns a matrix raised to the power m over a semiring."""
        if m == 0:
            return self.one()
        if m % 2 == 0:
            return self.pwr(self.mul(A, A), m // 2)
        else:
            return self.mul(A, self.pwr(A, m - 1))

    def calc_poly(self, p, A):
        """Given a matrix A and a polynomial p over a semiring. Returns p(A)."""
        d = len(p) - 1
        C = self.zero()
        D = self.one()
        for i in range(d + 1):
            C = self.sum(C, self.mul_by_coef(p[i], D))
            if i != d:
                D = self.mul(D, A)

        return C
