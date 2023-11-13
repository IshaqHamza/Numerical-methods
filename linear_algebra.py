from numpy.linalg import *
from math import *


class matrix:
    def __init__(self, rows):
        self.rows = rows
        self.dimension = (len(rows), len(rows[0]))

    def __add__(self, other):
        return matrix([[self.rows[i][j] + other.rows[i][j] for j in range(self.dimension[1])] for i in range(self.dimension[0])])
    
    def __sub__(self, other):
        return matrix([[self.rows[i][j] - other.rows[i][j] for j in range(self.dimension[1])] for i in range(self.dimension[0])])
    
    def __mul__(self, other):
        if type(other) == matrix:
            return matrix([[sum([self.rows[i][k]*other.rows[k][j] for k in range(self.dimension[1])]) for j in range(other.dimension[1])] for i in range(self.dimension[0])])
        return matrix([[other*self.rows[i][j] for j in range(self.dimension[1])] for i in range(self.dimension[0])])
    
    def __abs__(self):
        """returns determinant of matrix"""
        return det(self.rows)
    
    def transpose(self):
        return matrix([[self.rows[j][i] for j in range(self.dimension[0])] for i in range(self.dimension[1])])
    
    def inverse(self):
        return matrix(inv(self.rows))
    
    def __call__(self, vector):
        return vector([sum([self.rows[i][j]*vector.coordinates[j] for j in range(self.dimension[1])]) for i in range(self.dimension[0])])
    
    def __str__(self):
        return str(self.rows)
    
    def __repr__(self):
        return str(self.rows)
    
    def __pow__(self, n):
        if n == 0:
            return matrix([[1 if i == j else 0 for j in range(self.dimension[1])] for i in range(self.dimension[0])])
        if n == 1:
            return self
        if n < 0:
            return self.inverse()**(-n)
        return self*(self**(n-1))
    
    def __eq__(self, other):
        return self.rows == other.rows
    
    def __ne__(self, other):
        return self.rows != other.rows
    
    def is_symmetric(self):
        return self == self.transpose()
    
    def is_orthogonal(self):
        return self*self.transpose() == matrix([[1 if i == j else 0 for j in range(self.dimension[0])] for i in range(self.dimension[0])])
    
    def is_idempotent(self):
        return self*self == self
    
    def is_nilpotent(self):
        return self**self.dimension[0] == matrix([[0 for j in range(self.dimension[1])] for i in range(self.dimension[0])])
    
    def is_singular(self):
        return abs(self) == 0
    
    def print(self):
        for row in self.rows:
            print(row)


class vector(matrix):
    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.dimension = len(coordinates)
        
    def __abs__(self):
        return sqrt(sum([coordinate**2 for coordinate in self.coordinates]))
    
    def inf_norm(self):
        return max([abs(coordinate) for coordinate in self.coordinates])
    
    def __add__(self, other):
        return vector([self.coordinates[i] + other.coordinates[i] for i in range(self.dimension)])
    
    def __sub__(self, other):
        return vector([self.coordinates[i] - other.coordinates[i] for i in range(self.dimension)])
    
    def __mul__(self, other):
        if type(other) == vector:
            return sum([self.coordinates[i]*other.coordinates[i] for i in range(self.dimension)])
        return vector([other*self.coordinates[i] for i in range(self.dimension)])
    
    def outer(self, other):
        return matrix([[self.coordinates[i]*other.coordinates[j] for j in range(other.dimension)] for i in range(self.dimension)])
    
    def __str__(self):
        return str(self.coordinates)
    
    def __repr__(self):
        return str(self.coordinates)
    
    def __eq__(self, other):
        return self.coordinates == other.coordinates
    
    def __ne__(self, other):
        return self.coordinates != other.coordinates
    
    def print(self):
        print(self.coordinates)
    
def gauss_elim(A: matrix, b: vector)->vector:
    """returns the solution to the system Ax = b using gaussian elimination"""

    n = A.dimension[0]
    m = A.dimension[1]

    if n != m:
        raise ValueError("Matrix must be square")

    if n != b.dimension:
        raise ValueError("Matrix and vector must have same dimensions")

    Ab = matrix([[A.rows[i][j] for j in range(m)] + [b.coordinates[i]] for i in range(n)])

    for i in range(n-1):
        # find first non-zero entry in row i
        j = i
        while j < n and Ab.rows[j][i] == 0:
            j += 1
        if j == n:
            raise ValueError("Matrix is singular")
        # swap rows i and j
        Ab.rows[i], Ab.rows[j] = Ab.rows[j], Ab.rows[i]
        # subtract multiples of row i from rows below i
        for j in range(i+1, n):
            Ab.rows[j] = [Ab.rows[j][k] - Ab.rows[i][k]*Ab.rows[j][i]/Ab.rows[i][i] for k in range(m+1)]
        
    # back substitution
    x = vector([0 for i in range(n)])
    for i in range(n-1, -1, -1):
        x.coordinates[i] = (Ab.rows[i][-1] - sum([Ab.rows[i][j]*x.coordinates[j] for j in range(i+1, n)]))/Ab.rows[i][i]
    
    return x


def gauss_jack(A: matrix, b: vector, x: vector, N: int = 1000000, tol = 0.001)->vector:
    """returns the solution to the system Ax = b using the Gauss-Jacobi method"""

    n = A.dimension[0]
    m = A.dimension[1]

    if n != m:
        raise ValueError("Matrix must be square")

    if n != b.dimension:
        raise ValueError("Matrix and vector must have same dimensions")

    if n != x.dimension:
        raise ValueError("Matrix and vector must have same dimensions")
    
    for i in range(n):
        if A.rows[i][i] == 0:
            raise ValueError("Matrix is singular")
        
    for _ in range(N):
        x_old = x
        x = vector([(b.coordinates[i] - sum([A.rows[i][j]*x.coordinates[j] for j in range(n) if j != i]))/A.rows[i][i] for i in range(n)])
        x.print()
        if abs(x - x_old) < tol:
            break

    return x

def gauss_seidel(A: matrix, b: vector, x: vector, N: int = 1000)->vector:
    """returns the solution to the system Ax = b using the Gauss-Seidel method"""

    n = A.dimension[0]
    m = A.dimension[1]

    if n != m:
        raise ValueError("Matrix must be square")

    if n != b.dimension:
        raise ValueError("Matrix and vector must have same dimensions")

    if n != x.dimension:
        raise ValueError("Matrix and vector must have same dimensions")

    for i in range(n):
        if A.rows[i][i] == 0:
            raise ValueError("Matrix is singular")

    for _ in range(N):
        for i in range(n):
            x.coordinates[i] = (b.coordinates[i] - sum([A.rows[i][j]*x.coordinates[j] for j in range(n) if j != i]))/A.rows[i][i]

    return x

# print(gauss_elim(matrix([[4, -1, 1], [2, 5, 2], [1, 2, 4]]), vector([8, 3, 11])))
# print(gauss_jack(matrix([[10, -1, 2, 0], [-1, 11, -1, 3], [2, -1, 10, -1], [0, 3, -1, 8]]), vector([6, 25, -11, 15]), vector([0, 0, 0, 0])))
# print(gauss_seidel(matrix([[10, -1, 2, 0], [-1, 11, -1, 3], [2, -1, 10, -1], [0, 3, -1, 8]]), vector([6, 25, -11, 15]), vector([0, 0, 0, 0])))
# print(gauss_jack(matrix([[1, 2, 3], [2, -1, 2], [3, 1, -2]]), vector([5, 1, -1]), vector([0, 0, 0]))) does not converge
# print(gauss_seidel(matrix([[2, 8, 3, 1], [0, 2, -1, 4], [7, -2, 1, 2], [-1, 0, 5, 2]]), vector([-2, 4, 3, 5]), vector([0, 0, 0, 0]))) does not converge

