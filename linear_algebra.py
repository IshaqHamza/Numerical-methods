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


class vector(matrix):
    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.dimension = len(coordinates)
        
    def __abs__(self):
        return sqrt(sum([coordinate**2 for coordinate in self.coordinates]))
    
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
    
    

    

    
    