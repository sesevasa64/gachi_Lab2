import numpy as np


class Vector(np.ndarray):
    def __mul__(self, other: 'Matrix'):
        res = self.dot(other)
        return res / res[-1]


def Vec(*args, **kwargs):
    res = np.array(*args, **kwargs).view(Vector)
    return np.append(res, 1)


class Matrix(np.ndarray):
    def __mul__(self, mat: 'Matrix'):
        return self.dot(mat)

    @staticmethod
    def ones(size) -> 'Matrix':
        return np.eye(size).view(Matrix)

    @staticmethod
    def zeros(size) -> 'Matrix':
        return np.zeros((size, size)).view(Matrix)


def Mat(*args, **kwargs) -> Matrix:
    return np.array(*args, **kwargs).view(Matrix)
