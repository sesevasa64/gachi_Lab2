from typing import *
from linalg import *
import matplotlib.pyplot as plt
from copy import deepcopy


# Массив точек
class Face:
    def __init__(self, vertex: List[Vector] = None):
        if vertex is None:
            self.vertex = []
        else:
            self.vertex = vertex

    def __getitem__(self, idx) -> Vector:
        return self.vertex[idx]

    def __setitem__(self, idx, vec: Vector):
        self.vertex[idx] = vec

    def __len__(self):
        return len(self.vertex)

    def __mul__(self, other: Matrix):
        res = deepcopy(self)
        for i in range(len(res)):
            res[i] *= other
        return res

    def __imul__(self, other: Matrix):
        for i in range(len(self)):
            self[i] = self[i] * other
        return self

    def __iadd__(self, vec: Vector):
        self.vertex.append(vec)
        return self

    def plot(self, ax, color):
        x = [v[0] for v in self.vertex]
        y = [v[1] for v in self.vertex]
        z = [v[2] for v in self.vertex]
        ax.plot_wireframe(x, y, np.array([z, z]), color=color)


# Массив граней
class Square:
    def __init__(self, square: Union[List[Face], 'Square'] = None):
        if square is None:
            self.square = []
        else:
            self.square = square

    def __getitem__(self, idx) -> Face:
        return self.square[idx]

    def __setitem__(self, idx, vec: Face):
        self.square[idx] = vec

    def __len__(self):
        return len(self.square)

    def __mul__(self, other: Matrix):
        res = deepcopy(self)
        for i in range(len(res)):
            res[i] *= other
        return res

    def __imul__(self, other: Matrix):
        for i in range(len(self)):
            self[i] *= other
        return self

    def __iadd__(self, other: Union[Face, 'Square']):
        self.square.append(other)
        return self

    def delete(self, idx):
        del self.square[idx]

    def plot(self, ax, color):
        for face in self.square:
            face.plot(ax, color)


# Массив плоскостей
class Figure:
    def __init__(self, name):
        self.planes = []
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        with open(name) as file:
            i = 0
            faces = []
            self.planes.append(Square())
            for line in file.read().splitlines()[1:]:
                if line == '-' or line[0] == '-':
                    f = Face(faces)
                    self.planes[i] += f
                    faces = []
                    if line != '-':
                        self.planes.append(Square())
                        i += 1
                else:
                    arr = np.fromstring(line, sep='\t')
                    res = np.append(arr, 1).view(Vector)  # Для однородных
                    faces.append(res)

    def __getitem__(self, idx) -> Square:
        return self.planes[idx]

    def __setitem__(self, idx, square: Square):
        self.planes[idx] = square

    def __len__(self):
        return len(self.planes)

    def __mul__(self, other: Matrix):
        res = deepcopy(self)
        for i in range(len(res)):
            res[i] *= other
        return res

    def __imul__(self, other: Matrix):
        for i in range(len(self)):
            self[i] *= other
        return self

    def __iadd__(self, other: Union[Face, Square]):
        self.planes.append(other)
        return self

    def save(self, name):
        with open(name, 'w') as file:
            for j in range(len(self)):
                plane = self.planes[j]
                file.write(f'-{self.colors[j]}-\n')
                for i in range(len(plane)):
                    face = plane[i]
                    for v in face:
                        file.write(f'{v[0]:5}\t{v[1]:5}\t{v[2]:5}\n')
                    if i != len(plane) - 1:
                        file.write('-\n')
            file.write('-\n')

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        for i in range(len(self.planes)):
            self.planes[i].plot(ax, self.colors[i])
        plt.show()
