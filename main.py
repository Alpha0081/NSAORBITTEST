import matplotlib.pyplot as plt
import numpy as np
import json

class Vector:
    def __init__(self, coords: np.array) -> None:
        self.__coords: np.array = coords
        self.__normalize()

    def __normalize(self) -> None:
        a: float = self.__get_length()
        self.__coords /= a
    
    def __get_length(self) -> float:
        sum: float = 0
        for coord in self.__coords:
            sum += coord**2
        return sum ** .5

    def get_coords(self) -> np.array:
        return self.__coords

    def get_angle(self, vector) -> float:
        x, y, z = self.__coords
        i, j, k = vector.get_coords()
        return np.arccos(x * i + y * j + z * k) * 180 / np.pi
    
    def draw(self, ax, color=None) -> None:
        x, y, z = self.__coords
        if color is not None:
            ax.plot([0,x], [0,y], [0,z], color = color)
        else:
            ax.plot([0,x], [0,y], [0,z])
        


class Basis:
    def __init__(self, vector1: Vector, vector2: Vector, vector3: Vector) -> None:
        self.__basis = np.array((vector1, vector2, vector3), dtype=Vector)
        print(vector1.get_angle(vector2))
        print(vector2.get_angle(vector3))
        print(vector1.get_angle(vector3))


    def get_matrix(self) -> np.array:
        vector1, vector2, vector3 = self.__basis
        return np.array((vector1.get_coords(), vector2.get_coords(), vector3.get_coords()), dtype=float)

    def draw(self, ax, color=None) -> None:
        for vector in self.__basis:
            vector.draw(ax, color)

class Orbit:
    def __init__(self, eccentricity: float, inclination: float, raan: float, pericenter: float, semi_major: float) -> None:
        self.__eccentricity: float = eccentricity
        self.__inclination: float = inclination
        self.__raan: float = raan
        self.__pericenter: float = pericenter
        self.__normal: Vector = self.__define_normal_vector()
        self.__semi_major: float = semi_major

    def __define_normal_vector(self) -> Vector:
        normal = np.zeros(3)
        normal[0] = np.cos((90 - self.__inclination) * np.pi / 180) * np.cos((self.__raan - 90) * np.pi / 180)
        normal[1] = np.cos((90 - self.__inclination) * np.pi / 180) * np.sin((self.__raan - 90) * np.pi / 180)
        normal[2] = np.sin((90 - self.__inclination) * np.pi / 180)
        return Vector(normal)

    def get_perpendicular(self) -> Vector:
        print(np.cos((self.__raan + 90) * np.pi / 180))
        print(np.sin((self.__raan + 90) * np.pi / 180))
        return Vector(np.array((np.cos((self.__raan + 90) * np.pi / 180), np.sin((self.__raan + 90) * np.pi / 180), 0)))

    def get_focus(self):
        pass

    def get_basis(self, vector1: Vector, vector2: Vector) -> Basis:
        """
        | i   j   k|
        | x1  y1 z1| 
        | x2  y2 z2|  
        """
        x1, y1, z1 = vector1.get_coords()
        x2, y2, z2 = vector2.get_coords()
        i = y1 * z2 - z1 * y2
        j = -(x1 * z2 - z1 * x2)
        k = x1 * y2 - y1 * x2
        return Basis(Vector(np.array((i, j, k))), vector1, vector2)

    def get_normal_vector(self) -> Vector:
        return self.__normal

    def get_line(self) -> Vector:
        return Vector(np.array((np.cos(self.__raan * np.pi / 180), np.sin(self.__raan * np.pi / 180), 0)))

    def get_pericenter_vector(self) -> Vector:
        basis: np.array = self.get_basis(self.__normal, self.get_line()).get_matrix().transpose()
        
        inverse = np.linalg.inv(basis)
        x4 = np.cos(self.__pericenter * np.pi / 180) 
        y4 = np.sin(self.__pericenter * np.pi / 180)
        vector = np.array((x4, y4, 0))
        print(np.dot(inverse, vector).transpose())
        return Vector(np.dot(inverse, vector).transpose())


    def get_trace(self):
        x = np.zeros(100)
        y = np.zeros(100)
        z = np.zeros(100)
        basis = self.get_basis(self.__normal, self.get_pericenter_vector()).get_matrix().transpose()
        e = self.__eccentricity
        a = self.__semi_major
        b = (1 - e**2) ** .5 * a
        x1 = np.linspace(-a, a, 50)
        y1 = (1 - x1**2 / a**2) **.5 * b
        y[:50] = y1
        y[50:] = -y1 
        x[:50] = x1
        x[50:] = x1 
        matrix = np.array((x, y, z))
        inverse2 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        inverse = np.linalg.inv(basis)
        #x, y = np.meshgrid(x, y)
        return np.dot(basis, matrix).transpose()

orbit = Orbit(0, inclination = 0 , raan = 50, pericenter = 90, semi_major = 2)
figure = plt.figure()

ax = figure.add_subplot(1, 1, 1, projection="3d")
#ax.set_axis_off()
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ax.plot([0, 0], [0,0], [0,1], color = "blue")
ax.plot([0,1],[0,0], [0,0], color = "blue")
ax.plot([0, 0], [0, 1],[0, 0], color = "blue")

orbit.get_basis(orbit.get_line(), orbit.get_normal_vector()).draw(ax, "green")

orbit.get_basis(orbit.get_normal_vector(), orbit.get_pericenter_vector()).draw(ax, "red")



t = orbit.get_trace()
for coords in t:
    ax.scatter(coords[0], coords[1], coords[2], color = "blue")

#x,y,z = orbit.get_basis()
#ax.scatter(x, y, z)
plt.show()