from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    TOLERANCE = 1e-10

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __getitem__(self,index):
        return self.coordinates[index]

    def __str__(self):
        return 'Value: {}'.format(self.coordinates)

    def __iter__(self):
        return iter(self.coordinates)

    def __eq__(self, str, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coord = [i + j for i, j in zip(self.coordinates, v.coordinates)]
        return Vector(new_coord)

    def minus(self, v):
        new_coord = [i - j for i, j in zip(self.coordinates, v.coordinates)]
        return Vector(new_coord)

    def times_scalar(self, c):
        new_coord = [Decimal(c) * i for i in self.coordinates]
        return Vector(new_coord)

    def magnitude(self):
        return sqrt(sum([i**2 for i in self.coordinates]))


    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0') / Decimal(magnitude))

        except ZeroDivisionError as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:

                raise Exception('cannot compute angle for zero vector')

            else:

                raise e
  
        
    def dot(self, v):
        return sum([i * j for i, j in zip(self.coordinates, v.coordinates)])

    def is_zero(self, tolerance = 1e-10):
        return self.magnitude() < tolerance

    def angle_with(self, v, in_degrees = False):
        try:
            u1 = self.normalized()
            u2 = v.normalized()

            d = round(u1.dot(u2), 4)
            radian_angle = acos(d)
            if in_degrees:
                return radian_angle * 180.0 /pi
            else:
                return radian_angle
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('cannot compute angle for zero vector')
            else:
                raise e

    def is_parallel_to(self, v):
        return(v.is_zero() or v.is_zero() or self.angle_with(v) == 0 or self.angle_with(v) == pi)

    def is_orthogonal_to(self, v, tolerance = 1e-10):
        return abs(self.dot(v)) < tolerance

    def component_parallel_to(self, b):
        ub = b.normalized()
        return ub.times_scalar(self.dot(ub))

    def component_orthogonal_to(self, b):
        projection = self.component_parallel_to(b)
        return self.minus(projection)

    def cross(self, v):
        if len(self.coordinates) > 3:
            return 'undefined'
        x1, y1, z1 = self.coordinates
        x2, y2, z2 = v.coordinates
        return Vector([y1*z2 - y2*z1, x2*z1 - x1*z2, x1*y2 - y1*x2])

    def parallelogram_area(self, v):
        cross_product = self.cross(v)
        return cross_product.magnitude()

    def triangle_area(self, v):
        return 0.5 * self.parallelogram_area(v)

    def get_coordinates(self):
        return self.coordinates

