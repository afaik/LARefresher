import poc_simpletest
from vector import Vector
from decimal import Decimal
from line import Line
from plane import Plane
def add_vectors(v1, v2):
    return Vector([i + j for i,j in zip(v1.get_coordinates(), v2.get_coordinates())])
def sub_vectors(v1, v2):
    return Vector([i - j for i,j in zip(v1.get_coordinates(), v2.get_coordinates())])
def scalar_mult(v, s):
    return Vector([i * s for i in v.get_coordinates()])

v1 = Vector([1, -4])
v2 = Vector([2, 5])
v_test = Vector([3,4])
tester = poc_simpletest.TestSuite()
tester.run_test(add_vectors(v1, v2).get_coordinates(), (3, 1), 'addition test')
tester.run_test(sub_vectors(v1, v2).get_coordinates(), (-1, -9), 'subtraction test')
tester.run_test(scalar_mult(v1, 2).get_coordinates(), (2, -8), 'scalar multiplication test')
tester.run_test(v_test.magnitude(), 5.0 , 'magnitude')
#tester.run_test(v_test.direction().get_coordinates(), (Decimal(3)/5.0, Decimal(4)/5.0), 'calc direction')

##v3 = vector.Vector([7.119, 8.215])
##v4 = vector.Vector([-8.223, 0.878])
##v5 = vector.Vector([8.218, -9.341])
##v6 = vector.Vector([-1.129, 2.111])
##v7 = vector.Vector([1.671, -1.012, -0.318])
##s = 7.41
##print add_vectors(v5, v6)
##print sub_vectors(v3, v4)
##print scalar_mult(v7, s)
##
##print v5.plus(v6)
##print v3.minus(v4)
##print v7.times_scalar(s)

##tester.report_results()



# print vector.Vector([7.887, 4.138]).dot(vector.Vector([-8.802, 6.776]))
# print vector.Vector([-5.955, -4.904, -1.874]).dot(vector.Vector([-4.496, -8.755, 7.103]))
# print vector.Vector([3.183, -7.627]).angle_with(vector.Vector([-2.668, 5.319]))
# print vector.Vector([7.35, 0.221, 5.188]).angle_with(vector.Vector([2.751, 8.259, 3.985]), True)
#
#
# print Vector([-7.579, -7.88]).is_parallel_to(Vector([22.737, 23.64]))
# print Vector([-7.579, -7.88]).is_orthogonal_to(Vector([22.737, 23.64]))
# v = [-2.029, 9.97, 4.172]
# w = [-9.231, -6.639, -7.245]
# print Vector(v).is_parallel_to(Vector(w))
# print Vector(w).is_orthogonal_to(Vector(w))
l1 = Line(Vector(['1.182', '5.562']), '6.744')
l2 = Line(Vector(['1.773', '8.343']), '9.525')
p1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46 )
p2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65 )
print p1.is_parallel_to(p2)
print 't1', p1 == p2


