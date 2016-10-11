from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)



    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]


    def multiply_coefficient_and_row(self, coefficient, row):
        target_plane = self[row]
        p = Plane(target_plane.normal_vector.times_scalar(coefficient), target_plane.constant_term*coefficient)

        self[row] = p


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = (k1 * coefficient) + k2
        new_plane = Plane(normal_vector = new_normal_vector, constant_term = new_constant_term)
        self[row_to_be_added_to] = new_plane

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        def clear_bellow_rows(sys, col, row):
            #row_to_add = Plane(sys.planes[row].normal_vector, sys.planes[row].constant_term)
            for arow in range(row + 1, len(system)):
                if sys[arow][col] != 0:
                    coefficient = (sys[arow][col] * -1) / sys[row][col]
                    sys.add_multiple_times_row_to_row(coefficient, row, arow)


        system = deepcopy(self)
        j = 0
        n = system.dimension
        for i in range(len(system)):
            current_plane = system[i]
            while (j < n):
                c = MyDecimal(current_plane[j])
                if not c.is_zero():
                    clear_bellow_rows(system, j, i)
                    j += 1

                else: # then the current is zero coefficient, find next non zero
                    nonzero_indices = system.indices_of_first_nonzero_terms_in_each_row()
                    for r in range(i + 1, len(nonzero_indices)):
                        # found one, swap, clear bellow and stop the loop
                        if nonzero_indices[r] == j:
                            system.swap_rows(i, r)
                            clear_bellow_rows(system, j, i)
                            j += 1
                            break
                break



        return system

    # def compute_triangular_form(self):
    #     system = deepcopy(self)
    #     nonzero_indices = system.indices_of_first_nonzero_terms_in_each_row()
    #     # swap rows first
    #     for i in range(len(system)):
    #         p = system[i]
    #         if i < system.dimension:
    #             if p.normal_vector.coordinates[i] == 0:
    #                 if i in nonzero_indices:
    #                     next_suitable_nonzero_row = nonzero_indices.index(i)
    #                     if next_suitable_nonzero_row != -1:
    #                         system.swap_rows(i, next_suitable_nonzero_row)
    #                 else:
    #                     raise Exception
    #
    #     #add mutliples of rows
    #     for i in range(1, len(system)):
    #         p = system[i]
    #         if p.normal_vector.coordinates[i-1] != 0:
    #             c = p.normal_vector.coordinates[i-1] * -1
    #             system.add_multiple_times_row_to_row(c, i-1, i)
    #
    #
    #     return system

    def compute_rref(self):
        tf = self.compute_triangular_form()
        m = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        for row in reversed(range(m)):
            current_pivot_index = pivot_indices[row]
            if current_pivot_index ==  -1:
                continue
            row_pivot = tf[row][current_pivot_index]

            if row_pivot != 1:
                # make the pivot 1
                tf.multiply_coefficient_and_row(MyDecimal(1.0)/row_pivot, row)
            adjusted_row_pivot_after_multiplication = tf[row][current_pivot_index]
            proceeding_row = 0
            while proceeding_row < row:
                value_above_pivot = tf[proceeding_row][current_pivot_index]
                if not MyDecimal(value_above_pivot).is_zero():
                    coefficient = (value_above_pivot * MyDecimal(-1.0)) / adjusted_row_pivot_after_multiplication
                    tf.add_multiple_times_row_to_row(coefficient, row, proceeding_row)
                proceeding_row += 1
        return tf

    def raise_exception_if_contradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    def raise_exception_if_too_few_pivots(self):
        non_zero_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if i >= 0 else 0 for i in non_zero_indices])
        if num_pivots < self.dimension:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def solve_linsystem(self):
        try:
            rref = self.compute_rref()

            rref.raise_exception_if_contradictory_equation()
            rref.raise_exception_if_too_few_pivots()

            solution_coordinates = [rref.planes[i].consant_term for i in range(rref.dimension)]
            return Vector(solution_coordinates)
        except Exception as e:
            if str(e) == self.NO_SOLUTIONS_MSG or str(e) == self.INF_SOLUTIONS_MSG:
                return str(e)
            else:
                raise e


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps
p1 = Plane(Vector([5.862, 1.178, -10.366]), '-8.15')
p2 = Plane(Vector([-2.931, -0.589, 5.183]), '-4.075')
s = LinearSystem([p1, p2])
print s.solve_linsystem()



#
# p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
# s = LinearSystem([p1,p2])
# r = s.compute_rref()
# if not (r[0] == p1 and
#         r[1] == Plane(constant_term='1')):
#     print 'test case 2 failed'
#
# p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
# p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
# s = LinearSystem([p1,p2,p3,p4])
# r = s.compute_rref()
# if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
#         r[1] == p2 and
#         r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
#         r[3] == Plane()):
#     print 'test case 3 failed'
#
# p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
# s = LinearSystem([p1,p2,p3])
# r = s.compute_rref()
# if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
#         r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
#         r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
#     print 'test case 4 failed'
# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
# s = LinearSystem([p1, p2])
# t = s.compute_triangular_form()
# if not (t[0] == p1 and
#                 t[1] == p2):
#     print 'test case 1 failed'
#
# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
# s = LinearSystem([p1, p2])
# t = s.compute_triangular_form()
# if not (t[0] == p1 and
#                 t[1] == Plane(constant_term='1')):
#     print 'test case 2 failed'
#
# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
# p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
# s = LinearSystem([p1, p2, p3, p4])
# t = s.compute_triangular_form()
# if not (t[0] == p1 and
#                 t[1] == p2 and
#                 t[2] == Plane(normal_vector=Vector(['0', '0', '-2']), constant_term='2') and
#                 t[3] == Plane()):
#     print 'test case 3 failed'
#
# p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
# s = LinearSystem([p1, p2, p3])
# t = s.compute_triangular_form()
# if not (t[0] == Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2') and
#                 t[1] == Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1') and
#                 t[2] == Plane(normal_vector=Vector(['0', '0', '-9']), constant_term='-2')):
#     print 'test case 4 failed'

#
# p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
# p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
# p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
# p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
#
# s = LinearSystem([p0,p1,p2,p3])
# print s[0]
# s.swap_rows(0,1)
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 1 failed'
#
# s.swap_rows(1,3)
# if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
#     print 'test case 2 failed'
#
# s.swap_rows(3,1)
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 3 failed'
#
# s.multiply_coefficient_and_row(1,0)
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 4 failed'
#
# s.multiply_coefficient_and_row(-1,2)
# if not (s[0] == p1 and
#         s[1] == p0 and
#         s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 5 failed'
#
# s.multiply_coefficient_and_row(10,1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
#         s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 6 failed'
#
# s.add_multiple_times_row_to_row(0,0,1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
#         s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 7 failed'
#
# s.add_multiple_times_row_to_row(1,0,1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
#         s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 8 failed'
#
# s.add_multiple_times_row_to_row(-1,1,0)
# if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
#         s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
#         s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 9 failed'
#
#
#
# print s.indices_of_first_nonzero_terms_in_each_row()
# print '{},{},{},{}'.format(s[0],s[1],s[2],s[3])
# print len(s)
# print s
#
# s[0] = p1
# print s
#
# print MyDecimal('1e-9').is_near_zero()
# print MyDecimal('1e-11').is_near_zero()
