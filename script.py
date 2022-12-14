from math import sqrt, tan, pi, acos, sin
from numbers import Real


eta = pi / 2


def cot(x):
    return 1 / tan(x)


def cot_deg(x):
    return cot(x * eta / 90)


def fmt(double):
    return f"{double:.3f}"


class SmallDeterminant(Exception):
    pass


class v2:
    def __init__(self, x, y):
        assert isinstance(x, Real)
        assert isinstance(y, Real)
        self.x = x
        self.y = y

    def dot(self, other):
        assert isinstance(other, v2)
        return self.x * other.x + self.y * other.y

    def norm(self):
        return sqrt(self.dot(self))

    def __add__(self, other):
        assert isinstance(other, v2)
        return v2(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        assert isinstance(other, v2)
        return v2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        assert isinstance(other, Real)
        return v2(self.x * other, self.y * other)

    def __rmul__(self, other):
        assert isinstance(other, Real)
        return v2(self.x * other, self.y * other)

    def __truediv__(self, other):
        assert isinstance(other, Real)
        return v2(self.x / other, self.y / other)

    def normalized(self):
        n = self.norm()
        return v2(self.x / n, self.y / n)

    def __repr__(self):
        return "(" + fmt(self.x) + ", " + fmt(self.y) + ")"


class v3:
    def __init__(self, x, y, z):
        assert isinstance(x, Real)
        assert isinstance(y, Real)
        assert isinstance(z, Real)
        self.x = x
        self.y = y
        self.z = z

    def dot(self, other):
        assert isinstance(other, v3)
        return self.x * other.x + self.y * other.y + self.z * other.z

    def norm(self):
        return sqrt(self.dot(self))

    def __add__(self, other):
        assert isinstance(other, v3)
        return v3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        assert isinstance(other, v3)
        return v3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        assert isinstance(other, Real)
        return v3(self.x * other, self.y * other, self.z * other)

    def __truediv__(self, other):
        assert isinstance(other, Real)
        return v3(self.x / other, self.y / other, self.z / other)

    def __rmul__(self, other):
        assert isinstance(other, Real)
        return v3(self.x * other, self.y * other, self.z * other)

    def normalized(self):
        n = self.norm()
        return v3(self.x / n, self.y / n, self.z / n)

    def __repr__(self):
        return "(" + fmt(self.x) + ", " + fmt(self.y) + ", " + fmt(self.z) + ")"

    def cross(self, other):
        return v3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x, 
        )

    def drop_x(self):
        return v2(self.y, self.z)

    def drop_y(self):
        return v2(self.x, self.z)

    def drop_z(self):
        return v2(self.x, self.y)


class m22:
    def __init__(self, c1, c2):
        assert isinstance(c1, v2)
        assert isinstance(c2, v2)
        self.a = c1.x
        self.b = c1.y
        self.c = c2.x
        self.d = c2.y

    def det(self):
        return self.a * self.d - self.b * self.c

    def inverse(self):
        D = self.det()
        if abs(D) < 0.001:
            raise SmallDeterminant
        return m22(v2(self.d/D, -self.b/D), v2(-self.c/D, self.a/D))

    def row1(self):
        return v2(self.a, self.c)

    def row2(self):
        return v2(self.b, self.d)

    def col1(self):
        return v2(self.a, self.b)

    def col2(self):
        return v2(self.c, self.d)

    def __mul__(self, other):
        if isinstance(other, m22):
            # self is on the left, other is on the right
            c1 = other.col1()
            c2 = other.col2()
            r1 = self.row1()
            r2 = self.row2()
            return m22(
                v2(c1.dot(r1), c1.dot(r2)),
                v2(c2.dot(r1), c2.dot(r2)),
            )

        elif isinstance(other, v2):
            return v2(self.row1().dot(other), self.row2().dot(other))

        elif isinstance(other, Real):
            return m22(
                self.col1() * other, 
                self.col2() * other,
            )

        else:
            return NotImplemented

    def __neg__(self):
        return self * (-1)

    def __pos__(self):
        return self

    def __repr__(self):
        return fmt(self.a) + " " + fmt(self.c) + "\n" + fmt(self.b) + " " + fmt(self.d)


class m33:
    def __init__(self, c1, c2, c3):
        assert isinstance(c1, v3)
        assert isinstance(c2, v3)
        assert isinstance(c3, v3)

        self.a = c1.x
        self.b = c1.y
        self.c = c1.z

        self.d = c2.x
        self.e = c2.y
        self.f = c2.z

        self.g = c3.x
        self.h = c3.y
        self.i = c3.z

    def a_minor(self):
        return m22(
            self.col2().drop_x(),
            self.col3().drop_x(),
        )

    def b_minor(self):
        return m22(
            self.col2().drop_y(),
            self.col3().drop_y(),
        )

    def c_minor(self):
        return m22(
            self.col2().drop_z(),
            self.col3().drop_z(),
        )

    def d_minor(self):
        return m22(
            self.col1().drop_x(),
            self.col3().drop_x(),
        )

    def e_minor(self):
        return m22(
            self.col1().drop_y(),
            self.col3().drop_y(),
        )

    def f_minor(self):
        return m22(
            self.col1().drop_z(),
            self.col3().drop_z(),
        )

    def g_minor(self):
        return m22(
            self.col1().drop_x(),
            self.col2().drop_x(),
        )

    def h_minor(self):
        return m22(
            self.col1().drop_y(),
            self.col2().drop_y(),
        )

    def i_minor(self):
        return m22(
            self.col1().drop_z(),
            self.col2().drop_z(),
        )

    def det(self):
        return \
            self.a * self.a_minor().det() - \
            self.d * self.d_minor().det() + \
            self.g * self.g_minor().det()

    def row1(self):
        return v3(self.a, self.d, self.g)

    def row2(self):
        return v3(self.b, self.e, self.h)

    def row3(self):
        return v3(self.c, self.f, self.i)

    def col1(self):
        return v3(self.a, self.b, self.c)

    def col2(self):
        return v3(self.d, self.e, self.f)

    def col3(self):
        return v3(self.g, self.h, self.i)

    def transpose(self):
        return m33(self.row1(), self.row2(), self.row3())

    def signed_minors_matrix(self):
        return m33(
            v3(+self.a_minor().det(), -self.b_minor().det(), +self.c_minor().det()),
            v3(-self.d_minor().det(), +self.e_minor().det(), -self.f_minor().det()),
            v3(+self.g_minor().det(), -self.h_minor().det(), +self.i_minor().det()),
        )

    def inverse(self):
        D = self.det()
        if abs(D) < 0.001:
            raise SmallDeterminant
        
        return self.transpose().signed_minors_matrix() / D

    def __truediv__(self, other):
        if isinstance(other, Real):
            return m33(
                self.col1() / other,
                self.col2() / other,
                self.col3() / other,
            )

        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, m33):
            # self is on the left, other is on the right
            r1 = self.row1()
            r2 = self.row2()
            r3 = self.row3()
            c1 = other.col1()
            c2 = other.col2()
            c3 = other.col3()
            return m33(
                v3(c1.dot(r1), c1.dot(r2), c1.dot(r3)),
                v3(c2.dot(r1), c2.dot(r2), c2.dot(r3)),
                v3(c3.dot(r1), c3.dot(r2), c3.dot(r3)),
            )

        elif isinstance(other, v3):
            return v3(
                self.row1().dot(other), 
                self.row2().dot(other),
                self.row3().dot(other),
            )

        elif isinstance(other, Real):
            return m33(
                self.col1() * other,
                self.col2() * other,
                self.col3() * other,
            )

        else:
            return NotImplemented

    def __neg__(self):
        return self * (-1)

    def __pos__(self):
        return self

    def __repr__(self):
        return \
            fmt(self.a) + " " + fmt(self.d) + " " + fmt(self.g) + "\n" + \
            fmt(self.b) + " " + fmt(self.e) + " " + fmt(self.h) + "\n" + \
            fmt(self.c) + " " + fmt(self.f) + " " + fmt(self.i)


def find_the_projection_of_the_A_B_and_C_D_midpoints_from_the_projections_of_A_B_C_D(A, B, C, D):
    # this function assumes that A, B, C, D go in clockwise order around the rectangle
    # (or counter-clockwise would work too, but circular order in any case)
    assert isinstance(A, v2)
    assert isinstance(B, v2)
    assert isinstance(C, v2)
    assert isinstance(D, v2)

    # system of equations
    # A + lambda * (B - A) = D + mu * (C - D)
    # (B - A, D - C) * (lambda, mu) = D - A

    m = m22(B - A, D - C)

    try:
        lambda_and_mu = m.inverse() * (D - A)
        O = A + (B - A) * lambda_and_mu.x
        O_verification = D + (C - D) * lambda_and_mu.y

        # print("")
        # print("O:           ", O)
        # print("compare with:", O_verification)

        OA = (A - O).norm()
        OB = (B - O).norm()
        OC = (C - O).norm()
        OD = (D - O).norm()

        OAB = 1 / (0.5 * ((1 / OA) + (1 / OB)))
        OCD = 1 / (0.5 * ((1 / OC) + (1 / OD)))

        AB_3d_midpoint_projection = O + (A - O).normalized() * OAB
        CD_3d_midpoint_projection = O + (C - O).normalized() * OCD

        return \
            AB_3d_midpoint_projection, \
            CD_3d_midpoint_projection

    except SmallDeterminant:
        return (A + B) / 2, (C + D) / 2


def recover_two_points_from_their_projections_and_the_projection_of_their_3d_midpoint_and_from_their_3d_distance(
        screen_A, 
        screen_B, 
        screen_M, 
        original_distance
    ):
    screen_A_3d = fake_camera_setup.screen_to_world(screen_A, 1)
    screen_B_3d = fake_camera_setup.screen_to_world(screen_B, 1)
    screen_M_3d = fake_camera_setup.screen_to_world(screen_M, 1)

    # let:
    #   - O denote the position of the camera in camera coordinates, i.e., (0, 0, 0)
    #   - alpha be the angle AOM = MOA
    #   - beta be the angle BOM = MOB

    # the following are distances:

    OA = screen_A_3d.norm()
    OB = screen_B_3d.norm()
    OM = screen_M_3d.norm()
    AM = (screen_A_3d - screen_M_3d).norm()
    BM = (screen_B_3d - screen_M_3d).norm()

    # the law of cosines states:
    # AM^2 = OA^2 + OM^2 - 2 OA OM cos(alpha)

    # from which:
    # 2 OA OM cos(alpha) = OA^2 + OM^2 - AM^2
    # cos(alpha) = (OA^2 + OM^2 - AM^2) / (2 OA OM)
    # alpha = arccos((OA^2 + OM^2 - AM^2) / (2 OA OM))

    # by symmetry:
    # beta = arccos((OB^2 + OM^2 - BM^2) / (2 OB OM))

    alpha = acos((OA**2 + OM**2 - AM**2) / (2 * OA * OM))
    beta  = acos((OB**2 + OM**2 - BM**2) / (2 * OB * OM))

    # let:
    #   - gamma be the angle between the ray starting at world_M 
    #     pointing to world_A and the ray starting at world_M and 
    #     pointing opposite to camera
    #   - gamma_prime the complement of gamma: the angle between 
    #     the ray starting at world_M and pointing to world_A and 
    #     the ray starting at world_M and pointing to camera
    # ...except that we don't know the "real gamma" to start with, 
    # we're hunting for it!

    # the following function returns the ratio of distances 
    # world_A <-> world_M and world_B <-> world_M for our imagined 
    # value of gamma (but true, already computed values of alpha
    # and beta); we'll know the right gamma when this ratio hits 1:
    def ratio_for_gamma(gamma):
        gamma_prime = pi - gamma
        alpha_prime = pi - gamma_prime - alpha
        beta_prime = pi - gamma - beta
        ratio = (sin(alpha) / sin(alpha_prime)) / (sin(beta) / sin(beta_prime))
        return ratio

    # hunt for gamma; starting conditions:
    gamma = eta
    last_ratio = ratio_for_gamma(gamma)
    step = 0.1
    num_reversals = 0
    num_iterations = 0

    # binary search:
    while num_reversals < 6:
        num_iterations += 1
        gamma1 = min(gamma + step, pi - 0.05)
        gamma2 = max(gamma - step, 0.05)
        ratio1 = ratio_for_gamma(gamma1)
        ratio2 = ratio_for_gamma(gamma2)

        if (abs(1 - ratio1) < abs(1 - ratio2)):
            gamma = gamma1
            new_ratio = ratio1

        else:
            gamma = gamma2
            new_ratio = ratio2

        if (1 - last_ratio) * (1 - new_ratio) < 0:
            num_reversals += 1
            step /= 10

        last_ratio = new_ratio

    print("")
    print(f"final ratio:  {last_ratio:.7f} ({num_iterations} iterations)")
    print("compare with: 1")

    if abs(1 - last_ratio) > 0.001:
        raise ValueError

    # the true gamma is now ours, as well as the true alpha_prime, beta_prime:
    gamma_prime = pi - gamma
    alpha_prime = pi - gamma_prime - alpha

    # actual distance from the camera to world_M, by law of sines:
    AM_for_real = original_distance / 2
    BM_for_real = original_distance / 2

    OM_for_real = sin(alpha_prime) * AM_for_real / sin(alpha)
    OA_for_real = sin(gamma_prime) * AM_for_real / sin(alpha)
    OB_for_real = sin(gamma) * BM_for_real / sin(beta)

    A_in_camera_coordinates_z_value = OA_for_real / OA  # because OA = screen_A_3d.norm() and because screen_A_3d.z == 1
    B_in_camera_coordinates_z_value = OB_for_real / OB  # because OB = screen_B_3d.norm() and because screen_B_3d.z == 1

    supposed_world_A = fake_camera_setup.screen_to_world(screen_A, A_in_camera_coordinates_z_value)
    supposed_world_B = fake_camera_setup.screen_to_world(screen_B, B_in_camera_coordinates_z_value)

    return supposed_world_A, supposed_world_B


def find_four_corners_from_projections_and_AB_CD_sidelengths(
    screen_A,
    screen_B,
    screen_C,
    screen_D,
    AB_real_world_sidelength,
    CD_real_world_sidelength, # as long as AB, CD are parallel in the real world, this function works (more general than rectangles)
):
    # BE CAREFUL THIS FUNCTION ASSUMES THAT A, B, C, D APPEAR EITHER
    # IN PURELY CLOCKWISE OR PURELY COUNTERCLOCKWISE ORDER
    screen_AB_mid, screen_CD_mid = find_the_projection_of_the_A_B_and_C_D_midpoints_from_the_projections_of_A_B_C_D(
        screen_A, 
        screen_B, 
        screen_C, 
        screen_D
    )

    world_A, world_B = recover_two_points_from_their_projections_and_the_projection_of_their_3d_midpoint_and_from_their_3d_distance(
        screen_A,
        screen_B,
        screen_AB_mid,
        AB_real_world_sidelength,
    )

    world_C, world_D = recover_two_points_from_their_projections_and_the_projection_of_their_3d_midpoint_and_from_their_3d_distance(
        screen_C,
        screen_D,
        screen_CD_mid,
        CD_real_world_sidelength,
    )

    return world_A, world_B, world_C, world_D


def cook_up_example(
    A_position,
    AB_direction,
    other_thing_that_will_turned_into_BC_direction,
    AB_length,
    AC_length
):
    unit_vector_from_world_A_to_world_B = AB_direction.normalized()
    unit_vector_from_world_B_to_world_C = \
        other_thing_that_will_turned_into_BC_direction.cross(unit_vector_from_world_A_to_world_B).normalized()

    world_A = A_position
    world_B = world_A + unit_vector_from_world_A_to_world_B * AB_length
    world_C = world_B + unit_vector_from_world_B_to_world_C * AC_length
    world_D = world_A + unit_vector_from_world_B_to_world_C * AC_length
    
    screen_A = fake_camera_setup.world_to_screen(world_A)
    screen_B = fake_camera_setup.world_to_screen(world_B)
    screen_C = fake_camera_setup.world_to_screen(world_C)
    screen_D = fake_camera_setup.world_to_screen(world_D)

    return (screen_A, screen_B, screen_C, screen_D), (world_A, world_B, world_C, world_D)


class fake_camera_setup:
    camera_pos = v3(0, 0, 0)
    camera_x = v3(1, 0, 0)                                              # you must take care that...
    camera_y = v3(0, 1, 0)                                              # ...these three vectors...
    camera_z = v3(0, 0, 1)                                              # ...are orthonormal, or else incorrect results will be obtained
    camera_frame = m33(camera_x, camera_y, camera_z)                    # camera_x, camera_y, camera_z are columns of the matrix
    device_screen_ax = 45                                               # degrees (half-aperture of camera in x-direction)
    device_screen_hw = 1000                                             # pixels (half-width of phone screen)
    screen_distance = device_screen_hw * cot_deg(device_screen_ax)      # pixels

    @classmethod
    def world_to_screen(cls, in_world):
        assert isinstance(in_world, v3)
        in_camera = cls.camera_frame.transpose() * (in_world - cls.camera_pos)  # (like applying Unity's worldToCameraMatrix)
        return (in_camera * cls.screen_distance / in_camera.z).drop_z()         # (like applying Unity's projectionMatrix)

    @classmethod
    def screen_to_world(cls, on_screen, z):
        assert isinstance(on_screen, v2)
        assert isinstance(z, Real)
        in_camera = v3(on_screen.x, on_screen.y, cls.screen_distance) * (z / cls.screen_distance)  # (like un-applying Unity's projectionMatrix)
        return cls.camera_pos + cls.camera_frame * in_camera                                       # (like un-applying Unity's worldToCameraMatrix)


QR_code_size = 2


screen_ABCD, actual_ABCD = cook_up_example(
    v3(100, 100, 1000),
    v3(-2, 1, 1),
    v3(1, 1, 1),
    QR_code_size,
    QR_code_size,
)

reconstructed_ABCD = find_four_corners_from_projections_and_AB_CD_sidelengths(
    screen_ABCD[0],
    screen_ABCD[1],
    screen_ABCD[2],
    screen_ABCD[3],
    QR_code_size,
    QR_code_size
)


print("")
for i in range(4):
    print(f"actual_ABCD[{i}]:", actual_ABCD[i])

print("")
print("compare with:")
for i in range(4):
    print(f"reconstructed_ABCD[{i}]:", reconstructed_ABCD[i])


