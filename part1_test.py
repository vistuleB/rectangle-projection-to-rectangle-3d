from math import sqrt, tan, pi
from numbers import Real


eta = pi / 2


def cot(x):
    return 1 / tan(x)


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
        return "(" + fmt(self.x) + "," + fmt(self.y) + ")"


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
        return "(" + fmt(self.x) + "," + fmt(self.y) + "," + fmt(self.z) + ")"

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


def compute_M(A, B, C, D):
    assert isinstance(A, v2)
    assert isinstance(B, v2)
    assert isinstance(C, v2)
    assert isinstance(D, v2)

    m = m22(B - A, C - D)

    try:
        lambdaMu = m.inverse() * (C - A)
        O = A + (B - A) * lambdaMu.x
        OA = (A - O).norm()
        OB = (B - O).norm()

        # old fundamental equation:
        #     OA / OM = OM / OB
        # ==> OM^2 = OA * OB
        # ==> OM = sqrt(OA * OB)
        # ==> M = O + OM * (A - O) / OA = O + (A - O) * sqrt(OA * OB) / OA = O + (A - O) * sqrt(OB / OA)

        # corrected fundamental equation:
        # OM = 1 / (0.5 * ((1 / OA) + (1 / OB)))

        OM = 1 / (0.5 * ((1 / OA) + (1 / OB)))
        return O + (A - O).normalized() * OM

    except SmallDeterminant:
        return (A + B) / 2


# make up camera data
camera_pos = v3(0, 0, 0)
camera_z = v3(0, 0, 1)
camera_y = v3(0, 1, 0)
camera_x = v3(1, 0, 0)

# (this is actually the identity matrix:)
camera_frame = m33(camera_x, camera_y, camera_z)

# a good camera_frame should be a rotation matrix, i.e., the columns should have length 1
# and be orthogonal; i.e., the following should print out the identity matrix:
# print(camera_frame * camera_frame.transpose())

device_screen_ax = 45   # degrees (half-aperture of camera in x-direction)
device_screen_hw = 1000  # pixels (half-width of phone screen)

# distance from camera to device screen, measured in pixels
screen_distance = device_screen_hw * cot(device_screen_ax * eta / 90)


def project_to_screen(v):
    assert isinstance(v, v3)
    components = camera_frame.transpose() * (v - camera_pos)  # camera_frame.transpose() is actually the identity right now, so right now this is just v - camera_pos
    scaled = (screen_distance / components.z) * components
    return scaled.drop_z()


real_world_AB_half_sidelength = 30  # cm, or whatever real-world units we're using
real_world_BC_half_sidelength = 20  # cm, or whatever real-world units we're using

unit_vector_from_barA_to_barB = v3(-2, 1, 1).normalized()  # choose any numbers; for next line as well!
unit_vector_from_barB_to_barC = v3(1, 1, 1).\
    cross(unit_vector_from_barA_to_barB).normalized()  # using 'cross' is a way to get a vector at 90°, because a.cross(b) is perpendicular to both a and b

barA = v3(100, 100, 1000)
barB = barA + unit_vector_from_barA_to_barB * real_world_AB_half_sidelength * 2
barC = barB + unit_vector_from_barB_to_barC * real_world_BC_half_sidelength * 2
barD = barA + unit_vector_from_barB_to_barC * real_world_BC_half_sidelength * 2

A = project_to_screen(barA)
B = project_to_screen(barB)
C = project_to_screen(barC)
D = project_to_screen(barD)

M = compute_M(A, B, C, D)

print("M:", M)


barM = (barA + barB) / 2

print("true M:", project_to_screen((barA + barB) / 2))
