from math import sqrt

def _2d_norm(x, y):
    return sqrt(x * x + y * y)


# some made-up numbers:
(xA, yA) = (2, 3)
(xB, yB) = (4, 5)
(xC, yC) = (5, 6)
(xD, yD) = (7, 7)


a = xB - xA
b = yB - yA
c = xC - xD
d = yC - yD


det = a * d - b * c


if (abs(det) < 0.001):
    xM = (xA + xB) / 2
    yM = (yA + yB) / 2

else:
    def _2_2_matrix_inverse(A, B, C, D):
        DET = A * D - B * C
        return (D / DET, -B / DET, -C / DET, A / DET)

    def _2_2_matrix_times_vector(A, B, C, D, X, Y):
        return (A * X + C * Y, B * X + D * Y)

    (a_, b_, c_, d_) = _2_2_matrix_inverse(a, b, c, d)

    e = xC - xA
    f = yC - yA

    # note: I cannot write 'lambda' because 'lambda' is a keyword in python; otherwise I would write 'lambda'

    (lamda, mu) = _2_2_matrix_times_vector(a_, b_, c_, d_, e, f)
    
    xO = xA + lamda * (xB - xA)
    yO = yA + lamda * (yB - yA)

    xO_verification = xC + mu * (xD - xC)
    yO_verification = yC + mu * (yD - yC)

    print(xO, xO_verification)
    print(yO, yO_verification)

    # fundamental equation:
    #     OA / OM = OM / OB
    # ==> OM^2 = OA * OB
    # ==> OM = sqrt(OA * OB)
    # ==> M = O + OM * (A - O) / OA = O + (A - O) * sqrt(OA * OB) / OA = O + (A - O) * sqrt(OB / OA)

    OA = _2d_norm(xA - xO, yA - yO)
    OB = _2d_norm(xB - xO, yB - yO)

    multiplier = sqrt(OB / OA)

    xM = xO + multiplier * (xA - xO)
    yM = yO + multiplier * (yA - yO)


print(xM, yM)