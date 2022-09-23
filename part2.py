from math import pi, tan, sqrt, acos, sin


def _2d_norm(x, y):
    return sqrt(x * x + y * y)


def cot(x):
    return 1 / tan(x)

eta = pi / 2

# made up device-related numbers

device_screen_hw = 100   # half-width of the screen in pixels
device_screen_ax = 45    # half-aperture of the camera in the width direction, in degrees

# made up rectangle dimensions

real_world_AB_half_sidelength = 30   # the '30' could be cm or whatever other units the real world uses; final answer will be in terms of these units

# made up 2d pixel coordinates of the projection of the rectangle, assuming that (0, 0) is the center of the screen

xA = 10
yA = -12

xB = 3
yB = 33

xM = xA + (xB - xA) * 0.45
yM = yA + (yB - yA) * 0.45

# distance from the camera to the imaginary projection screen, measured in pixels (the same dimensions as the device_screen_hw)

fov = device_screen_hw * cot(device_screen_ax * eta / 90)

FA = sqrt(fov * fov + xA * xA + yA * yA)
FB = sqrt(fov * fov + xB * xB + yB * yB)
FM = sqrt(fov * fov + xM * xM + yM * yM)

AM = _2d_norm(xA - xM, yA - yM)
BM = _2d_norm(xB - xM, yB - yM)

# from law of cosines
#     AM^2 = FA^2 + FM^2 - 2 FA FM cos(alpha)
# ==> 2 FA FM cos(alpha) = FA^2 + FM^2 - AM^2
# ==> cos(alpha) = (FA^2 + FM^2 - AM^2) / (2 FA FM)
# ==> alpha = arccos((FA^2 + FM^2 - AM^2) / (2 FA FM))

alpha = acos((FA**2 + FM**2 - AM**2) / (2 * FA * FM))
beta  = acos((FB**2 + FM**2 - BM**2) / (2 * FB * FM))

print(alpha * 90 / eta)
print(beta * 90 / eta)


def ratio_for_gamma(g):
    g_prime = pi - g
    alpha_prime = pi - g_prime - alpha
    beta_prime = pi - g - beta
    ratio = (sin(alpha) / sin(alpha_prime)) / (sin(beta) / sin(beta_prime))
    print(ratio)
    return ratio


print(ratio_for_gamma(eta - 0.1))

sign = 1 if alpha > beta else -1

step = 0.1
num_reversals = 0

gamma = eta

while num_reversals < 6:
    while ratio_for_gamma(gamma) < 1:
        gamma += sign * step
    while ratio_for_gamma(gamma) > 1:
        gamma -= sign * step
    num_reversals += 1
    step /= 10

alpha_prime = gamma - alpha
real_world_distance_to_AB_midpoint = real_world_AB_half_sidelength * sin(alpha_prime) / sin(alpha)
