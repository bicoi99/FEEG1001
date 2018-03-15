
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from aclab3 import *
from aclab1 import *


def one_dim_opt(x0, alpha, target_cl):
    """find the minimum of fixed_lift function"""
    opt_out = minimize(fixed_lift, x0, args=(alpha, target_cl))
    return opt_out


def parametric_aerofoil_four_var(w):
    """create an aerofoil based on array of weights w"""
    # control points
    p = np.array([[1, 0.0], [0.5, 0.08], [0.0, -0.05]])
    q = np.array([[0.0, 0.1], [0.4, 0.2], [1, 0.0]])

    # weights
    zp = np.array([1, w[0], w[1], 1])
    zq = np.array([1, w[2], w[3], 1])

    # calculate degree
    n = len(p) - 1
    m = len(q) - 1

    # connection point
    connection_point = (m / (m + n)) * q[0] + (n / (m + n)) * p[n]

    # add connection point
    pp = np.vstack((p, connection_point))
    qq = np.vstack((connection_point, q))

    # calculate 2 curves
    lower = rational_bezier(pp, zp)
    upper = rational_bezier(qq, zq)

    # join 2 curves together
    points = np.vstack([lower[:101, :], upper[1:, :]])

    return points


def fixed_lift_four_var(w, alpha, target_cl):
    """returns a scalar value which represent how close your cl is to target cl"""
    points = parametric_aerofoil_four_var(w)
    result = vortex_panel(points, alpha, 1)
    plt.clf()
    f = (result[0] - target_cl) ** 2
    return f


def four_dim_opt(w0, weight_limits, alpha, target_cl):
    """find the minimum of the fixed_lift_four_var function"""
    opt_out = minimize(fixed_lift_four_var, w0, bounds=weight_limits, args=(alpha, target_cl))
    return opt_out
