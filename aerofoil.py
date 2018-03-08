
from aclab1 import rational_bezier
import numpy as np
from aclabtools import vortex_panel
import matplotlib.pyplot as plt


def bezier_spline_aerofoil():
    """function that creates aerofoil"""
    # control points
    p = np.array([[1, 0.0], [0.5, 0.08], [0.0, -0.05]])
    q = np.array([[0.0, 0.1], [0.4, 0.2], [1, 0.0]])

    # weights
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, 1, 1])

    # calculate degree
    n = len(p) - 1
    m = len(q) - 1

    # calculate connection points by setting derivative of start and end point to be equal at connection point
    # (equation 19 in notes)
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


def alpha_sweep(points, alpha_min, alpha_max, n_alpha):
    """function that sweeps through different angles and produce lift coefficient graph"""
    alpha = np.linspace(alpha_min, alpha_max, n_alpha)
    cl = []
    cl_lin = []

    for i in alpha:
        result = vortex_panel(points, i, 1)
        cl.append(result[0])
        plt.cla()

    cl_function = np.poly1d(np.polyfit(alpha, cl, 1))
    for i in alpha:
        cl_lin.append(cl_function(i))

    plt.plot(alpha, cl_lin, 'r')  # cl_lin is your 1st order polynomial prediction
    plt.plot(alpha, cl, 'ko')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$C_L$')
    return cl


def parametric_aerofoil(w):
    """Function that creates an aerofoil based on a parameter w"""
    # control points
    p = np.array([[1, 0.0], [0.5, 0.08], [0.0, -0.05]])
    q = np.array([[0.0, 0.1], [0.4, 0.2], [1, 0.0]])

    # weights
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, w, 1])

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


def parameter_sweep(w_array):
    """sweep through a range of w values and return corresponding cl values"""
    cl_list = []
    for i in w_array:
        points = parametric_aerofoil(i)
        result = vortex_panel(points, 0, 1)
        cl_list.append(result[0])
        plt.clf()
    cl = np.array(cl_list)
    return cl


def fixed_lift(w, alpha, target_cl):
    """returns a scalar value of how close the lift of current aerofoil is to the target lift"""
    points = parametric_aerofoil(w)
    result = vortex_panel(points, alpha, 1)
    plt.clf()
    a = (result[0] - target_cl) ** 2
    return a


def fixed_lift_sweep(w_array, alpha, target_cl):
    """returns a scalar value of how close the lift of current aerofoil is to the target lift"""
    a_list = []

    for i in w_array:
        points = parametric_aerofoil(i)
        result = vortex_panel(points, alpha, 1)
        plt.clf()
        a_list.append((result[0] - target_cl) ** 2)

    a_array = np.array(a_list)
    return a_array
