
from aclab1 import *
from aclabtools import *
from scipy.optimize import minimize


def parametric_aerofoil2(u2):
    """Create an aerofoil with input u2 as the upper spline 2nd control point"""
    # control points
    p = np.array([[1, 0.0], [0.5, 0.08], [0.0, -0.05]])
    q = np.array([[0.0, 0.1], [u2, 0.2], [1, 0.0]])

    # weights
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, 1, 1])

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


def find_alpha(target_cl, points, guess):
    """Find an angle of attack, alpha, of a given airfoil which has the target cl"""
    def find_alpha_fixed_lift(alpha_f, points_f, target_cl_f):
        result = vortex_panel(points_f, alpha_f, 0)
        a = (result[0] - target_cl_f) ** 2
        return a
    opt_out = minimize(find_alpha_fixed_lift, guess, args=(points, target_cl))
    alpha = float(opt_out.x)
    cl = vortex_panel(points, alpha, 0)[0]
    return alpha, cl


def drag(u2, target_cl, U, c, plotting):
    """Return cl/cd and cd"""
    nu = 1.45e-5
    Re = U * c / nu
    points = parametric_aerofoil2(u2)
    alpha = find_alpha(target_cl, points, 5)[0]
    result = viscous_solver(points, alpha, Re, nu, plotting)
    cl, cd = result[0], result[1]
    return cl/cd, cd
