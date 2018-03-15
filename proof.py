
from aclab1 import *
from aclabtools import *
import matplotlib.pyplot as plt

target_cl = (0.7746 / (0.9 * 0.4)) / 0.8
u2 = 0.5


# function to generate airfoil
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
    plt.cla()

    # join 2 curves together
    points = np.vstack([lower[:101, :], upper[1:, :]])

    return points


points = parametric_aerofoil2(u2)
alpha_plot = np.linspace(-5, 45, 20)
cl_plot = []
for i in alpha_plot:
    cl_plot.append(vortex_panel(points, i, 0)[0])

plt.plot(alpha_plot, cl_plot, 'k-', label='Cl plot for this airfoil geometry')  # Cl against alpha plot
plt.plot(np.linspace(-5, 45, 50), target_cl * np.ones(50), label='target Cl line')  # target CL line
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_L$')
plt.legend()
plt.show()
