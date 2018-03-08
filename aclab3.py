import numpy as np
import matplotlib.pyplot as pl
from aclab1 import rational_bezier
from aclabtools import vortex_panel

def parametric_aerofoil(w):

    p = np.array([[1, 0],[0.5, 0.08],[0.0, -0.05]])
    q = np.array([[0, 0.1],[0.4, 0.2],[1, 0.0]])
    # weights
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, w, 1])
    
    # calculate degree
    n = np.float(np.shape(p)[0])
    m = np.float(np.shape(p)[0])
    
    # calculate connection point
    q_start = p_end = (n / (m + n)) * p[-1, :] + (m / (m + n)) * q[0, :]
    
    # and add to control points
    pp = np.vstack([p, p_end])
    qq = np.vstack([q_start, q])
    
    # calculate two curves
    lower = rational_bezier(pp, zp)
    upper = rational_bezier(qq, zq)
    
    points=np.vstack([lower[0:100],upper])

    return points
    
def fixed_lift(w, alpha, target_cl):
    points = parametric_aerofoil(w)
    [cl, cp, xc] = vortex_panel(points, alpha, 0)
    a = (cl - target_cl)**2
    return a