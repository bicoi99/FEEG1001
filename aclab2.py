import numpy as np
import matplotlib.pyplot as pl
from aclab1 import rational_bezier
from aclabtools import vortex_panel
        
def bezier_spline_aerofoil():
    
    # define the control points:
    #p = np.array(...)
    #q = np.array(...)
    ### BEGIN SOLUTION
    p = np.array([[1, 0.0], [0.5, 0.08], [0.0, -0.05]])
    q = np.array([[0.0, 0.1], [0.4, 0.2], [1, 0.0]])
    ### END SOLUTION
    
    # weights
    zp = np.array([1, 1, 1, 1])
    zq = np.array([1, 1, 1, 1])
    
    # calculate degree
    #n = ...
    #m = ...
    ### BEGIN SOLUTION
    n = np.float(np.shape(p)[0])
    m = np.float(np.shape(p)[0])
    ### END SOLUTION
    
    # calculate connection point
    # q_start = p_end = ...
    ### BEGIN SOLUTION
    q_start = p_end = (n/(m+n))*p[-1,:]+(m/(m+n))*q[0,:]
    ### END SOLUTION
    
    # and add to control points
    pp = np.vstack([p, p_end])
    qq = np.vstack([q_start, q])
    
    # calculate two curves 
    # lower = ...
    # upper = ...
    ### BEGIN SOLUTION
    lower = rational_bezier(pp, zp)
    upper = rational_bezier(qq, zq)
    ### END SOLUTION
    
    # and join together (removing repeat point at leading edge):
    ### BEGIN SOLUTION
    points = np.vstack([lower[0:100], upper])
    ### END SOLUTION

    return points

def alpha_sweep(points, alpha_min, alpha_max, n_alpha):
    
    alpha = np.linspace(alpha_min, alpha_max, n_alpha)
    
    ### BEGIN SOLUTION
    # initialise arrays
    cl = np.zeros([n_alpha])
    for i in range(0, n_alpha):
        [cl[i], cp, xc] = vortex_panel(points, alpha[i], 0)
    z = np.polyfit(alpha, cl , 1)
    p = np.poly1d(z)    
    cl_lin = [p(x) for x in alpha]
    ### END SOLUTION
    
    # plot results
    pl.plot(alpha, cl_lin, 'r') #cl_lin is your 1st order polynomial prediction
    pl.plot(alpha, cl, 'ko')
    pl.xlabel(r'$\alpha$')
    pl.ylabel(r'$C_L$')
    return cl
