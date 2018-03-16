
from lapsim import simulator
import numpy as np
from scipy.optimize import minimize
from aerofoil_optimization import four_dim_opt, parametric_aerofoil_four_var
from aclab1 import *
from aclabtools import *

# Optimization of Cd and Cl
limits = ((0.0, 1.5), (0.1, 0.5))
x0 = np.array([0.5, 0.5])
opt_out = minimize(simulator, x0, method='L-BFGS-B', args=([0]), bounds=limits, options={'eps': 0.01})
print(opt_out.x[0])

# Wing geometry to achieve optimized Cl
w0 = [0.6, 0.6, 0.6, 0.6]
alpha = 0
target = opt_out.x[0]
w_limit = ((0.2, 2.0), (0.2, 0.8), (0.2, 0.8), (0.2, 2.0))
w_opt_out = four_dim_opt(w0, w_limit, alpha, target)
parametric_aerofoil_four_var(w_opt_out.x)


