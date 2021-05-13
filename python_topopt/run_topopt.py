# pip install cvxpy nlopt topopt
# https://github.com/zfergus/topopt

import time
import numpy
from topopt.boundary_conditions import MBBBeamBoundaryConditions, CantileverBoundaryConditions
from topopt.problems import ComplianceProblem
from topopt.solvers import TopOptSolver
from topopt.filters import DensityBasedFilter, SensitivityBasedFilter
from topopt.guis import GUI
from termcolor import cprint

nelx, nely = 720, 240 # Number of elements in the x and y
volfrac = 0.3  # Volume fraction for constraints
penal = 3.0  # Penalty for SIMP
rmin = 2.0  # Filter radius

start_time = time.time()
# Initial solution
x = volfrac * numpy.ones(nely * nelx, dtype=float)

# Boundary conditions defining the loads and fixed points
# bc = MBBBeamBoundaryConditions(nelx, nely)
bc = CantileverBoundaryConditions(nelx, nely)

# Problem to optimize given objective and constraints
problem = ComplianceProblem(bc, penal)
topopt_filter = DensityBasedFilter(nelx, nely, rmin)
# topopt_filter = SensitivityBasedFilter(nelx, nely, rmin)

solver = TopOptSolver(problem, volfrac, topopt_filter, maxeval=1000)
solver.xtol_abs = 1e-3
solver.ftol_abs = 1e-3
solver.ftol_rel = 1e-3

x_opt = solver.optimize(x)

cprint('Total time (problem setup + opt): {}'.format(time.time() - start_time), 'cyan')

gui = GUI(problem, "Topology Optimization Example")
gui.update(x_opt)
input('Press enter to exit')
