using TopOpt, LinearAlgebra

E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
rmin = 3.0 # filter radius
V = 0.5 # volume fraction
xmin = 0.001 # minimum density
problem_size = (60, 20)
p = 3.0

problem = PointLoadCantilever(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)
#problem = HalfMBB(Val{:Linear}, (60, 20), (1.0, 1.0), E, v, f)

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)
TopOpt.setpenalty!(solver, p)

cheqfilter = DensityFilter(solver, rmin = rmin)
stress = TopOpt.MicroVonMisesStress(solver)

function obj(x)
    return sum(cheqfilter(x)) / length(x)
end
function constr(x)
    return norm(stress(cheqfilter(x)), 5) - 3.0
end

options = MMAOptions(maxiter=1000, tol = Tolerance(kkt = 1e-3))
x0 = fill(0.5, prod(problem_size))

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)
r = Nonconvex.optimize(m, MMA87(), x0, options = options)
obj(r.minimizer)
constr(r.minimizer)
