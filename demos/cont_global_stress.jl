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
#problem = HalfMBB(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)
TopOpt.setpenalty!(solver, p)

cheqfilter = DensityFilter(solver, rmin = rmin)
stress = TopOpt.MicroVonMisesStress(solver)

# Global stress aggregation

function obj(x)
    return sum(cheqfilter(x)) / length(x)
end
function constr(x)
    return norm(stress(cheqfilter(x)), 5) - 1.0
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = MMAOptions(maxiter=1000, tol = Tolerance(kkt = 1e-3))
x0 = fill(1.0, prod(problem_size))
r = Nonconvex.optimize(m, MMA87(), x0, options = options)

obj(r.minimizer)
constr(r.minimizer)

# Local stress constraints

function constr(x)
    return stress(cheqfilter(x)) .- 1.0
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = AugLagOptions(maxiter=1000, rtol = 1e-4)
r = Nonconvex.optimize(m, AugLag(), x0, options = options)

obj(r.minimizer)
maximum(constr(r.minimizer))
