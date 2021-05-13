module ContComplianceDemo1

using Makie, TopOpt, LinearAlgebra, StatsFuns
using TopOpt.TopOptProblems.Visualization: visualize

E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
rmin = 4.0 # filter radius
xmin = 0.0001 # minimum density
problem_size = (60, 20)
V = 0.5 # maximum volume fraction
p = 4.0 # penalty
x0 = fill(V, prod(problem_size)) # initial design

problem = HalfMBB(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)
cheqfilter = DensityFilter(solver, rmin = rmin)
comp = TopOpt.Compliance(problem, solver)

function obj(x)
    return comp(cheqfilter(x))
end
function constr(x)
    return sum(cheqfilter(x)) / length(x) - V
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = MMAOptions(
    maxiter=1000, tol = Tolerance(kkt = 1e-4, f = 1e-4),
)
TopOpt.setpenalty!(solver, p)
@time r = Nonconvex.optimize(m, MMA87(),x0, options = options);

obj(r.minimizer)
constr(r.minimizer)
topology = cheqfilter(r.minimizer);
fig = visualize(problem; topology = topology)
Makie.display(fig)

end
