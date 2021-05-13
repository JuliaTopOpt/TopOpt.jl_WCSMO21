using Makie, TopOpt, LinearAlgebra, StatsFuns
using TopOpt.TopOptProblems.Visualization: visualize
Nonconvex.show_residuals[] = true

E = 1.0 # Young’s modulus
v = 0.3 # Poisson’s ratio
f = 1.0 # downward force
rmin = 4.0 # filter radius
xmin = 0.001 # minimum density
problem_size = (160, 40)
x0 = fill(1.0, prod(problem_size))
p = 3.0
stress_threshold = 1.0

problem = PointLoadCantilever(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)
#problem = HalfMBB(Val{:Linear}, problem_size, (1.0, 1.0), E, v, f)

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)

cheqfilter = DensityFilter(solver, rmin = rmin)
stress = TopOpt.MicroVonMisesStress(solver)
comp = TopOpt.Compliance(problem, solver)

function obj(x)
    return sum(cheqfilter(x)) / length(x)
end
function constr1(x)
    return norm(stress(cheqfilter(x)), 30) - stress_threshold
end
function constr2(x)
    return comp(cheqfilter(x)) - 500
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr1)
Nonconvex.add_ineq_constraint!(m, constr2)

options = MMAOptions(maxiter=1000, tol = Tolerance(kkt = 1e-4, f = 1e-4))
TopOpt.setpenalty!(solver, p)
r1 = Nonconvex.optimize(
    m, MMA87(dualoptimizer = ConjugateGradient()),
    x0, options = options,
)

obj(r1.minimizer)
constr1(r1.minimizer)
constr2(r1.minimizer)
maximum(stress(cheqfilter(r1.minimizer)))
topology1 = cheqfilter(r1.minimizer);
fig = visualize(problem; topology = topology1)
Makie.display(fig)
