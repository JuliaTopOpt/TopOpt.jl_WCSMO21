module TrussComplianceDemo2D2

using Makie, TopOpt, LinearAlgebra, StatsFuns
using TopOpt.TrussTopOptProblems.TrussVisualization: visualize

# 2D
ndim = 2
node_points, elements, mats, crosssecs, fixities, load_cases = 
    parse_truss_json("./tim_$(ndim)d.json");
ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
loads = load_cases["0"]
problem = TrussProblem(
    Val{:Linear}, node_points, elements,
    loads, fixities, mats, crosssecs,
);

xmin = 0.0001 # minimum density
x0 = fill(1.0, ncells) # initial design
p = 4.0 # penalty
compliance_threshold = 5.0 # maximum compliance

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)
comp = TopOpt.Compliance(problem, solver)

function obj(x)
    return sum(x) / length(x)
end
function constr(x)
    return comp(x) - compliance_threshold
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = MMAOptions(
    maxiter=1000, tol = Tolerance(kkt = 1e-4, f = 1e-4),
)
TopOpt.setpenalty!(solver, p)
@time r = Nonconvex.optimize(
    m, MMA87(dualoptimizer = ConjugateGradient()),
    x0, options = options,
);

@show obj(r.minimizer)
@show constr(r.minimizer)
fig = visualize(
    problem, solver.u; topology = r.minimizer,
)
Makie.display(fig)

end