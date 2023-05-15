using HiGHS, JuMP, JuMPModelPlotting

# The paint factory problem
pf = Model(HiGHS.Optimizer)
@variable(pf, x[1:2] >=0)
@constraint(pf, M1_avail, 6*x[1] + 4*x[2] <= 24)
@constraint(pf, M2_avail, x[1] + 2*x[2] <= 6)
@constraint(pf, x_diff, -x[1] + x[2] <= 1)
@constraint(pf, x2_ub, x[2] <= 2)
@objective(pf, Max, 5*x[1] + 4*x[2])

# Has not been solved, no objective function line will be plotted
plt = plot_model(pf)

optimize!(pf)
# Has been solved, objective function line will be plotted
plot_model(pf)

# Make a version of the model with integer variables
pf_int = Model(HiGHS.Optimizer)
@variable(pf_int, x >= 0)
@variable(pf_int, y >= 0, Int)
@constraint(pf_int, M1_avail, 6*x + 4*y <= 24)
@constraint(pf_int, M2_avail, x + 2*y <= 6)
@constraint(pf_int, x_diff, -x + y <= 1)
@constraint(pf_int, x2_ub, y <= 2)
@objective(pf_int, Max, 5*x + 4*y)

set_silent(pf_int)
optimize!(pf_int)
plot_model(pf_int)