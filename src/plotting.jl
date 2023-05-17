using JuMP, Polyhedra, Plots, SparseArrays


function plot_feasibleregion(model)
    
    # Check if model has exactly two decision variables
    if num_variables(model) != 2
        throw(ArgumentError("Model is not 2-dimensional"))
    end
    
    if num_nonlinear_constraints(model) > 0
        throw(ArgumentError("Model has nonlinear constraints"))
    end
    

    # Convert the feasible region of the model into a polyhedron object
    LP_relax = copy(model)
    undo_relax = relax_integrality(LP_relax)
    ph = polyhedron(LP_relax)
    # undo_relax() # Undo the relaxation of integer variables

    # Obtain the extreme points and rays from the V-representation
    extremepoints = collect(points(vrep(ph)))
    extremerays = rays(vrep(ph))

    # From the extreme points, we can calculate the plotting window. 
    # We leave some empty space around the polyhedron itself 
    min_vals = minimum(hcat(extremepoints...), dims = 2)
    max_vals = maximum(hcat(extremepoints...), dims = 2)

    width = max_vals[1] - min_vals[1]
    height = max_vals[2] - min_vals[2]

    xlim = [min_vals[1]-width/10, max_vals[1]+width/4]
    ylim = [min_vals[2]-height/10, max_vals[2]+height/4]


    # If there are extreme rays, the feasible region is unbounded
    # and we add extra constraints outside the plotting area to enable plotting the feasible region
    if !isempty(extremerays)
        @info("Feasible region is unbounded")
        # TODO: Only add the necessary extra constraint(s)
        ph = ph ∩ HalfSpace(SparseArrays.SparseVector([1, 0]), xlim[2]+1) ∩ HalfSpace(SparseArrays.SparseVector([0, 1]), ylim[2]+1)
    end

    # Create empty plot
    plt = plot(xlim=xlim, ylim=ylim, xlabel="x1", ylabel="x2")
    
    # Plot the polyhedral feasible region
    plot!(ph, color=:blue, opacity=0.3, label="LP feasible region", legend=true)
    return plt
end

function add_linear_constraints!(plt, model)

    # Get the (exactly two) variables from the model
    x = all_variables(model)
    
    xlim = [lim for lim in xlims(plt)]
    
    # Add lines Ax=b to the plot 
    for con in all_constraints.(model, AffExpr, [MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64}])[1]
        # Get the constraint coefficients A and b
        a1 = normalized_coefficient(con, x[1])
        a2 = normalized_coefficient(con, x[2])
        b = normalized_rhs(con)
        # Plot the lines
        if a1!=0 && a2!=0
            # a1x1 + a2x2 = b => x2 = (b-a1x1)/a2
            plot!(plt, xlim, [(-a1*xlim .+ b)./ a2], label="$(name(con))") 
        elseif a1!=0 && a2==0
            vline!(plt, [b], label="$(name(con))")
        elseif a1==0 && a2!=0
            hline!(plt, [b], label="$(name(con))")
        end
    end
    # Plot x bounds
    if has_lower_bound(x[1])
        vline!(plt, [lower_bound(x[1])], label="x1=$(lower_bound(x[1]))")
    end
    if has_upper_bound(x[1])
        vline!(plt, [upper_bound(x[1])], label="x1=$(upper_bound(x[1]))")
    end
    if has_lower_bound(x[2])
        hline!(plt, [lower_bound(x[2])], label="x2=$(lower_bound(x[2]))")
    end
    if has_upper_bound(x[2])
        hline!(plt, [upper_bound(x[2])], label="x2=$(upper_bound(x[2]))")
    end
    return plt
end

function add_integrality!(plt, model)
    x = all_variables(model)
    LP_relax = copy(model)
    undo_relax = relax_integrality(LP_relax)
    ph = polyhedron(LP_relax)
    extremepoints = collect(points(vrep(ph)))

    xlim = [lim for lim in xlims(plt)]
    ylim = [lim for lim in ylims(plt)]
    
    # Check if both variables or one of them is integer and plot the integer feasible points if that is the case 
    if (is_integer(x[1]) || is_binary(x[1])) && (is_integer(x[2]) || is_binary(x[2]))
        # Both integer, collect the integer feasible points and plot them in a scatter plot

        feasiblepoints_x = []
        feasiblepoints_y = []
        for x in ceil(xlim[1]):floor(xlim[2]), y in ceil(ylim[1]):floor(ylim[2])
            if in([x,y], hrep(ph))
                push!(feasiblepoints_x, x) 
                push!(feasiblepoints_y, y)
            end
        end
        scatter!(plt, feasiblepoints_x, feasiblepoints_y, color=:black, label="Integer feasible region")
    elseif (is_integer(x[1]) || is_binary(x[1])) && !(is_integer(x[2]) || is_binary(x[2]))
        # Only x1 is integer, go through all the feasible values of x1 and for each one,
        # plot a line parallel to the x2-axis representing the feasible line segment

        flag = true # Used for plotting labels

        for x in ceil(xlim[1]):floor(xlim[2])
            # Get all extreme points that are on the left- or right-hand side of the given x-value
            p⁻ = filter(p -> p[1] <= x, extremepoints)
            p⁺ = filter(p -> p[1] >= x, extremepoints)

            # Make sure we don't have an empty collection of points, that would result in an error
            if !(isempty(p⁻) || isempty(p⁺))
                # At the value x, collect all y-values on the lines between the extreme points in p⁻ and p⁺
                # TODO: this code uses x2 and y somewhat interchangeably, might want to make it clearer
                yvals = []
                for p1 in p⁻, p2 in p⁺
                    if !(p2[1] ≈ p1[1])
                        yval = (p2[2] - p1[2])*(x - p1[1])/(p2[1] - p1[1]) + p1[2]
                        push!(yvals, yval)
                    end
                end

                # If we're plotting for the first time, flag=true and we give a label, otherwise no label
                lbl = (flag ? "Integer feasible region" : false) 
                if minimum(yvals) ≈ maximum(yvals) # Just a single point
                    scatter!(plt, [x], [minimum(yvals)], ms=3, color=:black, label=lbl)
                else # A line segment
                    plot!(plt, [x, x], [minimum(yvals), maximum(yvals)], lw = 3, color=:black, label=lbl)
                end
                flag = false 
            end
        end
    elseif !(is_integer(x[1]) || is_binary(x[1])) && (is_integer(x[2]) || is_binary(x[2]))
        
        # Only x2 is integer, go through all the feasible values of x2 and for each one,
        # plot a line parallel to the x1-axis representing the feasible line segment

        flag = true # Used for plotting labels

        for y in ceil(ylim[1]):floor(ylim[2])
            # Get all extreme points that are on the bottom- or top side of the given y-value 
            p⁻ = filter(p -> p[2] <= y, extremepoints)
            p⁺ = filter(p -> p[2] >= y, extremepoints)

            # Make sure we don't have an empty collection of points, that would result in an error
            if !(isempty(p⁻) || isempty(p⁺))
                # At the value y, collect all x-values on the lines between the extreme points in p⁻ and p⁺
                # TODO: this code uses x2 and y somewhat interchangeably, might want to make it clearer
                xvals = []
                for p1 in p⁻, p2 in p⁺
                    if !(p2[2] ≈ p1[2]) 
                        xval = (p2[1] - p1[1])*(y - p1[2])/(p2[2] - p1[2]) + p1[1]
                        push!(xvals, xval)
                    end
                end

                # If we're plotting for the first time, flag=true and we give a label, otherwise no label
                lbl = (flag ? "Integer feasible region" : false)
                if minimum(xvals) ≈ maximum(xvals) # Just a single point
                    scatter!(plt, [minimum(xvals)], [y], ms=3, color=:black, label=lbl)
                else # A line segment
                    plot!(plt, [minimum(xvals), maximum(xvals)], [y, y], lw = 3, color=:black, label=lbl)
                end
                flag = false
            end
        end
    end 
    return plt
end

function add_solution!(plt, model)
    # If the model has been solved and we have objective values,
    # get the optimal value and plot the line corresponding to that
    if has_values(model)
        x = all_variables(model)
        objfun = objective_function(model, AffExpr)
        c1 = coefficient(objfun, x[1])
        c2 = coefficient(objfun, x[2])
        z_opt = objective_value(model)
        xlim = [lim for lim in xlims(plt)]
        plot!(plt, xlim, [(-c1*xlim .+ z_opt)./ c2], linestyle=:dash, linewidth=2, color=:black, label="Objective function")
    end

    # Return the plot
    return plt
end

function plot_model(model)
   plt = plot_feasibleregion(model)
   add_linear_constraints!(plt, model)
   add_integrality!(plt, model)
   add_solution!(plt, model)
   return plt
end 
