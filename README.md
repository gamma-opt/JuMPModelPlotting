# JuMPModelPlotting.jl
 Plotting functions for linear JuMP models with exactly two variables.

The main function is `plot_model(model)`, which takes a linear [JuMP](https://jump.dev/) model with exactly two variables and returns a graphical representation of the model. This graphical representation consists of:
- The LP feasible region, using the [Polyhedra](https://github.com/JuliaPolyhedra/Polyhedra.jl) package
- The affine constraints $Ax=b$
- The integer feasible region, if one or both of the variables are integer/binary
- The objective level curve corresponding to the optimal solution, if the model has been solved and has a solution available

The implementation of the package has been split into functions corresponding to the four bullet points above, in order to make the code somewhat modular. The `plot_model(model)` -function calls these four functions, resulting in the final plot. 

## Example:

Given a model 

$$
\begin{align*}
    \max \ & x_1 + 2x_2 \\
    \text{s.t. } & 6x_1 + 4x_2 && \leq 24 \\
    & x_1 + 2x_2 && \leq 6 \\
    & -x_1 + x_2 && \leq 1 \\
    & x_2 && \leq 2 \\
    & x_1 \geq 0 \\
    & x_2 \text{ integer},
\end{align*}
$$

the resulting plot is

![example plot](https://github.com/solliolli/JuMPModelPlotting/blob/main/examples/paintfactory_int.png)
