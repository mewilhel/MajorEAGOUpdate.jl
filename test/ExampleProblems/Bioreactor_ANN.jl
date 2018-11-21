# Box constraints for input variables
xLBD = [0.623, 0.093, 0.259, 6.56, 1114,  0.013, 0.127, 0.004]
xUBD = [5.89,  0.5,   1,     90,   25000, 0.149, 0.889, 0.049]

# Weights associated with the hidden layer
W = [ 0.54 -1.97 0.09 -2.14 1.01 -0.58 0.45 0.26;
     -0.81 -0.74 0.63 -1.60 -0.56 -1.05 1.23 0.93;
     -0.11 -0.38 -1.19 0.43 1.21 2.78 -0.06 0.40]

# Weights associated with the output layer
D = [-0.91 0.11 0.52]

# Bias associated with the hidden layer
B1 = [-2.698 0.012 2.926]

# Bias associated with the output layer
B2 = -0.46

# Model construction
model = Model(with_optimizer(EAGO.Optimizer))
@variable(model, xLBD[i] <= x[i in 1:8] <= xUBD[i])
@NLexpression(model, prop[i=1:3], B1[i] + sum(W[i,j]*x[j] for j in 1:8))
@NLobjective(model, Max, B2 + sum(D[i]*(2/(1+exp(-2*prop[i]))) for i in 1:3))

# Solves the model
optimize!(model)

# Access functions for the solution
ObjectiveValue = JuMP.objective_value(model)
Feasibility = JuMP.primal_status(model)
Solution = JuMP.value.(x)
