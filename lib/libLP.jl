module libLP

using Grid: CoordInterpGrid, BCnan, InterpLinear
using JuMP
using Clp

function cost{T<:Real}(x::T,y::T)
	return 0.5*(x-y)^2
end


function cost_matrix{T<:Real}(x::Array{T}, y::Array{T})
	CM = zeros(T,(length(y),length(x)))
	for row = 1:length(y), col = 1:length(x)
		CM[row,col] = cost(x[col],y[row])
	end
	return CM
end


function get_threshold{T<:Real}(Gamma::AbstractMatrix{T})
	"""
	Determines the threshold under which elements
	of Gamma are considered null.
	
	See 'A Numerical Method to solve Optimal Transport
	Problems with Coulomb Cost' from Benamou, Carlier,
	Nenna
	"""
	# Maximums over each dimension of Gamma
	Max = T[]
	Min = T[]
	
	for i = 1:ndims(Gamma)
		Max = maximum(Gamma,i)
		push!(Min, minimum(Max))
	end
		
	coeff = 0.5
	threshold = coeff * minimum(Min)
	return threshold
end


function get_nb_active_elements{T<:Real}(M::AbstractMatrix{T}, threshold)
	"""
	Return the number of elements in M which are
	superior to the threshold.
	"""
	threshold = threshold * ones(M)
	J = M .> threshold
	return sum(J)
end


function linear_interpolation{T<:Real}(x_range,y_range, Gamma::Matrix{T})
	I = CoordInterpGrid((x_range,y_range),Gamma,BCnan,InterpLinear)
	return I
end


function eval_interpolant{T<:Real}(x::Array{T}, y::Array{T}, I)
	z = zeros(T,(length(y),length(x)))
	for row = 1:length(y), col = 1:length(x)
		z[row,col]= I[x[col],y[row]]
	end
	return z
end


function to_sparse!{T<:Real}(M::Matrix{T}, threshold::T)
	zero = zeros(M)
	J = M .< threshold
	M[J] = zero[J]
	return sparse(M)
end


function fill_data(x::Array{Float64}, y::Array{Float64}, Gamma::SparseMatrixCSC{Float64,Int64}, eps::Float64)
	row = findnz(Gamma)[1]	# row indicies of non nul elements
	col = findnz(Gamma)[2]	# column indicies of non nul elements
	data = ones(Float64, length(row))
	for i in 1:length(row)
		data[i] = exp(-cost(x[col[i]],y[row[i]]) / eps)
	end
	return data
end

function solve_LP{T<:Real}(mu::Array{T},nu::Array{T},C::Matrix{T})
	Nx = size(mu)[1]
	Ny = size(nu)[1]
	
	m = Model(solver=ClpSolver(LogLevel=4))
	@defVar(m,phi[1:Nx])
	@defVar(m,psi[1:Ny])

	for i in 1:Ny, j in 1:Nx
		@addConstraint(m, phi[j] + psi[i] <= C[i,j])
	end

	@setObjective(m, Max, sum{phi[i]*mu[i], i=1:Nx} + sum{psi[j]*nu[j], j=1:Ny})
	stat = solve(m)
	println("Solve status: ", stat)
	
	return getValue(phi),getValue(psi)
end

end
