module libIPFP

using Grid: CoordInterpGrid, BCnan, InterpLinear

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
	I = CoordInterpGrid((y_range,x_range),Gamma,BCnan,InterpLinear)
	return I
end


function eval_interpolant{T<:Real}(x::Array{T}, y::Array{T}, I)
	Nx = length(x)
	Ny = length(y)
	z = zeros(T,Ny,Nx)
	for i=1:Ny, j=1:Nx
		z[i,j]= I[y[i],x[j]]
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


function solve_IPFP{T<:Real}(mu::Vector{T}, nu::Vector{T}, C::Matrix{T}, epsilon::T)

	G_bar = exp(-C/epsilon)
	b = ones(mu)
	a = ones(nu)
	error = 1
	error_min = 1e-3
	count = 0
	
	while error > error_min

		b = mu./(G_bar'*a)
		an = nu./(G_bar*b)
		#println(b)
		error = sum(abs(an-a))/sum(a)
		println("error at step ", count, " = ", error)
		count = count + 1
		a = an
	end
	
	G = diagm(a)*G_bar*diagm(b)
	
	return G,a,b
end


function solve_IPFP{T<:Real}(mu::Vector{T}, nu::Vector{T}, C::SparseMatrixCSC{T,Int64}, epsilon::T)

	G_bar = exp(-C/epsilon)
	b = ones(mu)
	a = ones(nu)
	error = 1
	error_min = 1e-3
	count = 0
	
	while error > error_min
	
		b = mu./(G_bar'*a)
		an = nu./(G_bar*b)
		
		error = 1e-4 #sum(abs(an-a))/sum(a)
		println("error at step ", count, " = ", error)
		count = count + 1
		a = an
	end
	G = spdiagm(a)*G_bar*spdiagm(b)

	return G,a,b
end


function solve_IPFP{T<:Real}(x::Vector{T},y::Vector{T},mu::Vector{T},nu::Vector{T},epsilon::T)

	b = ones(mu)
	a = ones(nu)
	error = 1
	error_min = 1e-3
	count = 0

	while error > error_min
		
		#sa = zeros(x)
		col_gbar = zeros(y)
		for j in 1:length(x)
			for i in 1:length(y)
				col_gbar[i] = exp(-1./epsilon*cost(x[j],y[i]))
			end
			#sa[j] = dot(col_gbar,a)
			b[j] = mu[j] / dot(col_gbar,a)
		end
		
		an = zeros(a)
		#sb = zeros(y)
		lin_gbar = zeros(x)
		for i in 1:length(y)
			for j in 1:length(x)
				lin_gbar[j] = exp(-1./epsilon*cost(x[j],y[i]))
			end
			#sb[i] = dot(lin_gbar,b)
			an[i] = nu[i] / dot(lin_gbar,b)
		end
		
		error = sum(abs(an-a))/sum(a)
		println("error at step ", count, " = ", error)
		count = count + 1
		a = an
	end
	
	G = diagm(a)*G_bar*diagm(b)
	
	return G,a,b
end
end
