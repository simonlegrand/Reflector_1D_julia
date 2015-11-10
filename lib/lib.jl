module lib

using ApproXD
using Grid: CoordInterpGrid, BCnan, InterpQuadratic, InterpLinear, InterpCubic

function cost{T<:Real}(x::T,y::T)
	return 0.5*(x-y)^2
end

function cost_log{T<:Real}(x::Array{T},y::Array{T})
	return -log(1-dot(vec(x),vec(y)))
end

function cost_matrix{T<:Real}(x::Array{T}, y::Array{T})
	CM = zeros(T,(length(y),length(x)))
	for row = 1:length(y), col = 1:length(x)
		CM[row,col] = cost(x[col],y[row])
	end
	return CM
end

function cost_matrix_log{T<:Real}(x::Array{T}, y::Array{T})
	CM = zeros(T,(size(y)[1],size(x)[1]))
	println(x[5,:])
	for row = 1:size(y)[1], col = 1:size(x)[1]
		CM[row,col] = cost_log(x[col,:],y[row,:])
	end
	return CM
end


function derivates(x,u):
	"""
	This functions calculates u'(x). Datas must be equally spaced.
	"""
	n = length(u)
	if (length(x)!=n)
		println("Error: x and u must have the same number of elements!")
	end
	
	der = zeros(u)
	h = x[2] - x[1]
	for i in 2:n-1
		der[i] = (u[i+1] - u[i-1]) / (2*h)
	end
	der[1] = (u[2] - u[1]) / h
	der[n] = (u[n] - u[n-1]) / h
	
	return der
end

function interpolation_1d{T<:Real}(x_range, f::Array{T})
	I = CoordInterpGrid(x_range,f,BCnan,InterpQuadratic)
	return I
end

function inverse_transform(func, a, b, N, n_sample=2000)
	"""
	Invert Transform Method implementation
	Returns N random numbers following the func
	probability law.
	f has to be positive on [a,b]
	"""
	x = [linrange(a, b, n_sample)]

	F = ones(n_sample)
	dx = (b-a)/convert(Float64, n_sample)
	for i in 1:n_sample
		F[i] = quadgk(func, a, a + i*dx)[1]
	end
	F /= maximum(F)

	g = Array{Float64,1}[]          # grid container
	push!(g,F)
	interp = Lininterp(x,g)
	r = rand(N)
	invert = zeros(Float64,N)
	for i in 1:N
		invert[i] = getValue(interp,r[i])[1]
	end

	return invert
end

end
