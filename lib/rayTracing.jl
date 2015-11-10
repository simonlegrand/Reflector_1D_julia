""" Module containing functions for the ray tracer"""
module rayTracing

using lib
using PyPlot
using Grid

function reflection(points, I, s1)
	"""
	Computes directions s2 of reflected rays.
	Direct application of Snell & Descartes law.
	
	Parameters
	----------
	s1 : array (1,2)
		direction of incident rays
	points : array (N,2)
		points where the interpolant gradient is evaluated
	I : interpolant of the reflector
		
	Returns
	-------
	s2 : array (N,2)
		Direction of reflected rays(normalized).
	"""	
	gx = I[points]
#	gx = zeros(size(points)[1])
#	for i in 1:size(points)[1]
#		gx[i] = Grid.valgrad(I,points[i])[2]
#		#println(gx)
#	end
	ny = ones(points)
	n = cat(2, gx, -ny)
	for i in 1:size(n)[1]
		n[i,:] = n[i,:] / norm(n[i,:])
	end
	inner_product = s1[1] * n[:,1] + s1[2] * n[:,2]
	
	s2 = zeros(n)
	for i in 1:size(s2)[1]
		s2[i,:] = s1 - 2*inner_product[i]*n[i,:]
	end

	return s2
end

function ray_tracer{T<:Real}(s1::Array{T,2},source_density, s_box::Array{T,2}, t_box::Array{T,2}, interpol, base::Array{T,2}, niter::Int64=None):
	"""
	This function computes the simulation of reflection on the reflector
	and plot the image produced on the target screen.
	
	Parameters
	----------
	s1 : array (N,2)
		direction of incident rays
	s_box : tuple[2]
		enclosing square box of the source support
		[xmin, xmax]
	t_box : tuple[]
		enclosing square box of the target support
		[ymin, ymax]
	interpol : Cubic interpolant of the reflector
	base : [1]e_xi : Direct orthonormal 
		   basis of the target plan
		   [2]n_plan : Normal vector to the target plan
	niter : Number of loop of the ray tracer
	"""	
	e_xi = base[1,:]
	n_plan = base[2,:]
	
	M = None
	if niter == None
		niter = 1
	end
	for i in 1:niter
		nray = 200000
		# Generate source point uniformely
		#points = s_box[0] + (s_box[1] - s_box[0])*np.random.rand(nray)
		points = lib.inverse_transform(source_density, s_box[1], s_box[2], nray)
		#plt[:hist](points,100)
		s2 = reflection(points, interpol, s1)

		##### New polar coordinates #####
		# psi is the inclination with respect
		# to -n_plan
		d = norm(n_plan)
		scal_prod = zeros(nray)
		psi = zeros(nray)
		
		for k in 1:nray
			scal_prod[k] = -s2[k,1]*(n_plan/d)[1] - s2[k,2]*(n_plan/d)[2]
			psi[k] = acos(scal_prod[k])
		end
		
		J = s2[:,2] .< zeros(nray)
		psi[J] = -psi[J]
		
		##### Planar coordinates #####
		# computation of intersection of reflected rays
		# on the target plan and selection of points
		# inside t_box
		xi = d * tan(psi)
		xi_min = ones(nray) * t_box[1]
		xi_max = ones(nray) * t_box[2]
		K = (xi.<=xi_max) & (xi.>=xi_min)
		xi = xi[K]
		
		Miter = fill_sparse_vector(xi, t_box)
		if M == None
			M = Miter
		else
			M += Miter
		end
	
		println("it $(i): $(i*nray) rays thrown")
	end
	M = 255.0*M/maximum(M)
	
	return M
end


function fill_sparse_vector(x,box)

	h = box[2] - box[1] 
	n_linepix = 256
	nmesh = length(x)
	dx = h/n_linepix
	
	I = floor(Int64,(x-box[1])/dx + 1)
	M = zeros(n_linepix,1)
	for k in 1:length(I)
		#println(I[k])
		M[I[k]] += 1
	end
				
	return M
end

end
