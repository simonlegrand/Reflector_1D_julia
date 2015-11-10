module geometry

function planar_to_gradient(eta, xi, base, s1=None)
	"""
	This function computes the surface derivatives of the reflector
	for incident rays s1 and impact points of reflected rays in (eta,xi)
	Parameters
	----------
	eta : 1D array
		Coordinate eta on the target plane
	xi : 1D array
		Coordinate xi on the target plane
	base : [0]e_eta , [1]e_xi : Direct orthonormal 
		   basis of the target plan
		   [2]n_plan : Normal vector to the target plan
		   Its norm equals distance from plan to reflector.
	s1 : (1,3) array
		Incident ray direction

	Returns
	-------
	p,q : 1D arrays
		surface derivatives of the reflector
		
	See Also
	--------
	Inverse Methods for Illumination Optics, Corien Prins, chapter 5.3.1
    """
	e_eta = base[1]
	e_xi =base[2]
	n = base[3]
	
	if s1==None
		s1 = [0.,0.,1.]	
	else
		s1 = s1 / norm(s1)
	end


	# Distance target plan/reflector
	d = norm(n)
	n = n / d

	# Reflected rays
	# The reflector is considered ponctual and
	# as the origin of the coordinate system
	s2 = zeros((len(eta),3))
	s2[:,1] = eta*e_eta[1] + xi*e_xi[1] - d*n[1]
	s2[:,2] = eta*e_eta[2] + xi*e_xi[2] - d*n[2]
	s2[:,3] = eta*e_eta[3] + xi*e_xi[3] - d*n[3]

	s2 = s2 / norm(s2)

	p = -(s2[:,1] - s1[1])/(s2[:,3] - s1[3])
	q = -(s2[:,2] - s1[2])/(s2[:,3] - s1[3])

	return p, q
end


function planar_to_gradient_boundaries(target_plane_box,target_plane_base)
	"""
	Computes the p-domain boundaries from the target box
	"""
	e_xi = target_plane_base[1,:]
	n_plan = target_plane_base[2,:]
	
	s1 = [0. 1.]

	s2_bound = zeros(Float64, (2,2))
	s2_bound[1,1] = target_plane_box[1]*e_xi[1] - n_plan[1]
	s2_bound[1,2] = target_plane_box[1]*e_xi[2] - n_plan[2]
	s2_bound[1,:] /= norm(s2_bound[1,:])
	s2_bound[2,1] = target_plane_box[2]*e_xi[1] - n_plan[1]
	s2_bound[2,2] = target_plane_box[2]*e_xi[2] - n_plan[2]
	s2_bound[2,:] /= norm(s2_bound[2,:])

	pmin = -(s2_bound[1,1] - s1[1])/(s2_bound[1,2] - s1[2])
	pmax = -(s2_bound[2,1] - s1[1])/(s2_bound[2,2] - s1[2])
	
	return pmin,pmax
end
	
function planar_to_gradient(xi, base, s1=None)
	"""
	This function computes the surface derivative of the 1D reflector
	for incident rays s1 and impact points of reflected rays xi.
	Parameters
	----------
	xi : 1D array
		Coordinate xi on the target
	base : [1]e_xi : Basis of the target
		   [2]n_plan : Normal vector to the target
		   Its norm equals distance from plan to reflector.
	s1 : (1,2) array
		Incident rays direction

	Returns
	-------
	p : 1D arrays
		surface derivative of the reflector
		
	See Also
	--------
	Inverse Methods for Illumination Optics, Corien Prins, chapter 5.3.1
    """
	e_xi = base[1,:]
	n = base[2,:]
	println(n)
	if s1==None
		s1 = [0.,1.]	
	else
		s1 = s1 / norm(s1)
	end

	# Distance target/reflector
	d = norm(n)
	n = n / d

	# Reflected rays
	# The reflector is considered ponctual and
	# at the origin of the coordinate system
	s2 = zeros(Float64,(length(xi),2))
	for i=1:length(xi)
		s2[i,1] = xi[i]*e_xi[1] - d*n[1]
		s2[i,2] = xi[i]*e_xi[2] - d*n[2]
	end

	p = zeros(size(s2,1))
	for i=1:size(s2,1)
		s2[i,:] = s2[i,:] / norm(s2[i,:])
		p[i] = -(s2[i,1] - s1[1])/(s2[i,2] - s1[2])
	end
	
	return p
end

function planar_to_directions_boundaries(target_plane_box,target_plane_base)
	"""
	Computes the rays angles boundaries from the target box
	"""
	e_xi = target_plane_base[1,:]
	n_plan = target_plane_base[2,:]

#	s2_bound = zeros(Float64, (2,2))
#	s2_bound[1,1] = target_plane_box[1]*e_xi[1] - n_plan[1]
#	s2_bound[1,2] = target_plane_box[1]*e_xi[2] - n_plan[2]
#	s2_bound[1,:] /= norm(s2_bound[1,:])
#	s2_bound[2,1] = target_plane_box[2]*e_xi[1] - n_plan[1]
#	s2_bound[2,2] = target_plane_box[2]*e_xi[2] - n_plan[2]
#	s2_bound[2,:] /= norm(s2_bound[2,:])

	psi_min = pi/2. - atan(target_plane_box[2]/norm(n_plan))
	psi_max = pi/2. - atan(target_plane_box[1]/norm(n_plan))
	return psi_min,psi_max
end

function density_change_variable(p,target_plane_base,target_plane_density)

	e_xi = target_plane_base[1,:]
	n_plan = target_plane_base[2,:]
	d = norm(n_plan)
	n_plan = n_plan / d
	
	s2 = zeros((length(p),2))
	s2[:,1] = 1./(p.*p + 1) .* 2.*p
	s2[:,2] = 1./(p.*p + 1) .* (p.*p - 1)
	
	xi = -d * (e_xi[1]*s2[:,1] + e_xi[2]*s2[:,2])
	xi ./= n_plan[1]*s2[:,1] + n_plan[2]*s2[:,2]
	
	G_tilde = (xi.*xi + d*d)/d .* target_plane_density(xi)
	g = 2 * G_tilde ./ (p.*p + 1)
	
	return g
end
end
