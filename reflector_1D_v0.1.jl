using PyPlot
# Path to IPFP libraries
reload("./lib/libIPFP")
reload("./lib/plots")
reload("./lib/libLP")
# Path to geometry function related to tehe reflector pb
reload("./lib/geometry.jl")
reload("./lib/lib.jl")
reload("./lib/rayTracing.jl")

#### Marginals definition ####
function source_density(x)	
	source = exp(-x.*x)
	return source
end


function target_plane_density{T<:Real}(xi::Array{T})
	target = ones(xi)
	for i = 1:length(xi)
		target[i] = 0.5*exp(-10*(xi[i]-1.5)*(xi[i]-1.5)) + exp(-(xi[i]+1.5)*(xi[i]+1.5))
	end
	return target
end


I = true # IPFP resolution
LP = false # Linear programmation resolution


#### Parameters ####
Nx = 100
source_box = [-2. 2.]
x_range = linrange(source_box[1],source_box[2],Nx)
x = [x_range]
Np = 200
target_plane_box = [-4. 4.]
e_xi = [0. 1.]
n_plan = [-20. 0.]
target_plane_base = [e_xi;n_plan]
# For plot
xi = linspace(target_plane_box[1], target_plane_box[2], Np)
tar_plane_dens = target_plane_density(xi)
tar_plane_dens ./= sum(tar_plane_dens) 

#### Boundaries of the p domain ####
(pmin,pmax) = geometry.planar_to_gradient_boundaries(
						target_plane_box,target_plane_base)
p = linspace(pmin,pmax,Np)
g = geometry.density_change_variable(p,target_plane_base,target_plane_density)


# Marginals
mu = source_density(x)
mu ./= sum(mu)
nu = g
nu ./= sum(nu)

C = lib.cost_matrix(x,p)

if I
	#### IPFP resolution ####
	tic()
	epsilon = 0.006
	Gamma,a,b = libIPFP.solve_IPFP(mu,nu,C,epsilon)
	phi = epsilon*log(b)
	psi = epsilon*log(a)
	println("OT resolution")
	toc()
	plots.plot_everything(x,mu,xi,tar_plane_dens,u,Gamma)
end

if LP
	#### Linear programming resolution ####
	phi,psi = libLP.solve_LP(mu,nu,C)
	println("OT resolution")
end

#u0 = 0.5*x.*x - phi0
u = 0.5*x.*x - phi
#plot(x,u,"b.",x,u0,"r.")
grad_u = lib.derivates(x, u)
der_interpol = lib.interpolation_1d(x_range, grad_u)

#x_prime = [linrange(source_box[1],source_box[2],10000)]
#grad = zeros(10000)
#for i in 1:10000
#	grad[i] = der_interpol[x_prime[i]]
#end
#plot(x_prime,grad,x,u)

#### Ray tracing ####
s1 = [0. 1.]
M = rayTracing.ray_tracer(s1,source_density, source_box, target_plane_box, der_interpol, target_plane_base, 10)
x = linspace(target_plane_box[1],target_plane_box[2],size(M)[1])
subplot(224)
plot(x, M, "b-", ms=2)

