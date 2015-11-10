module plots

using PyPlot


function plot_marginals(x, mu, y, nu)
	subplot(221)
	plot(x, mu, "b.", y, nu, "r.")
	title("Marginals : blue=source  red=target")
end

function plot_polar_marginal(theta, mu)
	fig = figure("Polar density",figsize=(10,10)) # Create a new figure
	ax = axes(polar="true") # Create a polar axis
	title("Polar density")
	p = plot(theta,mu,linestyle="-",marker="None")
end


function plot_transport_plan(x, y, Gamma)
	if issparse(Gamma)
		Gamma = full(Gamma)
	end
	subplot(223)
	title("Transport plan")
	xlabel("x")
	ylabel("y")
	# Plot Gamma transpose because imshow inverts abs/ord
	imshow(Gamma, origin="lower", extent=[minimum(x),maximum(x),minimum(y),maximum(y)])
end


function plot_potential(x, potential)
	subplot(222)
	xlabel("x")
	ylabel("psi")
	title("Kantorovich potential")
	plot(x,potential,"r.", linewidth=1)
end


function plot_everything(x, mu, y, nu, potential, Gamma)
	"""
	Plot all the graphs on the same figure
	"""
	plot_marginals(x, mu, y, nu)
	plot_potential(x, potential)
	plot_transport_plan(x, y, Gamma)
end

end
