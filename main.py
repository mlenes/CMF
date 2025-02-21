from options import get_options
import simulate
import make_grid
import visuals

def main():
	# Include the options from the parser
	opts = get_options()

	dy_arr = make_grid.vertical_dist(opts.N, opts.dy, opts.Nlayers, opts.layerfactor)
	init_u = make_grid.init_u_dirichlet(opts.N, opts.flowrate, opts.wallspeed, dy_arr)
	init_mu = make_grid.init_mu(opts.N, opts.mu)
	u = simulate.simulation_pressure(init_u, init_mu, dy_arr, opts.iterations, opts.p_grad)
	print(u)

if __name__ == "__main__":
	main()
#test
