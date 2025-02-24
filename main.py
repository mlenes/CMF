import numpy as np
from options import get_options
import make_grid
import assembly
import visuals

def main():
	# Include the options from the parser
	opts = get_options()
	mu = make_grid.init_mu(opts.N, opts.mu)
	A, b = assembly.assemble_matrix(opts.N, opts.L, mu, opts.p_grad, opts.bndry_type, opts.wallspeed, opts.wallgradient)
	u = np.linalg.solve(A, b)
	visuals.visualize(u)

if __name__ == "__main__":
	main()
