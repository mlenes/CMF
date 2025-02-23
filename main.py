import numpy as np
from options import get_options
import simulate
import make_grid
import visuals
import assembly

def main():
	# Include the options from the parser
	opts = get_options()
	mu = make_grid.init_mu(opts.N, opts.mu)
	A, b = assembly.assemble_matrix(opts.N, opts.L, mu, opts.p_grad, 0)
	u = np.linalg.solve(A, b)
	print(u)

if __name__ == "__main__":
	main()
