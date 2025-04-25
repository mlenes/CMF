import numpy as np
import options
import assembly
import visuals
import tools
import iterators
import matplotlib.pyplot as plt

def main():
	# Include the options from the parser
	opts = options.get_options()

	if opts.global_type == 'flowrate':
		p_grad = -4*opts.flowrate/opts.L*opts.mu
		u_ref = opts.flowrate / (opts.rho_ref*opts.L**2)
	else:
		p_grad = opts.p_grad
		u_ref = np.sqrt(-p_grad*opts.L/(2*opts.rho_ref)) 
    
	# Compute viscosity at the faces of the cells and make dimensionless
	mu_faces = tools.init_mu_faces(opts.N, opts.mu)/(u_ref*opts.L*opts.rho_ref)

	# Make p_grad and flowrate dimensionless
	p_grad /= (2*opts.rho_ref*u_ref)
	flowrate = opts.flowrate/(u_ref*opts.L**2*opts.rho_ref)

	# Assemble the A matrix for calculating
	A = assembly.assemble_A(opts.N, opts.L, mu_faces, False, 1)

	# Assemble the b vector
	b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad)

	# Compute u from linear equation Au = b
	u = np.linalg.solve(A, b)
    
	if opts.global_type=="flowrate":
		u, p_grad = iterators.iter_flowrate(opts.N, opts.L, mu_faces, p_grad, A, u, flowrate)

	if opts.flow_type == 'turbulent':
		u, p_grad = iterators.iter_u(opts.iterations, u, opts.N, opts.L, opts.Ks, mu_faces, opts.rel_factor, p_grad, opts.global_type, flowrate)
        
	visuals.visualize(u*u_ref)

if __name__ == "__main__":
	main()
