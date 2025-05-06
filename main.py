import numpy as np
import options
import assembly
import visuals
import tools
import iterators

def main():
	# Include the options from the parser
	opts = options.get_options()
	
	u_list = []
	for t in range(opts.sim_time):

		p_grad = tools.get_sin_p(opts.p_grad, opts.pressure_period, t)
		u_ref = tools.get_u_ref(p_grad, opts.L, opts.rho_ref)
		
		# Compute viscosity at the faces of the cells and make dimensionless
		mu_faces = tools.init_mu_faces(opts.N, opts.mu)

		# Make p_grad and viscosity dimensionless
		p_grad /= 2*opts.rho_ref*u_ref
		mu_faces /= u_ref*opts.L*opts.rho_ref

		# Assemble the A matrix for calculating
		A = assembly.assemble_A(opts.N, opts.L, mu_faces)

		# Assemble the b vector
		b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad)

		# Compute u from linear equation Au = b
		u = np.linalg.solve(A, b)
		
		if opts.flow_type == 'turbulent':
			u = iterators.iter_u(opts.iterations, u, opts.N, opts.L, opts.Ks, mu_faces, opts.rel_factor, p_grad)
        
		#gehardcode oeps
		y0 = np.array([int(0.01*opts.N), int(0.3*opts.N), int(0.5*opts.N)]) #nodenumber where we insert particle in channel
		
		M = opts.M / (opts.rho_ref*opts.L**3)  #non-dimensional weight of particle
		x_list = iterators.particle(y0,u,opts.dt,opts.tracktime,mu_faces[1],opts.D,M)
		
		visuals.particles(x_list*opts.L,y0*opts.L/opts.N)
		
		u_list.append(u*u_ref)
	u_list = np.array(u_list)
	visuals.visualize(u_list)

if __name__ == "__main__":
	main()
