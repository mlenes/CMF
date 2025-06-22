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
			u, u_prime = iterators.iter_u(opts.iterations, u, opts.N, opts.L, opts.Ks, mu_faces, opts.rel_factor, p_grad)
        
			print("mean flow found, now adding particles for time", t)
        
		#gehardcode oeps
		y0 = np.array([0.01, 0.3, 0.5]) #non-dimensional height where we insert particle in channel
		v0 = np.array([4, 0])  #non-dimensional starting velocity
		
		M = opts.M / (opts.rho_ref*opts.L**3)  #non-dimensional weight of particle
		g = 9.81 * opts.L/ u_ref**2  #non-dimensional gravitational acceleration
		x_list, y_list = iterators.particle(y0,v0,u,opts.dt,opts.tracktime,mu_faces[1],opts.D,M,g,opts.N,opts.L,u_prime)
		
		visuals.particles(x_list*opts.L,y_list*opts.L)
		print(u[1000])
		u_list.append(u*u_ref)
	u_list = np.array(u_list)
	visuals.visualize(u_list)

if __name__ == "__main__":
	main()
