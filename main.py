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
	
	u_list = []
	for t in range(opts.sim_time):

		p_grad = tools.get_sin_p(opts.p_grad, opts.pressure_period, t)
		u_ref = 6
		
		# Compute viscosity at the faces of the cells and make dimensionless
		mu_faces = tools.init_mu_faces(opts.N, opts.mu)

		# Make p_grad and viscosity dimensionless
		p_grad /= ((opts.rho_ref*u_ref**2)/opts.L)
		mu_faces /= (u_ref*opts.L*opts.rho_ref)

		# Assemble the A matrix for calculating
		A = assembly.assemble_A(opts.N, opts.L, mu_faces)

		# Assemble the b vector
		b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad)

		# Compute u from linear equation Au = b
		u = np.linalg.solve(A, b)
		
		if opts.flow_type == 'turbulent':
			u, u_prime, u_tau = iterators.iter_u(opts.iterations, u, opts.N, opts.L, opts.Ks, mu_faces, opts.rel_factor, p_grad,1)
        
			print("mean flow found, now adding particles for time", t)
        
		y_plus = tools.get_dy(opts.N,opts.L)*u_tau*u_ref*opts.rho_ref/opts.mu/opts.L
# 		print("y_plus=",y_plus)
		
		#gehardcode oeps
		y0 = np.array([0.1,0.3,0.5]) #non-dimensional height where we insert particle in channel
		v0 = np.array([0, 0])  #non-dimensional starting velocity
		
		M = opts.rho_p / opts.rho_ref *np.pi*(opts.D/opts.L)**3/6  #non-dimensional weight of particle
		g = 9.81 * opts.L/ u_ref**2  #non-dimensional gravitational acceleration
		x_list, y_list = iterators.particle(y0,v0,u,opts.dt,opts.tracktime,mu_faces[1],opts.D,M,g,opts.N,opts.L,u_prime)
		visuals.particles(x_list*opts.L,y_list*opts.L)
		u_list.append(u*u_ref)
		
		
	#make a log plot of the velocity in plus units
# 	uarray = np.array(u_list[0])
# 	uplus = (uarray/u_ref)/u_tau
# 	conversion = u_tau*u_ref*opts.rho_ref/opts.mu
# 	n=50 #how many points you want to plot
# 	visuals.wallprofile(opts.L,opts.N,n,conversion,uplus)

	#make a plot of the mixing length function
# 	plt.figure()
# 	channel=np.linspace(0,1,10000)
# 	mixing_length = tools.calc_mixing_length(N=10000, L=1, u_tau=0.01, mu=4e-7)
# 	y_plus = (channel * 0.01) / 4e-7
# 	A_plus = 25
# 	D = (1-np.exp(-y_plus/A_plus))
# 	mixing_lengthD = mixing_length*D
# 	plt.plot(channel[:-9000],mixing_length[:-9000],label="no damping")
# 	plt.plot(channel[:-9000],mixing_lengthD[:-9000],label="van Driest damping")
# 	plt.xlabel("y/H")
# 	plt.ylabel("$l_m$/H")
# 	plt.xscale('log')
# 	plt.yscale('log')
# 	plt.ylim(1e-6,1e-1)
# 	plt.xlim(1e-4,1e-1)
# 	plt.legend()
# 	plt.show()
		
	u_list = np.array(u_list)
	visuals.visualize(u_list)
	
	

if __name__ == "__main__":
	main()
