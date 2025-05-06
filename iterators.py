import numpy as np
import assembly
import tools

	
def iter_u(iterations, u, N, L, Ks, mu_faces, rel_factor, p_grad):
    for i in range(iterations):
        #each iteration first loop until velocity profile matches eddy viscosity
        u_corr = 10*np.mean(u)
        while np.abs(np.mean(u_corr)) > np.abs(1e-5*np.mean(u)):
            
            # Get wall constant
            u_tau = u[1] * 0.41 / (np.log((tools.get_dy(N, L)*32.6/2)/Ks))
            tau_w = u_tau**2 #non-dimensional, so rho=1
            wall_constant = tau_w/u[1]
            
            # Eddy viscosity
            l_m = tools.calc_mixing_length(N, L, u_tau, mu_faces[1])
            mu_eff = mu_faces + l_m**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(N, L)
            
            # Get turbulent A and b
            A = assembly.assemble_A(N, L, mu_eff, wall_constant)
            b = assembly.assemble_b(N, L, mu_eff, p_grad)
            u_new = np.linalg.solve(A, b)
            u_corr = u_new - u
            u += rel_factor*u_corr
            
    return u

def particle(y0, u, dt, tracktime, mu, D, M):
	x_list = []
	
	#loop for every particle we inserted
	for n in range(len(y0)):
		x = np.zeros(int(tracktime/dt))
		v = np.zeros(int(tracktime/dt))
	
		#update position and velocity by calculating force
		for i in range(1,len(x)):
			x[i] = x[i-1] + v[i-1]*dt
			f=0.1 #friction factor, heb ik random gekozen nu
			F = f*3*np.pi*mu*D*(u[y0[n]]-v[i-1]) #now only drag force, all input is non-dimensional
			v[i] = v[i-1] + 1/M * F * dt
		x_list.append(x)
		
	return np.array(x_list) #this is the position of all particles over time
