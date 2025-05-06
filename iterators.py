import numpy as np
import assembly
import tools
import matplotlib.pyplot as plt

def iter_flowrate(N, L, mu_faces, p_grad, A, u, flowrate):

	# Dummy initial value
	error_flowrate = 10
	
	while (np.abs(error_flowrate) > 1e-3*flowrate):
		# Find corrected pressure gradient using calculated flowrate and assuming linear relation
		correction_flowrate =  (np.sum(u)*L/N) / flowrate
		p_grad /= correction_flowrate
		
		b = assembly.assemble_b(N, L, mu_faces, p_grad)
		u = np.linalg.solve(A, b)
		
		error_flowrate = np.sum(u)*L/N - flowrate
		
	return u, p_grad
	
def iter_u(iterations, u, N, L, Ks, mu_faces, rel_factor, p_grad, global_type, flowrate):
    for i in range(iterations):
        #each iteration first loop until velocity profile matches eddy viscosity
        #if flowrate is given, then check converged velocity profile with flowrate and if needed correct
        u_corr = 10*np.mean(u)
        while np.abs(np.mean(u_corr)) > np.abs(1e-5*np.mean(u)):
            #u_tau = u[1] * 1/0.41 * np.log(1/0.031* tools.get_dy(opts.N, opts.L)/opts.Ks +1) #previous value
            
            # Get wall constant
            u_tau = u[1] * 0.41 / (np.log((tools.get_dy(N, L)*32.6/2)/Ks))
            tau_w = u_tau**2 #non-dimensional, so rho=1
            wall_constant = tau_w/u[1]
            
            # Eddy viscosity
            l_m = tools.calc_mixing_length(N, L, u_tau, mu_faces[1])
            mu_eff = mu_faces + l_m**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(N, L)
            
            # Get turbulent A and b
            A = assembly.assemble_A(N, L, mu_eff, True, wall_constant)
            b = assembly.assemble_b(N, L, mu_eff, p_grad)
            u_new = np.linalg.solve(A, b)
            u_corr = u_new - u
            u += rel_factor*u_corr
            
        if global_type=="flowrate":
        	# Make u correct with flowrate
            u, p_grad = iter_flowrate(N, L, mu_faces, p_grad, A, u, flowrate)
    # plt.plot(u[1:200],np.linspace(0,0.5,len(u[1:200])))    
    # plt.yscale('log')
    # plt.xlabel("u")
    # plt.ylabel("nodenumber log")
    return u, p_grad

def particle(y0,u,dt,tracktime,mu,D,M):
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
