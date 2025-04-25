import numpy as np
import assembly
import tools

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
            mu_eff = mu_faces + tools.calc_mixing_length(N, L)**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(N, L)
            
            # Get turbulent A and b
            A = assembly.assemble_A(N, L, mu_eff, True, wall_constant)
            b = assembly.assemble_b(N, L, mu_eff, p_grad)
            u_new = np.linalg.solve(A, b)
            u_corr = u_new - u
            u += rel_factor*u_corr
            
        if global_type=="flowrate":
        	# Make u correct with flowrate
            u, p_grad = iter_flowrate(N, L, mu_faces, p_grad, A, u, flowrate)
            
    return u, p_grad
