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
            
            # Eddy viscosity
            l_m = tools.calc_mixing_length(N, L, u_tau, mu_faces[1])
            uprime = l_m * np.abs(u[:-1] - u[1:])/tools.get_dy(N, L)
            mu_eff = mu_faces + uprime*l_m
            
            # Get turbulent A and b
            A = assembly.assemble_A(N, L, mu_eff, tau_w)
            b = assembly.assemble_b(N, L, mu_eff, p_grad, tau_w)
            u_new = np.linalg.solve(A, b)
            u_corr = u_new - u
            u += rel_factor*u_corr
            
    return u
