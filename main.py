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
    A = assembly.assemble_A(opts.N, opts.L, mu_faces, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top, False, 1)
    
    # Assemble the b vector
    b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
    
    # Compute u from linear equation Au = b
    u = np.linalg.solve(A, b)
    
    if opts.global_type=="flowrate":
    	u, p_grad = iterators.iter_flowrate(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top, A, u, flowrate)

    if opts.flow_type == 'turbulent':
        for i in range(opts.iterations):
            #each iteration first loop until velocity profile matches eddy viscosity
            #if flowrate is given, then check converged velocity profile with flowrate and if needed correct
            u_corr = 10
            while np.abs(np.mean(u_corr)) > 1e-5*np.mean(u):
                #calculate u_tau using newton raphson: u[1]/u_tau = 1/K * ln (y[1] * u_tau / opts.mu*opts.rho_ref) + 1/K * ln(32.6 / (u_tau*opts.Ks/ opts.mu*opts.rho_ref))
                #u_tau = u[1] * 1/0.41 * np.log(1/0.031* tools.get_dy(opts.N, opts.L)/opts.Ks +1) #some value for now 
                u_tau = u[1] * 0.41 / (np.log((tools.get_dy(opts.N, opts.L)*32.6/2)/opts.Ks))
                tau_w = u_tau**2 #non-dimensional, so rho=1
                wall_constant = tau_w/u[1]
                mu_eff = mu_faces + tools.calc_mixing_length(opts.N, opts.L, opts.bndry_bot, opts.bndry_top)**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(opts.N,opts.L)
                A = assembly.assemble_A(opts.N, opts.L, mu_eff, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top, True, wall_constant)
                b = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                u_new = np.linalg.solve(A, b)
                u_corr = u_new - u
                u += opts.rel_factor*u_corr
                
            if opts.global_type=="flowrate":
                u, p_grad = iterators.iter_flowrate(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top, A, u, flowrate)
        
    visuals.visualize(u*u_ref)

if __name__ == "__main__":
	main()
