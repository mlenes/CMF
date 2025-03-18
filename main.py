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
    else:
        p_grad = opts.p_grad
    u_tau = np.sqrt(-p_grad*opts.L/(2*opts.rho_ref)) #assuming rho=1000
    
    # Compute viscosity at the faces of the cells dimensionless
    mu_faces = tools.init_mu_faces(opts.N, opts.mu)/(u_tau*opts.L*opts.rho_ref)
    
    # Make p_grad and flowrate dimensionless
    p_grad /= (2*opts.rho_ref*u_tau)
    flowrate = opts.flowrate/(u_tau*opts.L**2*opts.rho_ref)
    
    # Assemble the A matrix for calculating
    A = assembly.assemble_A(opts.N, opts.L, mu_faces, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
    
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
                #assume rho=1, so rho=rho_water
                mu_eff = mu_faces + tools.calc_mixing_length(opts.N, opts.L, opts.bndry_bot, opts.bndry_top)**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(opts.N,opts.L)
                A = assembly.assemble_A(opts.N, opts.L, mu_eff, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                b = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                u_new = np.linalg.solve(A, b)
                u_corr = u_new - u
                u += opts.rel_factor*u_corr
                
            if opts.global_type=="flowrate":
                error_flowrate = 10
                while np.abs(error_flowrate) > 1e-3*flowrate:
                    #find corrected pressure gradient using calculated flowrate and assuming linear relation
                    correction_flowrate =  (np.sum(u)*opts.L/opts.N) / flowrate
                    p_grad /= correction_flowrate
                    b_new = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                    u = np.linalg.solve(A, b_new)
                    error_flowrate = np.sum(u)*opts.L/opts.N - flowrate 
        
    visuals.visualize(u*u_tau)

if __name__ == "__main__":
	main()
