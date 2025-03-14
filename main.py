import numpy as np
from options import get_options
import assembly
import visuals
import tools

def main():
	# Include the options from the parser
    opts = get_options()
    
    global_vector = np.array([1,0]) #default is given pressure gradient
    if opts.global_type=="flowrate":
        global_vector = np.array([0,1])
    p_grad =  np.inner(global_vector, [opts.p_grad , -4*opts.flowrate/opts.L*opts.mu])
    # Compute viscosity at the faces of the cells
    mu_faces = tools.init_mu_faces(opts.N, opts.mu)
    
    A = assembly.assemble_A(opts.N, opts.L, mu_faces, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
    b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
    u = np.linalg.solve(A, b)
    
    if opts.global_type=="flowrate":
        error_flowrate = 10
        while np.abs(error_flowrate) > 1e-3*opts.flowrate:
            #find corrected pressure gradient using calculated flowrate and assuming linear relation
            correction_flowrate =  (np.sum(u)*opts.L) / opts.flowrate
            p_grad /= correction_flowrate
            b_new = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
            u = np.linalg.solve(A, b_new)
            error_flowrate = np.sum(u)*opts.L - opts.flowrate 
            
            

    if opts.flow_type == 'turbulent':
        for i in range(opts.iterations):
            #each iteration first loop until velocity profile matches eddy viscosity
            #if flowrate is given, then check converged velocity profile with flowrate and if needed correct
            u_corr = 10
            while np.abs(np.mean(u_corr)) > 1e-5*np.mean(u):
                #assume rho=1
                mu_eff = mu_faces + tools.calc_mixing_length(opts.N, opts.L, opts.bndry_bot, opts.bndry_top)**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(opts.N,opts.L)
                A = assembly.assemble_A(opts.N, opts.L, mu_eff, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                b = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                u_new = np.linalg.solve(A, b)
                u_corr = u_new - u
                u += opts.rel_factor*u_corr
                
            if opts.global_type=="flowrate":
                error_flowrate = 10
                while np.abs(error_flowrate) > 1e-3*opts.flowrate:
                    #find corrected pressure gradient using calculated flowrate and assuming linear relation
                    correction_flowrate =  (np.sum(u)*opts.L) / opts.flowrate
                    p_grad /= correction_flowrate
                    b_new = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
                    u = np.linalg.solve(A, b_new)
                    error_flowrate = np.sum(u)*opts.L - opts.flowrate 
                    print(error_flowrate)


    else:
        print("Flow type not implemented")
        
    visuals.visualize(u)

    tools.calc_mixing_length(opts.N,opts.L, opts.bndry_bot, opts.bndry_top)

if __name__ == "__main__":
	main()
