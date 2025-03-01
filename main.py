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
    
    if opts.flow_type == 'laminar':
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
            
            

    elif opts.flow_type == 'turbulent':
        A = assembly.assemble_A(opts.N, opts.L, mu_faces, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
        b = assembly.assemble_b(opts.N, opts.L, mu_faces, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
        u = np.linalg.solve(A, b)
        for i in range(opts.iterations):
            #assume rho=1
            mu_eff = mu_faces + tools.calc_mixing_length(opts.N, opts.L, opts.bndry_bot, opts.bndry_top)**2 * np.abs(u[:-1] - u[1:])/tools.get_dy(opts.N,opts.L)
            A = assembly.assemble_A(opts.N, opts.L, mu_eff, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
            b = assembly.assemble_b(opts.N, opts.L, mu_eff, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
            u = np.linalg.solve(A, b)
            # visuals.visualize(u)


    else:
        print("Flow type not implemented")
        
    visuals.visualize(u)

    tools.calc_mixing_length(opts.N,opts.L, opts.bndry_bot, opts.bndry_top)

if __name__ == "__main__":
	main()
