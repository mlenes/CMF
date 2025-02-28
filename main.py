import numpy as np
from options import get_options
import assembly
import visuals

def main():
	# Include the options from the parser
    opts = get_options()
    
    global_vector = np.array([1,0]) #default is given pressure gradient
    if opts.global_type=="flowrate":
        global_vector = np.array([0,1])
    p_grad =  np.inner(global_vector, [opts.p_grad , 4*opts.flowrate/opts.L])

    if opts.flow_type == 'laminar':
        A = assembly.assemble_A(opts.N, opts.L, opts.mu, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
        b = assembly.assemble_b(opts.N, opts.L, opts.mu, p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
        u = np.linalg.solve(A, b)
        
        if opts.global_type=="flowrate":
            #flowrate and pressure are linearly related, so find corrected pressure gradient using calculated flowrate
            correction_flowrate =  (np.sum(u)*opts.L) / opts.flowrate
            corrected_p_grad = p_grad / correction_flowrate
            b_new = assembly.assemble_b(opts.N, opts.L, opts.mu, corrected_p_grad, opts.bndry_bot, opts.bndry_top, opts.bndry_val_bot, opts.bndry_val_top)
            u = np.linalg.solve(A, b_new)

    elif opts.flow_type == 'turbulent':
        pass

    else:
        print("Flow type not implemented")
        
    visuals.visualize(u)

if __name__ == "__main__":
	main()
