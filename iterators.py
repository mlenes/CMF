import numpy as np
import assembly

def iter_flowrate(N, L, mu_faces, p_grad, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top, A, u, flowrate):

	# Dummy initial value
	error_flowrate = 10
	
	while (np.abs(error_flowrate) > 1e-3*flowrate):
		# Find corrected pressure gradient using calculated flowrate and assuming linear relation
		correction_flowrate =  (np.sum(u)*L/N) / flowrate
		p_grad /= correction_flowrate
		
		b = assembly.assemble_b(N, L, mu_faces, p_grad, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top)
		u = np.linalg.solve(A, b)
		
		error_flowrate = np.sum(u)*L/N - flowrate
		
	return u, p_grad
