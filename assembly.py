import numpy as np
import tools

def assemble_A(N, L, mu_faces, tau_w=None):
	A = np.zeros((N+1, N+1),dtype=float) # Coefficient matrix
	dy = tools.get_dy(N, L)/L
 	
	for i in range(1, N):
		A[i, i-1] = mu_faces[i-1]/dy**2
		A[i , i+1] = mu_faces[i]/dy**2
		A[i,i] = -(mu_faces[i-1] + mu_faces[i])/dy**2
    	
	#bottom boundary condition
	A[0,0] = 1
	A[0,1] = 1
	if tau_w != None:
		A[1,0] = 1
		A[1,1] = -1
    
	#top boundary condition   
	A[-1,-1] = 1
	A[-1,-2] = 1
	if tau_w != None:
		A[-2,-1] = 1
		A[-2,-2] = -1
 	    
	return A

def assemble_b(N, L, mu_faces, p_grad, tau_w=None):
	dy = tools.get_dy(N, L)/L
	b = np.full(N+1, p_grad, dtype=float) # Right hand side
	
	#bot boundary condition
	b[0] = 0

	#top boundary condition
	b[-1] = 0

	if tau_w != None:
		# Enforcement of wall boundary condition
		b[1] = dy*tau_w/mu_faces[0]
		b[N-1] = dy*tau_w/mu_faces[N-1]

	return b
    


