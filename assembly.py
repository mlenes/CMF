import numpy as np
import tools

def assemble_A(N, L, mu_faces, wall_function, wall_constant):
	A = np.zeros((N+1, N+1),dtype=float) # Coefficient matrix
	dy = tools.get_dy(N, L)/L
 	
	for i in range(1, N):
		A[i, i-1] = mu_faces[i-1]/dy**2
		A[i , i+1] = mu_faces[i]/dy**2
		A[i,i] = -(mu_faces[i-1] + mu_faces[i])/dy**2
    	
	#bottom boundary condition
	A[0,0] = 1
	A[0,1] = 1
	if wall_function == True:
		A[1,0] = 0
		A[1,1] = -mu_faces[1]/dy**2 - wall_constant/dy
    
	#top boundary condition   
	A[-1,-1] = 1
	A[-1,-2] = 1
	if wall_function == True:
		A[-2,-1] = 0
		A[-2,-2] = -mu_faces[-2]/dy**2 - wall_constant/dy
 	    
	return A

def assemble_b(N, L, mu_faces, p_grad):
	b = np.full(N+1, p_grad, dtype=float) # Right hand side
	
	#bot boundary condition
	b[0] = 0
		
	#top boundary condition
	b[-1] = 0

	return b
    


