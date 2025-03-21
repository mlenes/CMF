import numpy as np
import tools

def assemble_A(N, L, mu_faces, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top, wall_function,wall_constant):
    	
   	A = np.zeros((N+1, N+1),dtype=float) # Coefficient matrix
   
   	dy = tools.get_dy(N, L)
   	
   	for i in range(1, N):
   		A[i, i-1] = mu_faces[i-1]/dy**2
   		A[i , i+1] = mu_faces[i]/dy**2
   		A[i,i] = -(mu_faces[i-1] + mu_faces[i])/dy**2
   	
   	#bottom boundary condition
   	A[0,0] = 1
   	
    if bndry_bot == 'dirichlet':
        A[0,1] = 1
        if wall_function == True:
            A[1,0] = 0
            A[1,1] = -mu_faces[1]/dy**2 - wall_constant/dy
   	elif bndry_bot in ['neumann', 'stress']:
   		A[0,1] = -1
   	else:
   		print("Bot boundary condition not implemented")
   
   	#top boundary condition
    A[1,-1] = 1
       
    if bndry_top == 'dirichlet':
        A[-1,-2] = 1
        if wall_function == True:
            A[-2,-1] = 0
            A[-2,-2] = -mu_faces[-2]/dy**2 - wall_constant/dy
    elif bndry_top in ['neumann', 'stress']:
    	A[-1,-2] = -1
    else:
        print("Top boundary condition not implemented")
    
    if wall_function == True:
        A[1,0] = 0
        A[1,1] = -mu_faces[1]/dy**2 - wall_constant/dy
    
    	return A

def assemble_b(N, L, mu_faces, p_grad, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top):
	b = np.full(N+1, p_grad, dtype=float) # Right hand side
	
	dy = tools.get_dy(N, L)
	
	#bot boundary condition
	if bndry_bot == 'dirichlet':
		b[0] = 2*bndry_val_bot
	elif bndry_bot == 'neumann':
		b[0] = -dy*bndry_val_bot
	elif bndry_bot == 'stress':
		b[0] = -dy*bndry_val_bot/(mu_faces[0])
	else:
		print("Boundary type not implemented")
		
	#top boundary condition
	if bndry_top == 'dirichlet':
		b[-1] = 2*bndry_val_top
	elif bndry_top == 'neumann':
		b[-1] = -dy*bndry_val_top
	elif bndry_top == 'stress':
		b[-1] = -dy*bndry_val_top/(mu_faces[0])
	else:
		print("Boundary type not implemented")
	
	return b
    


