import numpy as np

def assemble_A(N, L, mu, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top):
	
	A = np.zeros((N+1, N+1),dtype=float) # Coefficient matrix

	dy = get_dy(N, L)
	# Compute viscosity at the faces of the cells
	mu_faces = init_mu_faces(N, mu)
	
	for i in range(1, N):
		A[i, i-1] = mu_faces[i-1]/dy**2
		A[i , i+1] = mu_faces[i]/dy**2
		A[i,i] = -(mu_faces[i-1] + mu_faces[i])/dy**2
	
	c1_bot, c1_top = get_c1_coeff(bndry_bot, bndry_top, bndry_val_bot, bndry_val_top)
	
	#bottom boundary condition
	A[0,0] = 1
	A[0,1] = c1_bot
	A[-1,-1] = 1
	A[-1,-2] = c1_top
	
	return A

def assemble_b(N, L, mu, p_grad, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top):
	b = np.full(N+1, p_grad, dtype=float) # Right hand side
	
	dy = get_dy(N, L)
	# Compute viscosity at the faces of the cells
	mu_faces = init_mu_faces(N, mu)
	
	b[0], b[-1] = get_c2_coeff(bndry_bot, bndry_top, bndry_val_bot, bndry_val_top, dy, mu_faces)
	return b

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces
    
def get_c1_coeff(bndry_bot, bndry_top, bndry_val_bot, bndry_val_top):
	if bndry_bot == 'dirichlet':
		c1_bot = 1
	elif bndry_bot in ['neumann', 'stress']:
		c1_bot = -1
	else:
		print("Boundary type not implemented")
		
	if bndry_top == 'dirichlet':
		c1_top = 1
	elif bndry_top in ['neumann', 'stress']:
		c1_top = -1
	else:
		print("Boundary type not implemented")

	return c1_bot, c1_top
    
def get_c2_coeff(bndry_bot, bndry_top, bndry_val_bot, bndry_val_top, dy, mu_faces):
	if bndry_bot == 'dirichlet':
		c2_bot = 2*bndry_val_bot
	elif bndry_bot == 'neumann':
		c2_bot = -dy*bndry_val_bot
	elif bndry_bot == 'stress':
		c2_bot = -dy*bndry_val_bot/(mu_faces[0])
	else:
		print("Boundary type not implemented")
		
	if bndry_top == 'dirichlet':
		c2_top = 2*bndry_val_top
	elif bndry_top == 'neumann':
		c2_top = dy*bndry_val_top
	elif bndry_top == 'stress':
		c2_top = dy*bndry_val_top/(mu_faces[-1])
	else:
		print("Boundary type not implemented")

	return c2_bot, c2_top
