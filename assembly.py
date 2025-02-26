import numpy as np

def assemble_matrix(N, L, mu, p_grad, bndry_type, u_wall, grad_wall):
	dy = L/(N-1)
	A = np.zeros((N-2, N-2),dtype=float) # Coefficient matrix
	b = np.full(N-2, p_grad, dtype=float) # Right hand side

	# Compute viscosity at the faces of the cells
	mu_faces = (mu[:-1] + mu[1:])/2
	
	for i in range(1, N-2):
		if i > 0:
			A[i, i-1] = mu_faces[i]/dy**2
		if i < N-3:
			A[i , i+1] = mu_faces[i+1]/dy**2
		A[i,i] = -(mu_faces[i] + mu_faces[i+1])/dy**2

	if bndry_type == 'dirichlet':
		A[0, :] = 0
		A[0,0] = 1

		A[-1, :] = 0
		A[-1,-1] = 1

		b[0] = u_wall
		b[-1] = u_wall

	if bndry_type == 'neumann':
		A[0, :] = 0
		A[0,0] = 1
		A[0,1] = -1

		A[-1, :] = 0
		A[-1,-1] = 1
		A[-1,-2] = -1

		b[0] = grad_wall*dy
		b[-1] = grad_wall*dy

		A[N//2, :] = 0
		A[N//2, N//2] = 1  # Set velocity reference at middle
		b[N//2] = 0  # Reference velocity set to zero (could be any value)

	return A, b
