import numpy as np

def assemble_matrix(N, L, mu, p_grad, u_wall):
	dy = L/(N-1)
	A = np.zeros((N-2, N-2)) # Coefficient matrix
	b = np.full(N-2, -p_grad) # Right hand side

	# Compute viscosity at the faces of the cells
	mu_faces = (mu[:-1] + mu[1:])/2
	
	for i in range(N-2):
		if i > 0:
			A[i, i-1] = mu_faces[i]/dy**2
		if i < N-3:
			A[i , i+1] = mu_faces[i+1]/dy**2
		A[i,i] = -mu_faces[i]/dy**2 + -mu_faces[i+1]/dy**2

	return A, b
