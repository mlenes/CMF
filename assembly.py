import numpy as np
import tools

def assemble_A(N, L, mu_faces, bndry_bot, bndry_top, bndry_val_bot, bndry_val_top, wall_function, wall_constant):
    A = np.zeros((N + 1, N + 1), dtype=float)  # Coefficient matrix
    dy = tools.get_dy(N, L)

    # Interior nodes (central difference)
    for i in range(1, N):
        A[i, i - 1] = mu_faces[i - 1] / dy**2
        A[i, i + 1] = mu_faces[i] / dy**2
        A[i, i] = -(mu_faces[i - 1] + mu_faces[i]) / dy**2

    # Bottom boundary (i = 0)
    if bndry_bot == 'dirichlet':
        A[0, 0] = 1.0  # Strongly impose Dirichlet
    elif bndry_bot in ['neumann', 'stress']:
        A[0, 0] = -1.0 / dy
        A[0, 1] = 1.0 / dy
    else:
        raise NotImplementedError("Bottom boundary condition not implemented.")

    # Apply wall function if needed at i = 1
    if wall_function and bndry_bot == 'dirichlet':
        A[1, :] = 0  # Clear row
        A[1, 0] = -wall_constant / dy
        A[1, 1] = mu_faces[1] / dy**2 + wall_constant / dy
        A[1, 2] = -mu_faces[1] / dy**2

    # Top boundary (i = N)
    if bndry_top == 'dirichlet':
        A[N, N] = 1.0  # Strongly impose Dirichlet
    elif bndry_top in ['neumann', 'stress']:
        A[N, N] = 1.0 / dy
        A[N, N - 1] = -1.0 / dy
    else:
        raise NotImplementedError("Top boundary condition not implemented.")

    # Apply wall function if needed at i = N-1
    if wall_function and bndry_top == 'dirichlet':
        A[N - 1, :] = 0  # Clear row
        A[N - 1, N] = -wall_constant / dy
        A[N - 1, N - 1] = mu_faces[N - 1] / dy**2 + wall_constant / dy
        A[N - 1, N - 2] = -mu_faces[N - 1] / dy**2

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
    


