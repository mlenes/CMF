import numpy as np

def simulation(init_u, init_mu, dy_arr, iterations, p_grad):
    """
    Arguments
    init_u : np.ndarray
        initial velocity in the x direction for iteration n
    init_mu : np.ndarray
        initial viscosity for iteration n
    dy_arr : np.ndarray
        difference in y position between two grid points
    iterations : int
        amount of iterations to perform
    p_grad : float
        Prescribed gradient in the x direction
    """
    # Number of gridpoints to simulate (+2 ghost points)
    N = init_u.shape[0]
    
    # Initialize the velocity profiles
    u = np.zeros((interations, N))
    u[0,:] = init_u
    
    mu = np.zeros((interactions, N))
    mu[0,:] = init_mu
    
    for n in iterations:
        for i in range(1,N):
            u[n,i] = (1/2)*(dy_arr[i]**2)*p_grad/mu[n,i] - (1/2)*u[n,i+1] - (1/2)*u[n,i-1]
    return u
            
        
