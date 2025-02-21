import numpy as np

def simulation_pressure(init_u, init_mu, dy_arr, iterations, p_grad):
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
    u = np.zeros((iterations, N))
    u[0,:] = init_u
    
    mu = np.zeros((iterations, N))
    mu[0,:] = init_mu
    
    for n in range(iterations):
        for i in range(1,N):
            u[n,i] = (1/2)*(dy_arr[i]**2)*p_grad/mu[n,i] - (1/2)*u[n,i+1] - (1/2)*u[n,i-1]
            #hier moeten we dan u in de ghost cells updaten, maar hoe is afhankelijk van 
            #het type boundary condition, voelt stom om hier weer allemaal if statements te maken?
    return u


def simulation_flowrate(init_u, init_mu, dy_arr, iterations, flowrate):
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
    u = np.zeros([iterations, N])
    u[0,:] = init_u
    
    mu = np.zeros([iterations, N])
    mu[0,:] = init_mu
    
    p_grad = np.zeros(iterations)
    
    #ik heb geprobeerd de flow rate boundary condition toe te passen, ik denk dat je dan de p_grad niet zomaar kan berekenen?
    for n in range(iterations):
        for i in range(1,N):
            p_grad [n+1] = (u[n,i] + (1/2)*u[n,i+1] + (1/2)*u[n,i-1]) *2*mu[n,i]/(dy_arr[i]**2)
            #ik denk dat de p_grad die hieruit komt nog niet klopt?
            u[n+1,i] = (1/2)*(dy_arr[i]**2)*p_grad [n+1,i]/mu[n] - (1/2)*u[n,i+1] - (1/2)*u[n,i-1]
            #door de foute p_grad telt de snelheid niet op tot de flowrate, dus dat kunnen we testen en fixen?
            Q =  np.sum(u[n+1,:-1]*dy_arr)  #flow rate berekenen, is shady met snelheden op grenzen, welke lengte neem je?
            errorQ = Q - flowrate
            #hopelijk een goede correctie maken op de snelheid, dan klopt het weer met de flow rate
            #volgende ronde is de p_grad dan hopelijk ook beter?
            u[n+1,i]+= np.ones(N)*errorQ/N
            #hier moeten we dan u in de ghost cells updaten, maar hoe is afhankelijk van 
            #het type boundary condition, voelt stom om hier weer allemaal if statements te maken?
    return u
            
        
