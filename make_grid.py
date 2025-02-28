# make 1D grid
import numpy as np

def init_mu(N,mu):
    #create array of viscosity at gridpoints
    viscosity=np.ones(N+1)*mu
    return viscosity


            
    
