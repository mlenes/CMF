import numpy as np
import matplotlib.pyplot as plt

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces

def calc_mixing_length(N,L, bndry_bot, bndry_top):
    #calculate mixing length at faces
    rel_height = np.linspace(0,L,N) / L
    if bndry_bot == "dirichlet" and bndry_top == "dirichlet":
        mixing_length = L/2 * (np.full(N,0.14) - np.full(N,0.08) * (np.ones(N)-2*rel_height)**2 - np.full(N,0.06) * (np.ones(N)-2*rel_height)**4)
    else:
        mixing_length = L * (np.full(N,0.14) - np.full(N,0.08) * (np.ones(N)-rel_height)**2 - np.full(N,0.06) * (np.ones(N)-rel_height)**4)
    return mixing_length
    

