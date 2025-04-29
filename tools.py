import numpy as np
import matplotlib.pyplot as plt

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces

def calc_mixing_length(N,L):
    #calculate mixing length at faces (non-dimensional)
    rel_height = np.linspace(0,L,N) / L
    mixing_length = L/2 * (np.full(N,0.14) - np.full(N,0.08) * (np.ones(N)-2*rel_height)**2 - np.full(N,0.06) * (np.ones(N)-2*rel_height)**4)
    return mixing_length/L
    
def get_sin_p(p0, T, t):
	# Pressure difference over time with a sine
	return A*np.sin(2*np.pi*t/T)
    

