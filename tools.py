import numpy as np
import matplotlib.pyplot as plt

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces

def calc_mixing_length(N,L,u_tau,mu):
    #calculate mixing length at faces (non-dimensional)
    first_half = np.linspace(0,L/2,int(N/2)) / L
    rel_height = np.pad(first_half,(0,int(N/2)),'symmetric')
    mixing_length = L/2 * (np.full(N,0.14) - np.full(N,0.08) * (np.ones(N)-2*rel_height)**2 - np.full(N,0.06) * (np.ones(N)-2*rel_height)**4)
    
    #van Driest damping
    y_plus = (rel_height * u_tau) / mu
    A_plus = 25
    D = (1-np.exp(-y_plus/A_plus))
    
    return mixing_length*D/L
    
def get_sin_p(p0, T, t):
	# Pressure difference over time with a sine
	return 0.5*p0*np.cos(2*np.pi*t/T) + 0.5*p0
    

