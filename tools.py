import numpy as np

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces