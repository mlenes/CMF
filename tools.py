import numpy as np

def get_dy(N, L):
	return L/N

def init_mu_faces(N,mu):
	#create array of viscosity at gridpoints
	viscosity=np.ones(N+1)*mu
	mu_faces = (viscosity[:-1] + viscosity[1:])/2
	return mu_faces

def calc_mixing_length(N,L,u_tau,mu):
    #calculate mixing length at faces (non-dimensional)
    first_half = np.linspace(0,L/2,int(N/2)) / (L/2)
    rel_height = np.pad(first_half,(0,int(N/2)),'symmetric')
    mixing_length = 1/2 * (np.full(N,0.14) - np.full(N,0.08) * (np.ones(N)-rel_height)**2 - np.full(N,0.06) * (np.ones(N)-rel_height)**4)
    #van Driest damping
    y_plus = (rel_height * u_tau) / mu
    A_plus = 25
    D = (1-np.exp(-y_plus/A_plus))
    # D=1
    return mixing_length*D
    
def get_sin_p(p0, T, t):
	# Pressure difference over time with a sine
	return 0.5*p0*np.cos(2*np.pi*t/T) + 0.5*p0

def Gaussian(sigma,x):
    prob = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-x**2/(2*sigma**2))
    return prob
    
def choose_vprime(u_prime):
    loop = True
    while loop==True:
        if u_prime==0:
            return 0
            break
        random1 = 100*u_prime*(2*np.random.rand(1)-1)
        prob = Gaussian(u_prime,random1)
        random2 = 1/(u_prime*np.sqrt(2*np.pi))*np.random.rand(1)
        if prob>random2:
            return random1
            break
    
    
	
    

