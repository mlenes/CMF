import matplotlib.pyplot as plt
import numpy as np

def lam_profile(y, p_grad, mu, L):
    return (1/(2*mu)) * p_grad * (y**2 - L*y)

def visualize(u, p_grad, mu, L):
    y = np.linspace(0, L, 100)/L
    print(lam_profile(y[50], p_grad, mu, L))
    plt.plot(lam_profile(y, p_grad, mu, L), y, label="Analytical Solution")
    plt.plot(u, np.linspace(0,1,len(u)), label="Numerical Solution")
    plt.legend()
    plt.xlabel("u [m/s]")
    plt.ylabel("y/H [-]")
    plt.ylim(0,1)
    plt.xlim(0)
    plt.savefig("u_profile.png")
