import matplotlib.pyplot as plt
import numpy as np

def visualize(u):
    plt.plot(u[1:-1], np.linspace(0,1,len(u)-2))
    plt.xlabel("u [m/s]")
    plt.ylabel("y/H [-]")
    plt.show()
