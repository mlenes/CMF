import matplotlib.pyplot as plt
import numpy as np

def visualize(u):
    plt.plot(u, np.linspace(0,1,len(u)))
    plt.xlabel("u [m/s]")
    plt.ylabel("y/H [-]")
    plt.ylim(0,1)
    plt.xlim(0)
    plt.show()
