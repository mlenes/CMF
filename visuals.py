import matplotlib.pyplot as plt
import numpy as np

def visualize(u):
	plt.plot(u, range(len(u)))
	plt.xlabel("u [m/s]")
	plt.show()