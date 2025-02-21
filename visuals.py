import matplotlib.pyplot as plt
import numpy as np


def plot_profile(vertical_distance, u):
	y_pos = np.cumsum(vertical_distance)
	plt.plot(u[1:-1], y_pos[:-1], '.')
	plt.xlabel('Velocity')
	plt.ylabel('y-position')
	
	plt.show()