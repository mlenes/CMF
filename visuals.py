import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

def visualize(u_list):
	u_list = np.array(u_list)
	fig, ax = plt.subplots()
	line = ax.plot(u_list[0], np.linspace(0,1,len(u_list[0])))[0]
	ax.set(xlim=[0, 1.1*np.max(u_list)], ylim=[0,1], xlabel="u [m/s]", ylabel="y/H [-]")
	
	def update(frame):
		line.set_xdata(u_list[frame])
		return line
		
	ani = animation.FuncAnimation(fig=fig, func=update, frames=len(u_list), interval=100)
	plt.show()

def particles(x_list,y_list):
	for i in range(len(x_list)):
		plt.plot(x_list[i],y_list[i], ".", label='particle')
	plt.xlabel("x [m]")
	plt.ylabel("y [m]")
	plt.legend()
	plt.show()
    
def wallprofile(L,N,n,conversion,uplus):
	channel=np.linspace(0,(n-1)/N,n)
	plt.figure()
	plt.plot(channel*conversion,uplus[1:n+1],".")	
	plt.xlabel('y+')
	plt.ylabel("u+")
	plt.xscale('log')
	plt.show()
