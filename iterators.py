import numpy as np
import assembly
import tools

	
def iter_u(iterations, u, N, L, Ks, mu_faces, rel_factor, p_grad, wall_constant=None):
    for i in range(iterations):
        #each iteration first loop until velocity profile matches eddy viscosity
        u_corr = 10*np.mean(u)

        while np.abs(np.mean(u_corr)) > np.abs(1e-4*np.mean(u)):
            tau_wall = mu_faces[0]*(u[1]-u[0])/(tools.get_dy(N, L)/L)
            u_tau = np.sqrt(tau_wall)
            
            # If you use a wall function, use log function for tau_wall
            if wall_constant!=None:
                u_tau = u[1] * 0.41 / (np.log((tools.get_dy(N, L)*32.6/2)/Ks))            
                tau_w = u_tau**2 #non-dimensional, so rho=1
                wall_constant = tau_w/u[1]
            
            # Eddy viscosity
            l_m = tools.calc_mixing_length(N, L, u_tau, mu_faces[1])
            u_prime = l_m * np.abs(u[:-1] - u[1:])/(tools.get_dy(N, L)/L)
            mu_eff = mu_faces + u_prime*l_m
            
            # Get turbulent A and b
            A = assembly.assemble_A(N, L, mu_eff, wall_constant)
            b = assembly.assemble_b(N, L, mu_eff, p_grad)
            u_new = np.linalg.solve(A, b)
            u_corr = u_new - u
            u += rel_factor*u_corr      
    return u, u_prime, u_tau

def particle(y0, v0, u, dt, tracktime, mu, D, M, g, N, L, u_prime):
	x_list = []
	y_list = []
	sim_steps = int(tracktime/dt)
	vprime = np.zeros(len(y0))
	eddytime = np.zeros(len(y0))
# 	eddytime = np.ones(len(y0))*tracktime
	f=0.01 #friction factor, heb ik random gekozen nu
	
	#loop for every particle we inserted
	for n in range(len(y0)):
		print("start particle",n)
		x = np.zeros(sim_steps)
		y = np.zeros(sim_steps)
		y[0] = y0[n]
		vx = np.zeros(sim_steps)
		vx[0] = v0[0]
		vy = np.zeros(sim_steps)
		vy[0] = v0[1]
		
	
		#update position and velocity by calculating force
		for i in range(1,len(x)):
			x[i] = x[i-1] + vx[i-1]*dt
			y[i] = y[i-1] + vy[i-1]*dt
			gridcell = int(y[i-1]*N)%N
			
			
			#if eddy died or particle passed through eddy, make a new one
			if dt*i>eddytime[n]:
				dUdy = np.abs(u[gridcell] - u[gridcell+1])/(tools.get_dy(N, L)/L)
				vprime[n] = tools.choose_vprime(np.sqrt(u_prime[gridcell]*dUdy))
				
				T_e = 0.15/dUdy #eddy lifetime
				T_r = (np.abs(vprime[n])/dUdy) / np.max([np.abs(vy[i-1]),np.abs(u[gridcell]-vx[i-1])]) #residence time
				eddytime[n]+= np.min([T_e,T_r]) #check in how much time this eddy will last				
				print(vprime[n],eddytime[n])

			Fx = f*3*np.pi*mu*D*(u[gridcell]+vprime[n]-vx[i-1]) #now only drag force, all input is non-dimensional
			
			vx[i] = vx[i-1] + 1/(M+0.5*np.pi*D**3/6) * Fx * dt	#added mass and non-dimensionally rho_fluid=1
				
			Fy = (np.pi*D**3/6-M) * g + f*3*np.pi*mu*D*(vprime[n]-vy[i-1]) #stokes drag and gravity
		 
			vy[i] = vy[i-1] + 1/(M+0.5*np.pi*D**3/6) * Fy * dt #added mass
			
			if y[i]<0.5*D or y[i]>(1-0.5*D):
				vy[i] = -vy[i-1]
		x_list.append(x)
		y_list.append(y)
		
	return np.array(x_list), np.array(y_list) #this is the position of all particles over time
