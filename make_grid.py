# make 1D grid
import numpy as np

def vertical_dist(N,dy,Nlayers,layerfactor):
    #create array with value of dy between all grid points
    #for refinement towards edges of the domain: Nlayers is amount of affected cells,
    #layerfactor is factor between adjecent cells
    #Nlayers should thus be at least twice as small as N
    vertical_distance=np.ones(N-1)*dy
    if Nlayers>0:
        for i in range(Nlayers-1):
            vertical_distance[i]=vertical_distance[i]/(layerfactor**(Nlayers-i))
            vertical_distance[-(i+1)]=vertical_distance[-(i+1)]/(layerfactor**(Nlayers-i))
    return vertical_distance


def init_mu(N,mu):
    #create array of viscosity at gridpoints
    viscosity=np.ones(N)*mu
    return viscosity


def init_u_dirichlet(N,flowrate,wallspeed,vertical_distance):
    #create array of velocities at gridpoints
    #Dirichelt boundary conditions dictate value of first and last cell, we set last boundary to 0
    u=flowrate/np.sum(vertical_distance[:-1])
    velocity=np.ones(N)*u
    velocity[0]=wallspeed+wallspeed-u #velocity is defined a center ghost cell??
    velocity[-1]=-u
    return velocity

def init_u_neumann(N,flowrate,wallgradient,vertical_distance):
    #create array of initial velocities at gridpoints
    #Neumann boundary conditions dictate gradient of velocity at first and last cell
    u=flowrate/np.sum(vertical_distance[:-1])
    velocity=np.ones(N)*u  
    velocity[0]=velocity[1]-vertical_distance[0]*wallgradient  #velocity at center ghost cell??
    velocity[-1]=velocity[-2]+vertical_distance[-1]*wallgradient
    return velocity


            
    
