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


def init_u_dirichlet(N,u,Dirichlet):
    #create array of velocities at gridpoints
    #Dirichelt boundary conditions dictate value of first and last cell
    velocity=np.ones(N)*u
    velocity[0]=Dirichlet
    velocity[-1]=0
    return velocity

def init_u_neumann(N,u,Neumann,vertical_distance):
    #create array of initial velocities at gridpoints
    #Neumann boundary conditions dictate gradient of velocity at first and last cell
    velocity=np.ones(N)*u  
    velocity[0]=velocity[1]-vertical_distance[0]*Neumann


            
    
