# make 1D grid
import numpy as np

def vertical_dist(N,dy,Nlayers,layerfactor):
    vertical_distance=np.ones(N-1)*dy
    if Nlayers>0:
        for i in range(Nlayers-1):
            vertical_distance[i]=vertical_distance[i]/(layerfactor**(Nlayers-i))
            vertical_distance[-(i+1)]=vertical_distance[-(i+1)]/(layerfactor**(Nlayers-i))
    return vertical_distance


            
    
