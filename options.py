import argparse

def get_options():
    """
    Use this in the terminal like python3 main.py --density 3.0
    """
    parser = argparse.ArgumentParser(prog='Simulation of multiphase flow')
    parser.add_argument("--iterations", help='Number of iterations', type=int, default = 30)
    parser.add_argument("--p_grad", help='prescribed gradient in x direction [Pa]', type=float, default=-0.8)
    parser.add_argument("--rho_ref", help='reference density [kg/m3]', type=float, default=1000)
    parser.add_argument("--N", help='Amount of grid points', type=int, default=100)
    parser.add_argument("--mu", help='Molecular viscosity [Pas]', type=float, default=0.001)
    parser.add_argument("--Ks", help='Sandgrain roughness [m]', type=float, default=0.001)
    parser.add_argument("--global_type", help='Type of global boundary condition', type=str, default='pressure')
    parser.add_argument("--flowrate", help='Flow rate global boundary condition [kg/s]', type=float, default=1000)
    parser.add_argument("--L", help="Width of the channel [m]", type=float, default=1)
    parser.add_argument("--flow_type", help="Either laminar or turbulent", type=str, default='turbulent')
    parser.add_argument("--rel_factor", help='Relaxation factor for updating turbulent velocity', type=float, default=0.5)
    parser.add_argument("--sim_time", help='Amount of time steps to simulate for', type=int, default=10)
    parser.add_argument("--pressure_period", help='The period of the cosine pressure in time', type=float, default=20)
    
    opts = parser.parse_args()
    return opts
