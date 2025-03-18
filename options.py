import argparse

def get_options():
    """
    Use this in the terminal like python3 main.py --density 3.0
    """
    parser = argparse.ArgumentParser(prog='Simulation of multiphase flow')
    parser.add_argument("--iterations", help='Number of iterations', type=int, default = 30)
    parser.add_argument("--p_grad", help='prescribed gradient in x direction [Pa]', type=float, default=-20)
    parser.add_argument("--rho_ref", help='reference density [kg/m3]', type=float, default=1000)
    parser.add_argument("--N", help='Amount of grid points', type=int, default=100)
    parser.add_argument("--mu", help='Molecular viscosity [Pas]', type=float, default=0.001)
    parser.add_argument("--global_type", help='Type of global boundary condition', type=str, default='p_grad')
    parser.add_argument("--flowrate", help='Flow rate global boundary condition [kg/s]', type=float, default=20)
    parser.add_argument("--L", help="Width of the channel [m]", type=float, default=1)
    parser.add_argument("--bndry_bot", help="Boundary condition type on bottom wall", type=str, default='dirichlet')
    parser.add_argument("--bndry_top", help="Boundary condition type on top wall", type=str, default='dirichlet')
    parser.add_argument("--bndry_val_bot", help="Boundary value on bottom wall", type=float, default=0)
    parser.add_argument("--bndry_val_top", help="Boundary value on top wall", type=float, default=0)
    parser.add_argument("--flow_type", help="Either laminar or turbulent", type=str, default='turbulent')
    parser.add_argument("--rel_factor", help='Relaxation factor for updating turbulent velocity', type=float, default=0.5)
    
    opts = parser.parse_args()
    return opts
