import argparse

def get_options():
    """
    Use this in the terminal like python3 main.py --density 3.0
    """
    parser = argparse.ArgumentParser(prog='Simulation of multiphase flow')
    parser.add_argument("--iterations", help='Number of iterations', type=int, default = 11)
    parser.add_argument("--p_grad", help='prescribed gradient in x direction', type=float, default=-0.8)
    parser.add_argument("--N", help='Amount of grid points', type=int, default=50)
    parser.add_argument("--mu", help='Molecular viscosity', type=float, default=0.001)
    parser.add_argument("--global_type", help='Type of global boundary condition', type=str, default='pressure')
    parser.add_argument("--flowrate", help='Flow rate global boundary condition', type=float, default=2)
    parser.add_argument("--L", help="Width of the channel", type=float, default=5)
    parser.add_argument("--bndry_bot", help="Boundary condition type on bottom wall", type=str, default='dirichlet')
    parser.add_argument("--bndry_top", help="Boundary condition type on top wall", type=str, default='dirichlet')
    parser.add_argument("--bndry_val_bot", help="Boundary value on bottom wall", type=float, default=0)
    parser.add_argument("--bndry_val_top", help="Boundary value on top wall", type=float, default=0)
    parser.add_argument("--flow_type", help="Either laminar or turbulent", type=str, default='turbulent')
    parser.add_argument("--rel_factor", help='Relaxation factor for updating turbulent velocity', type=float, default=0.5)
    
    opts = parser.parse_args()
    return opts
