import argparse

def get_options():
	"""
    Use this in the terminal like python3 main.py --density 3.0
    """
	parser = argparse.ArgumentParser(prog='Simulation of multiphase flow')
	parser.add_argument("--iterations", help='Number of iterations', type=int, default = 30)
	parser.add_argument("--p_grad", help='prescribed gradient in x direction [Pa]', type=float, default=-80)
	parser.add_argument("--rho_ref", help='reference density [kg/m3]', type=float, default=1000)
	parser.add_argument("--N", help='Amount of grid points', type=int, default=2000)
	parser.add_argument("--mu", help='Molecular viscosity [Pas]', type=float, default=0.001)
	parser.add_argument("--Ks", help='Sandgrain roughness [m]', type=float, default=0.001)
	parser.add_argument("--L", help="Width of the channel [m]", type=float, default=1)
	parser.add_argument("--flow_type", help="Either laminar or turbulent", type=str, default='turbulent')
	parser.add_argument("--rel_factor", help='Relaxation factor for updating turbulent velocity', type=float, default=0.5)
	parser.add_argument("--sim_time", help='Amount of time steps to simulate for', type=int, default=1)
	parser.add_argument("--dt", help='dt in particle tracking', type=float, default=0.0001)
	parser.add_argument("--tracktime", help='length of time to track particle for', type=float, default=0.02)
	parser.add_argument("--D", help='particle diameter as fraction of L [-]', type=float, default=0.01)
	parser.add_argument("--M", help='Particle mass', type=float, default=10e-7)
	parser.add_argument("--pressure_period", help='The period of the cosine pressure in time', type=float, default=20)
    
	opts = parser.parse_args()
	return opts
