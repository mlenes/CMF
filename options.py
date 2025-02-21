import argparse

def get_options():
	"""
	Use this in the terminal like python3 main.py --density 3.0
	"""
	parser = argparse.ArgumentParser(prog='Simulation of multiphase flow')
	parser.add_argument("--density", help='Density of the fluid', type=float)
	parser.add_argument("--iterations", help='Number of iterations', type=int, default = 100)
	parser.add_argument("--p_grad", help='prescribed gradient in x direction', type=float)
	parser.add_argument("--Nlayers", help='How many layers for grid', type=int, default=0)
	parser.add_argument("--layerfactor", help='How refined the grid is at the boundary', type=float, default=2)
	parser.add_argument("--dy", help='standard grid size', type=float, default=0.01)
	parser.add_argument("--N", help='Amount of grid points', type=int, default=50)
	parser.add_argument("--mu", help='Molecular viscosity', type=float, default=0.001)
	parser.add_argument("--wallspeed", help='speed of wall', type=int, default=0)
	parser.add_argument("--wallgradient", help='gradient of velocity at the wall', type=float, default=0)
	parser.add_argument("--flowrate", help='Flow rate global boundary condition', type=float, default=2)
	opts = parser.parse_args()
	return opts
