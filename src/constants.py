#!/usr/bin/python3

########################################################################
# parameters might be changed
########################################################################
# number of cycles (divisions) between each grid points to 
# reach the convergence
loops = 30

# number of maximum terms of series to test the convergence
end_n = 70

# number of outward cycles to get outward solution for orbital
outward_cycles = 60




########################################################################
#  constants
########################################################################
e      = 2.718281828459045       # natural log constant
pi     = 3.141592653589793       # pi constant
c_au   = 137.0359991               # spped of light


# maximum atomic number of atoms to  be included
max_Z = 118
# some other constants
half = 1.0/2.0
third = 1.0/3.0

# Bound state limit
E_init  = 0.0           # intitial estimated eigenvalue could be any
						# value less /or eq. zero for bound states
						   
eta     = 1.0E-10           # default relative error


# conversion factors

Hartree_to_eV = 27.2113962      # conversion factor from Hartree to eV

########################################################################
# finite nuclear size models (potential options 8 and 9)
########################################################################
# Bohr radius in femtometres (CODATA 2018); 1 fm = 1/a0_fm atomic units
a0_fm = 52917.7210903

# 10%-90% surface thickness t of the Fermi charge distribution (fm).
# The diffuseness is a = t/(4 ln 3) = 0.5234 fm for t = 2.3 fm
fermi_skin_t_fm = 2.3

########################################################################
## print errors and results of loops
# to print messages with errors

print_error_messages = True
#print_error_messages = False

# print output of outward and inward loops
# print iteration results

print_iteration_results = True

########################################################################
# WKB parameter for R infintiy
# defalut value is 35.0

R_WKB_parameter = 35

########################################################################


