#!/usr/bin/python3

from constants import e

# GRASP grid
class GRID:
	"""
	The GRID class defines the used grid for our numerical calculations
	#
	r_i = r_0 [ exp(t_i) - 1 ]
	t_i = ( i - 1 ) h
	#
	i = 1,2,3, .....,N
	# 
	and r_0 = T/Z
	#
	where 
	N = number of grid points  (default = 590)
	T = initial start point  (default = 1.1E-5 a.u.)
	h = step difference (default = 5.000E-2)
	#
	GRASP similariies
	T   == RNT
	h 	== H
	N 	== N
	"""
	def __init__(self, N=590, T=2.0E-5, H=5.0E-2, atomic_number=50):
		self.n 		= N
		self.T  	= T
		self.h 		= H
		self.Zz      = atomic_number
		
	def t(self):
		t = []
		for m in range(self.n):
			t1 = m * self.h
			t.append(t1)
		return t
	#
	def grid(self):
		t = self.t()
		r0 = float(self.T) / self.Zz
		r = []
		for k in range(self.n):
			r1 = r0 * ( (e**(t[k])) - 1 )
			r.append(r1)
		return r



#----------------#
class GRID_90:
	"""
	grid form represented as 
	r0 = 0.0
	r1 = 1.0E-8
	d = 1.06E0
	
	r_n = r_{n-1} * d
	"""
	def __init__(self,NP=3000,r_1=1.0E-8,D=1.012):
		self.np = NP
		self.r1 = r_1
		self.d = D
	#
	def grid_90(self):
		out = [0.0]
		out.append(self.r1)
		for i in range(2,self.np):
			rn = out[i-1]*self.d
			out.append(rn)
		return out
#----------------#
class GRID_C:
	"""
	NOT Included
	grid form represented as 
	NP = 500
	r0 = 1.0E-8
	h = 1.06E0
	
	r_n = r_{0} * h (e^{n} - 1)
	"""
	def __init__(self,NP=800,r_0=1.0E-4,h=1.06E0):
		self.np = NP
		self.r0 = r_0
		self.h = h
	#
	def grid_c(self):
		out = []
		for i in range(self.np):
			rn = self.r0*self.h*((e**i - 1 )**((i)*self.h) )
			out.append(rn)
		return out
