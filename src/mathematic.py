#!/usr/bin/python3


class MATH:
	"""
	solving some math functions:
	solving Tri Diagonal Matrix 
	see:
	https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	date: 29 June 2026
	# matrix form
	-                    -
	| b1 c1 0  0 ...... 0| |  x1 |    |  d1 |
	| a2 b2 c2 0 ...... 0| |  x2 |    |  d2 |
	| 0  a3 b3 c3 ..... 0| | . . |  = |  .. |
	|                    | | . . |    |  .. |
	|               c_n-1| |x_n-1|    |d_n-1|
	|0           a_n  b_n| | x_n |    |d_n  |
	-                    -
	input:
	 a = lower diagonal elements
	 b = diagonal elements
	 c = upper diagonal elements
	 d = is the free term
	output:
	 the x values in the matrix form
	"""
	def __init__(self, A, B, C, D):
		self.a 		= A
		self.b 	    = B
		self.c 		= C
		self.d      = D
	#
	def tri_diagonal_solver(self):
		# Check if there is data to process
		if len(self.d) == 0:
			return []
		# Check if there is only one equation
		if len(self.d) == 1:
			if self.b[0] != 0:
				return [self.d[0] / self.b[0]]
			else:
				return [0.0]
		#
		bB = list(self.b)
		dD = list(self.d)
		
		# Forward sweep
		for i in range(1, len(self.d)):
			if bB[i-1] == 0:
				# Avoid division by zero
				m = 0.0
			else:
				m = self.a[i-1] / bB[i-1]
			bB[i] = bB[i] - (m * self.c[i-1])
			dD[i] = dD[i] - (m * dD[i-1])
		
		# Check if the last diagonal is zero
		if bB[-1] == 0:
			return [0.0] * len(self.d)
		
		x = [0.0] * len(self.d)
		x[-1] = dD[-1] / bB[-1]
    
		# Back substitution
		for n in range(len(self.d)-2, -1, -1):
			if bB[n] == 0:
				x[n] = 0.0
			else:
				x[n] = (dD[n] - (self.c[n] * x[n+1])) / bB[n]
		
		return x
