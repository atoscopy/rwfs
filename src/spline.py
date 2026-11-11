#!/usr/bin/python3

from mathematic import MATH

class SPLINE:
	"""
	This class function calculates the coefficients of cubic spline
	interpolation of 
	 for each interval x_i , x_{i+1}
	 Fun = A + B*x + C*x^2 + D*x^3     
	 (absolute x, valid within each interval [x_i, x_{i+1}])
	#
	ref. 
	1- https://doi.org/10.1016/j.cpc.2019.02.011
	  paper named: radial A Fortran subroutine package for the 
	  solution of the radial Schrodinger and Dirac wave equations
	  appendex A
	 2- Numerical Methods in Engineering with Python
	    Book by Jaan Kiusalaas page 115
	#
	input:
		x = x data points
		y = y data points
		k0 = second derivatives at x0
		kn = second derivatives at xn
	output:
	 A,B,C,D coefficients
	 
	"""
	def __init__(self,xdata,ydata,k0=0.0,kn=0.0):
		self.x  = xdata
		self.y  = ydata
		self.k0 = k0
		self.kn = kn
	#
	def spline(self):
		n = len(self.x)
		
		# Check if we have enough data points
		if n < 2:
			# Return empty lists for insufficient data
			return [], [], [], []
		# standard, POSITIVE interval widths h_i = x_{i+1} - x_i
		h = [self.x[i+1]-self.x[i] for i in range(n-1)]

		if n == 2:
			A_coef = [self.y[0]]
			B_coef = [(self.y[1]-self.y[0]) / h[0]]
			C_coef = [0.0]
			D_coef = [0.0]
			return A_coef, B_coef, C_coef, D_coef

		n2 = n-2   # number of interior nodes / unknown curvatures

		# ---- build the tridiagonal system for the interior curvatures ----
		# standard formula (Kiusalaas): for interior node i=1..n-2 (row = i-1):
		#   h_{i-1} k_{i-1} + 2(h_{i-1}+h_i) k_i + h_i k_{i+1}
		#       = 6[ (y_{i+1}-y_i)/h_i - (y_i-y_{i-1})/h_{i-1} ]
		b = [2.0*(h[i]+h[i+1]) for i in range(n2)]

		a_diag = []   # sub-diagonal  (multiplies k_{i-1})
		c_diag = []   # super-diagonal (multiplies k_{i+1})
		d_rhs  = []
		for row in range(n2):
			i = row+1
			if row-1 >= 0:
				a_diag.append(h[i-1])
			if row+1 < n2:
				c_diag.append(h[i])
			s_i   = (self.y[i+1]-self.y[i]) / h[i]
			s_im1 = (self.y[i]-self.y[i-1]) / h[i-1]
			d_rhs.append(6.0*(s_i - s_im1))

		# boundary curvature contributions moved to RHS
		d_rhs[0]  -= h[0]  * self.k0
		d_rhs[-1] -= h[n2] * self.kn

		if len(a_diag) == 0 and len(c_diag) == 0 and n2 > 1:
			# degenerate fallback (shouldn't normally trigger)
			A_coef, B_coef, C_coef, D_coef = [], [], [], []
			for i in range(n-1):
				A_coef.append(self.y[i])
				B_coef.append((self.y[i+1]-self.y[i]) / h[i])
				C_coef.append(0.0)
				D_coef.append(0.0)
			return A_coef, B_coef, C_coef, D_coef

		sigma_sol = MATH(a_diag, b, c_diag, d_rhs).tri_diagonal_solver()

		if len(sigma_sol) == 0:
			A_coef, B_coef, C_coef, D_coef = [], [], [], []
			for i in range(n-1):
				A_coef.append(self.y[i])
				B_coef.append((self.y[i+1]-self.y[i]) / h[i])
				C_coef.append(0.0)
				D_coef.append(0.0)
			return A_coef, B_coef, C_coef, D_coef

		sigma = [self.k0] + list(sigma_sol) + [self.kn]

		# ---- convert (per-interval) shifted-basis coefficients to the
		#      ABSOLUTE-x basis used throughout the rest of the codebase:
		#      Fun(r) = A_i + B_i*r + C_i*r^2 + D_i*r^3   for r in [x_i, x_{i+1}]
		A_coef = [0.0]*(n-1)
		B_coef = [0.0]*(n-1)
		C_coef = [0.0]*(n-1)
		D_coef = [0.0]*(n-1)
		for i in range(n-1):
			hi = h[i]
			xi = self.x[i]
			# shifted-basis coefficients (Taylor form about x_i)
			As = self.y[i]
			Bs = (self.y[i+1]-self.y[i])/hi - hi*(2.0*sigma[i]+sigma[i+1])/6.0
			Cs = sigma[i]/2.0
			Ds = (sigma[i+1]-sigma[i])/(6.0*hi)
			# convert to absolute-x polynomial coefficients
			D_coef[i] = Ds
			C_coef[i] = Cs - 3.0*Ds*xi
			B_coef[i] = Bs - 2.0*xi*Cs + 3.0*xi*xi*Ds
			A_coef[i] = As - Bs*xi + Cs*xi*xi - Ds*xi*xi*xi

		return A_coef, B_coef, C_coef, D_coef

#----------------#	
class SPLINE_INTEGRATE:
	"""
	calculate the integeration of spline function
	input:
	lower = lower item of integeration
    upper = upper limit of integeration
    a,b,c,d = spline coefficients
    rad_points = radial grid points list
	"""
	def __init__(self,lower,upper,a,b,c,d,rad_points):
		self.lower = lower
		self.upper = upper
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.rad_points = rad_points
	#
	def integral(self):
		out = 0
		# Check if coefficients are empty
		if len(self.a) == 0:
			return 0.0
		
		for i in range(len(self.rad_points)-1):
			ra = self.rad_points[i]
			rb = self.rad_points[i+1]
			
			rb2 = rb*rb
			rb3 = rb2*rb
			rb4 = rb3*rb
			#
			ra2 = ra*ra
			ra3 = ra2*ra
			ra4 = ra3*ra
			
			# Use i if within bounds, otherwise use last index
			idx = i if i < len(self.a) else len(self.a) - 1
			
			integ1 = (rb*self.a[idx]) + ((rb2*self.b[idx])/2.0) + \
			 ((rb3*self.c[idx])/3.0) + ((rb4*self.d[idx])/4.0) - \
				( (ra*self.a[idx]) + ((ra2*self.b[idx])/2.0) + \
				 ((ra3*self.c[idx])/3.0) + ((ra4*self.d[idx])/4.0) )
			#
			out = out + integ1
		return out

#----------------#
class NORM_INTEGERAL:
	"""
	to get the value of 
	int_{0}^{r_max} ( P(r)^2 +Q(r)^2 ) dr
	R is r values for all the range 
	"""
	def __init__(self,Pr_matchedx,Qr_matchedx,R):
		self.Pr_matchedx = Pr_matchedx
		self.Qr_matchedx = Qr_matchedx
		self.R = R
	#
	def normalization_integeral(self):
		# first find
		# Y = P(r)^2 + Q(r)^2
		R_new = []
		Y = []
		for i in range(len(self.Pr_matchedx)):
			q2 = self.Qr_matchedx[i]*self.Qr_matchedx[i]
			p2 = self.Pr_matchedx[i]*self.Pr_matchedx[i]
			y1 = p2 + q2
			Y.append(y1)
			# r values
			R_new.append(self.R[i])
		# express Y interms of spline coeficients
		# new spline coefficients
		A1x,B1x,C1x,D1x = SPLINE(R_new,Y,0.0,0.0).spline()
		
		# Check if spline coefficients are valid
		if len(A1x) == 0:
			# Fallback: simple trapezoidal integration
			integx = 0.0
			for i in range(len(R_new)-1):
				integx += 0.5 * (Y[i] + Y[i+1]) * (R_new[i+1] - R_new[i])
			return integx, R_new, [0.0], [0.0], [0.0], [0.0]
		
		itg = SPLINE_INTEGRATE(0,R_new[-1],A1x,B1x,C1x,D1x,R_new)
		integx = itg.integral()
		return integx,R_new,A1x,B1x,C1x,D1x

#----------------#
class INTEG_VALUE:
	"""
	Get the integeral value
	#
	#   Function for correction
	#
	"""
	def __init__(self,Pr_mat,Qr_mat,R_value):
		self.Pr_mat = Pr_mat
		self.Qr_mat = Qr_mat
		self.R_value = R_value
	#
	def integeral_value(self):
		"""
		to get the value of 
		int_{0}^{R} ( P(r)^2 +Q(r)^2 ) dr
		R is r values for all the range 
		
		R_value = is the Pr_mat values
		"""
		# Check if data is valid
		if len(self.Pr_mat) == 0 or len(self.R_value) == 0:
			return 0.0, [], [], [], []
		
		# first find
		# Y = P(r)^2 + Q(r)^2
		Y = []
		for i in range(len(self.Pr_mat)):
			q2 = self.Qr_mat[i]*self.Qr_mat[i]
			p2 = self.Pr_mat[i]*self.Pr_mat[i]
			y1 = p2 + q2
			Y.append(y1)
		
		# express Y interms of spline coeficients
		# new spline coefficients
		sp = SPLINE(self.R_value, Y)
		A1x, B1x, C1x, D1x = sp.spline()
		
		# Check if spline coefficients are valid
		if len(A1x) == 0:
			# Fallback: simple trapezoidal integration
			integxm = 0.0
			for i in range(len(self.R_value)-1):
				integxm += 0.5 * (Y[i] + Y[i+1]) * (self.R_value[i+1] - self.R_value[i])
			return integxm, [0.0], [0.0], [0.0], [0.0]
		
		itl = SPLINE_INTEGRATE(self.R_value[0], self.R_value[-1], A1x,
					B1x, C1x, D1x, self.R_value)
		integxm = itl.integral()
		return integxm, A1x, B1x, C1x, D1x

#----------------#
class PQ_NORMALIZED:
	"""
	the normalized values of P(r) and Q(r)
	"""
	def __init__(self,Pr_matched,Qr_matched,integ):
		self.Pr_matched = Pr_matched
		self.Qr_matched = Qr_matched
		self.integ = integ
	#
	def pq_normalized(self):
		"""
		the normalized values of P(r) and Q(r)
		"""
		# Check if integ is valid
		if self.integ <= 0:
			return self.Pr_matched, self.Qr_matched
		
		# P and Q normalized 
		Pr_norm = []
		Qr_norm = []
		#
		#factor = 1.0 / self.integ
		# modified 
		factor = 1.0 / (self.integ**(1.0/2.0) )
		for i in range(len(self.Pr_matched)):
			p_norm = self.Pr_matched[i]*factor
			Pr_norm.append(p_norm)
			#
			q_norm = self.Qr_matched[i]*factor
			Qr_norm.append(q_norm)
		return Pr_norm,Qr_norm

#----------------#
class EXPECTATION_R:
	"""
	the expectation value of r = < r >
	< r > = integ_{- inf}^{ + inf} P^{*}(r) r P(r) dr 
	P^{*}(r) is the complex conjugate
	"""
	def __init__(self,Pr_matchedx,R):
		self.Pr_matchedx = Pr_matchedx
		self.R = R
	def r_expectation(self):
		"""
		to get the value of 
		int_{0}^{r_max} ( P(r) * r * P(r)  ) dr
		R is r values for all the range 
		"""
		# first find
		# Y = P(r) * r * P(r)
		R_new = []
		Y = []
		for i in range(len(self.Pr_matchedx)):
			y1 = self.Pr_matchedx[i]*self.R[i]*self.Pr_matchedx[i]
			Y.append(y1)
			# r values
			R_new.append(self.R[i])
		#--------------#
		# express Y interms of spline coeficients
		# new spline coefficients
		sp = SPLINE(R_new, Y)
		A1x, B1x, C1x, D1x = sp.spline()
		
		# Check if spline coefficients are valid
		if len(A1x) == 0:
			# Fallback: simple trapezoidal integration
			integ_r = 0.0
			for i in range(len(R_new)-1):
				integ_r += 0.5 * (Y[i] + Y[i+1]) * (R_new[i+1] - R_new[i])
			return integ_r
		
		itl = SPLINE_INTEGRATE(R_new[0], R_new[-1], A1x,
					B1x, C1x, D1x, R_new)
		integ_r = itl.integral()
		return integ_r

#----------------#

class EXPECTATION_R_general:
    """
    Compute the relativistic expectation value

        <r> = int [ (P(r)^2 + Q(r)^2) r] dr

    where P(r) and Q(r) are the large and small components
    of the Dirac radial wavefunction.
    """

    def __init__(self, P, Q, R):
        self.P = P
        self.Q = Q
        self.R = R

    def r_expectation(self):
        # Construct integrand
        Y = []
        R_new = []

        for p, q, r in zip(self.P, self.Q, self.R):
            Y.append((p*p + q*q) * r)
            R_new.append(r)

        # Cubic spline interpolation
        sp = SPLINE(R_new, Y)
        A, B, C, D = sp.spline()

        # Fallback if spline fails
        if len(A) == 0:
            integ = 0.0
            for i in range(len(R_new)-1):
                h = R_new[i+1] - R_new[i]
                integ += 0.5 * (Y[i] + Y[i+1]) * h
            return integ

        # Spline integration
        itl = SPLINE_INTEGRATE(R_new[0],R_new[-1],A,B,C,D,R_new)

        return itl.integral()


#----------------------------#

class EXPECTATION_R_GENERAL_fx:
    """
    Compute the relativistic expectation value

        <f(r)> = int [ (P(r)^2 + Q(r)^2) f(r) ] dr

    where

        P(r) : large component
        Q(r) : small component
        f(r) : arbitrary radial operator

    Examples
    --------
    <r>      : operator = lambda r: r
    <r^2>    : operator = lambda r: r**2
    <1/r>    : operator = lambda r: 1.0/r
    <r^k>    : operator = lambda r: r**k
    """

    def __init__(self, P, Q, R,
                 operator=lambda r: r,
                 skip_origin=True,
                 eps=1.0e-15):

        self.P = P
        self.Q = Q
        self.R = R
        self.operator = operator
        self.skip_origin = skip_origin
        self.eps = eps

    def r_expectation(self):

        Y = []
        R_new = []

        for p, q, r in zip(self.P, self.Q, self.R):

            # Skip the origin if requested
            if self.skip_origin and abs(r) < self.eps:
                continue

            # Evaluate the operator safely
            try:
                op = self.operator(r)
            except (ZeroDivisionError, OverflowError):
                continue

            Y.append((p*p + q*q) * op)
            R_new.append(r)

        if len(R_new) < 2:
            raise ValueError(
                "Not enough grid points to evaluate the expectation value."
            )

        # Cubic spline interpolation
        sp = SPLINE(R_new, Y)
        A, B, C, D = sp.spline()

        # Fallback to trapezoidal integration
        if len(A) == 0:
            integ = 0.0
            for i in range(len(R_new)-1):
                h = R_new[i+1] - R_new[i]
                integ += 0.5 * (Y[i] + Y[i+1]) * h
            return integ

        # Spline integration
        itl = SPLINE_INTEGRATE(R_new[0],R_new[-1],A,B,C,D,R_new)

        return itl.integral()

#-----------------------------#

# for interpolation of reading the potential function 
class Cubic_Spline_Interpolation_Function:
	def __init__(self, x, a, b, c, d):
		self.x = x
		self.a = a
		self.b = b
		self.c = c
		self.d = d

	def _find_interval(self, value):
		"""
		Return interval index such that x[i] <= value < x[i+1].
		"""
		n = len(self.x)

		if n == 0:
			return 0
		if value <= self.x[0]:
			return 0
		if value >= self.x[-1]:
			return n - 2 if n >= 2 else 0

		low = 0
		high = n - 1

		while high - low > 1:
			mid = (low + high) // 2
			if value < self.x[mid]:
				high = mid
			else:
				low = mid

		return low

	def __call__(self, x_new):
		# Check if data is valid
		if len(self.x) == 0 or len(self.a) == 0:
			return [0.0] * len(x_new) if isinstance(x_new, (list, tuple)) else 0.0

		# Single point
		if isinstance(x_new, (int, float)):
			i = self._find_interval(x_new)
			if i >= len(self.a):
				i = len(self.a) - 1
			if i < 0:
				i = 0
			# FIX: evaluate using ABSOLUTE x, matching PQ_MIDDLE/PQ_ZERO/
			# SPLINE_INTEGRATE's basis -- no shift by self.x[i].
			x = x_new
			return (self.a[i]
					+ self.b[i] * x
					+ self.c[i] * x * x
					+ self.d[i] * x * x * x)

		# List or tuple of points
		result = []
		for x in x_new:
			i = self._find_interval(x)
			if i >= len(self.a):
				i = len(self.a) - 1
			if i < 0:
				i = 0
			y = (self.a[i]
				 + self.b[i] * x
				 + self.c[i] * x * x
				 + self.d[i] * x * x * x)
			result.append(y)

		return result
