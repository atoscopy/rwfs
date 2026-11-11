#!/usr/bin/python3

import math
from constants import half,third,pi,e,a0_fm,fermi_skin_t_fm


# 1 Thomas- Fermi

class THOMAS_FERMI:
	"""
	The GRID class defines the used grid for our numerical calculations	
	# potential equal -rV(r)
	-rV(r) = Z * phi - z_inf ( sqrt( phi ) - 1 ) + V_0
	where
	z_inf = Z+1-N_e, 
	and V_0 = (r/r_TF) * Z^{1/3} * phi 
	INPUT:
		r = grid points in list format
		Z = atomic number
		N_e = number of electrons in atom/ion
	"""
	def __init__(self, r_grid, Z, N_e):
		self.r 		= r_grid
		self.Z 	    = Z
		self.N 		= N_e
	def z_inf(self):
		return self.Z + 1 - self.N
	def phi(self):
		phi = []
		const1 = (self.N -1.0)**third
		const2 = (1.0/4.0)*(((9*pi*pi)/2.0)**third)
		for i in range(len(self.r)):
			r = self.r[i]
			m = const1 * r / const2
			# T.F. potential as in GRASP 
			dm = 1.0+(1.81061*(m**half))+(0.60112*m)
			rt = 1.0+(1.81061*(m**half))+(1.39515*m)+\
			(0.77112*(m**(3.0/2.0)))+(0.21465*(m*m))+\
			(0.04793*(m**(5.0/2.0)))
			phi_m = (dm/rt)**2
			phi.append(phi_m)
		return phi
	# old function NOT used
	def V_node(self):
		r_TF = (1.0/4.0)*(((9*pi*pi)/(2.0*self.Z))**third)
		v0 = []
		phi = self.phi()
		for m in range(len(self.r)):
			v = (self.r[m]/r_TF) * (self.Z**(third)) * phi[m]
			v0.append(v)
		return v0
	def potential(self):
		rv = []
		v0  = self.V_node()
		phi = self.phi()
		z_inf = self.z_inf()
		for m in range(len(self.r)):
			rv0 = (self.Z*phi[m]) -( z_inf*phi[m]) +\
					z_inf 
			rv.append(-rv0)
		return rv

# 2- Coulomb potential

class COULOMB:
	"""
	The GRID class defines the used grid for our numerical calculations	
	
	V(r) = -Z/r
	
	INPUT:
		r = grid points in list format
		Z = atomic number
	"""
	def __init__(self, r_grid, Z):
		self.r 		= r_grid
		self.Z 	    = Z
	def potential(self):
		rv = []
		for m in range(len(self.r)):
			rv0 = -self.Z
			rv.append(rv0)
		return rv
	




# 3- Yukawa potential

class YUKAWA:
	"""
	The GRID class defines the used grid for our numerical calculations	
	Yukawa 
	
	V(r) = -Z/r e^{-alpha r} 
	
	INPUT:
		r = grid points in list format
		Z = atomic number
		alpha = screening constant 
	"""
	def __init__(self, r_grid, Z,  alpha):
		self.r 		= r_grid
		self.Z 	    = Z
		self.alpha  = alpha
	#-----#
	def potential(self):
		rv = []
		for m in range(len(self.r)):
			if float(self.r[m]) == 0.0:
				rv0 = 0.0
				rv.append(rv0)
			else:
				t1 = -self.Z/self.r[m]
				t2 = e**( -self.alpha *self.r[m] )
				#
				rv0 = self.r[m] * ( t1 * t2)
				rv.append(rv0)
		return rv
		
		
# 4- Multi Yukawa potential

class MULTI_YUKAWA:
	"""
	The GRID class defines the used grid for our numerical calculations	
	multi Yukawa 
	
	V(r) = -Z/r sum_{i}^{n_max} c_i e^{-alpha_i r} 
	
	INPUT:
		r = grid points in list format
		Z = atomic number
		terms = [(c_1,alpha_1),(c_2,alpha_2),(c_3,alpha_3), ....... ]
	"""
	def __init__(self, r_grid, Z ,n_max, terms):
		self.r 		= r_grid
		self.Z 	    = Z
		self.nmax  = n_max
		self.term  = terms
	#-----#
	def potential(self):
		rv = []
		for m in range(len(self.r)):
			sum_exp = 0.0
			for k in range(self.nmax):
				item  = self.term[k]
				c = float(item[0])
				a = float(item[1])
				tt = c * e**( -a *self.r[m] )
				sum_exp = sum_exp + tt
			#
			if float(self.r[m]) == 0.0:
				rv0 = 0.0
				rv.append(rv0)
			else:
				rv0 = self.r[m] * ( (-self.Z/self.r[m]) * sum_exp)
				rv.append(rv0)
		return rv
		
		
# 5- Screened potential

class SCREENED:
	"""
	The GRID class defines the used grid for our numerical calculations	
	multi Yukawa 
	
	V(r) = -1/r [ 1 + (Z-1)e^{-alpha r} ]
	
	INPUT:
		r = grid points in list format
		Z = atomic number
		alpha = screening parameters
		
	"""
	def __init__(self, r_grid, Z, alpha):
		self.r 		= r_grid
		self.Z 	    = Z
		self.alpha  = alpha
	#-----#
	def potential(self):
		rv = []
		for m in range(len(self.r)):
			if float(self.r[m]) == 0.0:
				rv0 = 0.0
				rv.append(rv0)
			else:
				t1 =  (self.Z - 1) * (e**(-self.alpha * self.r[m]) )
				rv0 = self.r[m] * ( (-1/self.r[m]) * (1 + t1) )
				rv.append(rv0)
		return rv

# 6- Screened Hydrogenic potential

class SCREENED_HYDROGENIC:
	"""
	The GRID class defines the used grid for our numerical calculations	
	multi Yukawa 
	
	V(r) = - Z_{eff}/r  == - (Z-alpha)/r
	
	INPUT:
		r = grid points in list format
		Z = atomic number
		alpha = screening parameters
		
	"""
	def __init__(self, r_grid, Z, alpha):
		self.r 		= r_grid
		self.Z 	    = Z
		self.alpha  = alpha
	#-----#
	def potential(self):
		rv = []
		for m in range(len(self.r)):
			if float(self.r[m]) == 0.0:
				rv0 = 0.0
				rv.append(rv0)
			else:
				zeff =  (self.Z - self.alpha) 
				rv0 = self.r[m] *  (- zeff /self.r[m])  
				rv.append(rv0)
		return rv
		
# 7- read from input file
class READ_pote_file_rv:
	"""
	read potential data from file with header
	
	r ,   v(r) in the first line 
	"""
	def __init__(self, file_name):
		self.file  = file_name
		
	#-----#
	def read(self):
		#
		r_new = []
		v_new = []
		#
		num_lev_mx = 0
		with open(self.file) as klo:  
			for line in klo:
				num_lev_mx += 1
				line1 = line.rstrip().split('\n')
				for aj in range(len(line1)):
					line2 = line1[aj].replace(" ",'\n')
					line3 = line2.rstrip().split()
					#
					if len(line3)> 0 and num_lev_mx> 1:
						r_new.append(float(line3[0]))
						v_new.append(float(line3[1]))
		#
		return r_new,v_new



########################################################################
# FINITE NUCLEAR SIZE MODELS
#   8 - uniformly charged sphere
#   9 - Fermi (two-parameter Woods-Saxon) charge distribution
########################################################################
#
# Atomic units (hbar = m_e = e = 1);  1 fm = 1/a0_fm a.u.
#
# A spherically symmetric nuclear charge density rho(r), normalized to
# the nuclear charge Z,
#
#     4 pi Int_0^inf rho(r) r^2 dr = Z                              (N)
#
# gives the potential energy of an electron
#
#     V(r) = - Int rho(r') / |r - r'| d^3r'
#          = - [ Q(r)/r + S(r) ]                                    (P)
#
#     Q(r) = 4 pi Int_0^r   rho(r') r'^2 dr'   (charge inside r)
#     S(r) = 4 pi Int_r^inf rho(r') r'   dr'   (potential of the
#                                               shell outside r)
#
# RWFS works with the function rV(r), which is what every potential()
# method returns:
#
#     rV(r) = - [ Q(r) + r S(r) ]                                   (RV)
#
# For r beyond the nucleus Q = Z and S = 0, so rV = -Z (point Coulomb
# tail), while at the origin rV(0) = 0, i.e. V(r) stays finite.
# Only the bare nuclear potential is generated: electron screening is
# NOT included (as for the Coulomb option 2).
#
#----------------#
def rms_radius_from_mass_number(A):
    """
    approximate rms nuclear charge radius (fm) from the mass number A
    (Johnson and Soff, At. Data Nucl. Data Tables 33, 405-446 (1985)):

        R_rms = 0.836 A^(1/3) + 0.570   fm,      valid for A > 9

    This is only a smooth global fit; use a measured radius when
    one is available.
    """
    return 0.836*(float(A)**third) + 0.570


#----------------#
def _gauss_legendre(n):
    """
    nodes x_i and weights w_i of the n-point Gauss-Legendre rule on
    [-1,1] (Newton iteration on the Legendre polynomial P_n)
    """
    xs = []
    ws = []
    for i in range(1,n+1):
        x = math.cos(math.pi*(i-0.25)/(n+0.5))
        for it in range(100):
            p0 = 1.0
            p1 = x
            for k in range(2,n+1):
                p0,p1 = p1,((2*k-1)*x*p1-(k-1)*p0)/k
            dp = n*(x*p1-p0)/(x*x-1.0)
            dx = p1/dp
            x = x - dx
            if abs(dx) < 1.0e-15:
                break
        p0 = 1.0
        p1 = x
        for k in range(2,n+1):
            p0,p1 = p1,((2*k-1)*x*p1-(k-1)*p0)/k
        dp = n*(x*p1-p0)/(x*x-1.0)
        xs.append(x)
        ws.append(2.0/((1.0-x*x)*dp*dp))
    return xs,ws

_GL_X,_GL_W = _gauss_legendre(8)


#----------------#
# 8 - uniformly charged sphere
class UNIFORM_SPHERE:
    """
    Uniformly charged sphere of radius R

        rho(r) = 3 Z / (4 pi R^3)     r <= R
               = 0                    r >  R

    The rms radius of the charge distribution is

        <r^2> = (3/5) R^2        ==>    R = sqrt(5/3) R_rms

    Integrating (P) gives

        V(r)  = - Z/(2R) [ 3 - (r/R)^2 ]        r <= R
              = - Z / r                         r >  R

        rV(r) = - Z [ 3r/(2R) - r^3/(2R^3) ]    r <= R
              = - Z                             r >  R

    rV and its first derivative are continuous at r = R.

    INPUT:
        r_grid   = grid points in list format (a.u.)
        Z        = atomic number
        R_rms_fm = rms nuclear charge radius (fm)
    """
    def __init__(self, r_grid, Z, R_rms_fm):
        self.r        = r_grid
        self.Z        = Z
        self.R_rms_fm = float(R_rms_fm)
        # radius of the sphere in a.u.
        self.R        = math.sqrt(5.0/3.0)*self.R_rms_fm/a0_fm
    #-----#
    def potential(self):
        rv = []
        Z = float(self.Z)
        R = self.R
        for x in self.r:
            if x >= R:
                rv.append(-Z)
            else:
                s = x/R
                rv.append(-Z*(1.5*s-0.5*s*s*s))
        return rv
    #-----#
    def constants(self):
        """ list of (name,value) pairs describing the model """
        return [("R_rms (fm)",self.R_rms_fm),\
                ("R sphere (fm)",self.R*a0_fm),\
                ("R sphere (a.u.)",self.R)]


#----------------#
# 9 - Fermi charge distribution
class FERMI_NUCLEUS:
    """
    Two-parameter Fermi (Woods-Saxon) charge distribution

        rho(r) = rho_0 / ( 1 + exp[(r-c)/a] ) = rho_0 f(r)          (F1)

    c = half-density radius,
    a = diffuseness, fixed by the 10%-90% surface thickness t,

        t = 4 a ln 3     ==>   a = t/(4 ln 3)                       (F2)

    with the standard t = 2.3 fm (a = 0.5234 fm).

    --- moments --------------------------------------------------------
    M_n = Int_0^inf r^n f(r) dr = - n! a^(n+1) Li_(n+1)( -exp(c/a) )
    and with the inversion formulae of the polylogarithm

        M_1 = c^2/2 + pi^2 a^2/6               - a^2   e_2
        M_2 = c^3/3 + pi^2 a^2 c/3             + 2 a^3 e_3
        M_4 = c^5/5 + 2 pi^2 a^2 c^3/3
                    + 7 pi^4 a^4 c/15          + 24 a^5 e_5         (F3)

        e_s = sum_{k=1}^inf (-1)^(k+1) exp(-k c/a) / k^s

    (the e_s terms are of order exp(-c/a) ~ 1e-5 and are kept).

    --- normalization (N) ----------------------------------------------
        rho_0 = Z / (4 pi M_2)                                      (F4)

    --- rms radius -----------------------------------------------------
        <r^2> = M_4 / M_2
              = (3/5) c^2 + (7/5) pi^2 a^2   (+ terms ~ exp(-c/a))  (F5)

    For a given R_rms the parameter c is obtained by solving
    M_4(c)/M_2(c) = R_rms^2; the starting value is

        c = sqrt( (5/3) R_rms^2 - (7/3) pi^2 a^2 )                  (F6)

    (this needs R_rms > sqrt(7/5) pi a = 1.95 fm for t = 2.3 fm; use
    the uniform sphere for lighter nuclei).

    --- potential ------------------------------------------------------
    Inserting (F1) and (F4) in (RV),

        rV(r) = - (Z/M_2) [ I_2(r) + r J_1(r) ]                     (F7)

        I_2(r) = Int_0^r   x^2 f(x) dx
        J_1(r) = Int_r^inf x   f(x) dx

    I_2 and J_1 are evaluated by 8-point Gauss-Legendre quadrature on
    sub-panels no wider than a/2 (the integrand varies on the scale a),
    accumulated forward (I_2) and backward (J_1) over the grid. Beyond
    r_cut = c + 60 a, f < 1e-26 and rV = -Z is used exactly.

    INPUT:
        r_grid   = grid points in list format (a.u., increasing)
        Z        = atomic number
        R_rms_fm = rms nuclear charge radius (fm)
        t_fm     = 10%-90% skin thickness (fm), default 2.3
    """
    def __init__(self, r_grid, Z, R_rms_fm, t_fm=fermi_skin_t_fm):
        self.r        = r_grid
        self.Z        = Z
        self.R_rms_fm = float(R_rms_fm)
        self.t_fm     = float(t_fm)
        # diffuseness a in a.u.   (F2)
        self.a        = self.t_fm/(4.0*math.log(3.0))/a0_fm
        # half-density radius c in a.u.  (F5,F6)
        self.c        = self._solve_c()
        self.M2,self.M4 = self._moments(self.c)
        self.rho0     = float(self.Z)/(4.0*pi*self.M2)      # (F4)
        self.r_cut    = self.c + 60.0*self.a
    #-----#
    def _eps(self,s,c):
        # e_s = sum_k (-1)^(k+1) exp(-k c/a)/k^s
        x = c/self.a
        total = 0.0
        for k in range(1,2000):
            term = math.exp(-k*x)/float(k**s)
            if term < 1.0e-20:
                break
            if k % 2 == 1:
                total += term
            else:
                total -= term
        return total
    #-----#
    def _moments(self,c):
        # M_2 and M_4 of (F3)
        a = self.a
        e3 = self._eps(3,c)
        e5 = self._eps(5,c)
        M2 = c**3/3.0 + pi*pi*a*a*c/3.0 + 2.0*a**3*e3
        M4 = c**5/5.0 + 2.0*pi*pi*a*a*c**3/3.0 + \
             7.0*pi**4*a**4*c/15.0 + 24.0*a**5*e5
        return M2,M4
    #-----#
    def _solve_c(self):
        # solve M_4(c)/M_2(c) = R_rms^2 by bisection (monotonic in c)
        a  = self.a
        R2 = (self.R_rms_fm/a0_fm)**2
        c2 = (5.0/3.0)*R2 - (7.0/3.0)*pi*pi*a*a
        if c2 <= 0.0:
            rmin = math.sqrt(7.0/5.0)*pi*a*a0_fm
            raise ValueError(" ERROR: R_rms = "+str(self.R_rms_fm)+\
              " fm is too small for a Fermi distribution with t = "+\
              str(self.t_fm)+" fm (need R_rms > "+str(round(rmin,3))+\
              " fm). Use the uniform sphere (option 8) instead.")
        c0 = math.sqrt(c2)
        lo = 0.0
        hi = 2.0*c0 + 10.0*a
        for it in range(200):
            mid = 0.5*(lo+hi)
            m2,m4 = self._moments(mid)
            if m4/m2 < R2:
                lo = mid
            else:
                hi = mid
        return 0.5*(lo+hi)
    #-----#
    def _fermi(self,x):
        u = (x-self.c)/self.a
        if u > 0.0:
            ex = math.exp(-u)
            return ex/(1.0+ex)
        return 1.0/(1.0+math.exp(u))
    #-----#
    def _panel(self,x0,x1):
        # (Int x^2 f dx , Int x f dx) over [x0,x1]
        if x1 > self.r_cut:
            x1 = self.r_cut
        if x0 >= x1:
            return 0.0,0.0
        nsub = int(math.ceil((x1-x0)/(0.5*self.a)))
        if nsub < 1:
            nsub = 1
        h = (x1-x0)/nsub
        s2 = 0.0
        s1 = 0.0
        for j in range(nsub):
            xm = x0 + (j+0.5)*h
            for xi,wi in zip(_GL_X,_GL_W):
                x  = xm + 0.5*h*xi
                fw = self._fermi(x)*wi
                s1 += fw*x
                s2 += fw*x*x
        return 0.5*h*s2,0.5*h*s1
    #-----#
    def potential(self):
        r = self.r
        n = len(r)
        Z = float(self.Z)
        # first grid point at or beyond r_cut
        m = n-1
        for i in range(n):
            if r[i] >= self.r_cut:
                m = i
                break
        # integrals over each grid interval
        d2 = [0.0]*n
        d1 = [0.0]*n
        for i in range(m):
            d2[i],d1[i] = self._panel(r[i],r[i+1])
        # I_2(r) forward, J_1(r) backward
        I2 = [0.0]*n
        for i in range(m):
            I2[i+1] = I2[i] + d2[i]
        J1 = [0.0]*n
        J1[m] = self._panel(r[m],self.r_cut)[1]
        for i in range(m-1,-1,-1):
            J1[i] = J1[i+1] + d1[i]
        # rV(r) of (F7)
        rv = []
        for i in range(n):
            if r[i] < self.r_cut:
                rv.append( -(Z/self.M2)*(I2[i]+r[i]*J1[i]) )
            else:
                rv.append(-Z)
        return rv
    #-----#
    def constants(self):
        """ list of (name,value) pairs describing the model """
        return [("R_rms (fm)",self.R_rms_fm),\
                ("t (fm)",self.t_fm),\
                ("c (fm)",self.c*a0_fm),\
                ("a (fm)",self.a*a0_fm)]
	


