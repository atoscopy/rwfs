#!/usr/bin/python3

from constants import c_au, half,R_WKB_parameter
from pq_middle import PQ_MIDDLE
from convergence import CONVEGE


class R_INF:
	"""
	defines the last grid point for each orbital defined by equation 44
	input:
	 l = orbital angular momentum
	 energy = initial energy estimate
	 r = grid points
	 A,B,C,D = spline coefficients
	 IOTP = last point from outward solution (turning point)
	 kappa = kappa quantum number for orbital
	"""
	def __init__(self,l,energy,r,A,B,C,D,IOTP,kappa):
		self.L1 = l
		self.E1 = energy
		self.R1 = r
		self.A1 = A
		self.B1 = B
		self.C1 = C
		self.D1 = D
		self.IOTP1 = IOTP
		self.kappa = kappa
	#		
	def r_infinity(self):
		sl2 = c_au*c_au
		par1 = (self.E1+(2*sl2))/(sl2)
		par2 = self.L1*(self.L1+1)
		#
		square_rmu  = []
		r_infty_idx = []
		r_inf = -1
		F_at_r_inf = None
		#
		F_idx       = []
		F_constant  = []
		#
		grid = len(self.R1)
		for i in range(grid-1,0,-1):
			m = i       # index of radial points 
			s = i-1     # spline index
			#
			r2 = self.R1[m]*self.R1[m]
			r3 = r2*self.R1[m]
			#
			V_i = ( self.A1[s]+(self.R1[m]*self.B1[s])+\
			      (r2*self.C1[s])+(r3*self.D1[s]) )/ self.R1[m]
			V_i_prime = ( -self.A1[s]+ (self.C1[s]*self.R1[m]) +\
			      (2*self.D1[s]*r3) ) /r2
			val1 = par1*( V_i - self.E1)
			val2 = par2/r2
			#
			inside1 = val1+val2  
			r_mu_sq = r2*inside1
			square_rmu.append(r_mu_sq)
			# F
			## NOTE
			#mg = 35.0*35.0
			mg = float(R_WKB_parameter)*float(R_WKB_parameter)
			#
			
			if r_mu_sq < mg and m >= self.IOTP1:
				r_infty_idx.append(m)
				if r_mu_sq > 0.0:
					F = 0.5*( self.B1[s] / (self.R1[m]*r_mu_sq) )-\
					   ( (0.25*r2*par1*V_i_prime)/r_mu_sq )-\
					   ( (r_mu_sq**half) /self.R1[m] )
					if r_inf == -1:          # first point where BOTH conditions hold
						r_inf = m
						F_at_r_inf = F
					F_idx.append(m)
					F_constant.append(F)
				
		#
		if r_inf > 0 and F_at_r_inf is not None:
			par3 = c_au /( self.E1+(2*c_au*c_au) )
			Pro = 1.0
			Qro = par3*( (self.kappa/self.R1[r_inf]) + F_at_r_inf )
			return r_inf,Pro,Qro
		else:
			print(" ERROR: (r_infinity) grid does not extend far enough for the")
			print(" asymptotic (WKB) matching condition — increase np/dr or rmax")
			return 0,0,0
		
#----------------#
class PQ_INWARD:
	"""
	Inward solution
	input:
	r_inf_idx = index of the last grid point for the orbital
	P0 , Q0 = the initial estimation of P,Q
	r = grid points
	A,B,C,D = spline coefficients
	IOTP = last point from outward solution (turning point)
	kappa = kappa quantum number for orbital
	energy = initial energy estimate 
	eta = estimated error
	"""
	def __init__(self,r_inf_idx,P0,Q0,r,A,B,C,D,IOTP,kappa,eta_error,E_init):
		self.r_infty_index = r_inf_idx
		self.P_r0_init = P0
		self.Q_r0_init = Q0
		self.R = r
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.IOTP = IOTP
		self.K = kappa
		self.eta = eta_error
		self.E = E_init
	#	  
	def pq_inward(self):
		Pr_inward = []
		Qr_inward = []
		# initial estimates
		itm_idx = self.r_infty_index
		#
		Pr_inward.append(self.P_r0_init)
		Qr_inward.append(self.Q_r0_init)
		# remain P,Q
		# P and Q values with the initial values
		n_inward = self.r_infty_index - self.IOTP
		r_inf = self.r_infty_index
		#
		for i in range(n_inward):
			# 1 values of P(x), Q(x) as P_rb,Q_rb
			initial   = r_inf - i              #  = i
			final     = r_inf - (i+1)          #  = j        
			#
			P_ra_in = Pr_inward[-1]
			Q_ra_in = Qr_inward[-1]
			#
			pqm = PQ_MIDDLE(initial,final,self.R,self.A,self.B,self.C,\
						self.D,self.IOTP,self.K,P_ra_in,Q_ra_in,\
						self.eta,self.E)
			P_rb,Q_rb,converg =  pqm.pq_middle()
			# check convergence
			if converg == 0:
				ran = self.R[initial]
				rbn = self.R[final]
				#
				gconv = CONVEGE(initial,final,ran,rbn,P_ra_in,\
						Q_ra_in,self.A,self.B,self.C,self.D,self.R,\
						self.IOTP,self.K,self.eta,self.E)
				P_rb,Q_rb,convergence = gconv.get_converg()
				#
				if convergence == 0:
					print(" ERROR : (pq_inward) check the fn.")
					break
			#
			Pr_inward.append(P_rb)
			Qr_inward.append(Q_rb)
			#
		# reverse order
		Pr_in = []
		Qr_in = []
		#
		count2 = len(Pr_inward)
		for t in range( count2 ):
			idx = count2-(t+1)  # to get reverse count for inward list
			Pr_in.append(Pr_inward[idx])
			Qr_in.append(Qr_inward[idx])
		# remove the points at turning points
		P_inw = Pr_in[0]
		Q_inw = Qr_in[0]
		P_inward = Pr_in[1:]
		Q_inward = Qr_in[1:]
		# end
		return P_inward,Q_inward,P_inw,Q_inw
		
#----------------#
class ENERGY_CORRECTION:
	"""
	This calculates the energy correction defined by equation 45
	"""
	def __init__(self,integr,P_outr,Q_outr,P_inwr,Q_inwr):
		self.integr = integr
		self.P_outr = P_outr
		self.Q_outr = Q_outr
		self.P_inwr = P_inwr
		self.Q_inwr = Q_inwr
	#
	def energy_correction(self):
		"""
		calculates the energy correction for binding energy 
		"""
		# energy correction
		# scaling of the Q inward
		fact1 = self.P_outr/self.P_inwr
		Q_inward = self.Q_inwr * fact1
		#
		integ_x = self.integr
		#
		MM = (self.P_outr*(self.Q_outr-Q_inward)) / integ_x   
		DE = c_au * MM 
		return DE
