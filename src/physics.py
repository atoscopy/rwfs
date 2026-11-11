#!/usr/bin/python3

from constants import c_au          # spped of light

class PHYSICS:
	"""
	this class defines some physical functions 
	such as:
	 angular_momentum = gets l from kappa
	 actual_nodes = gets node count of orbital
	"""
	def __init__(self,n,kappa):
		self.n = n
		self.k = kappa
	#
	def angular_momentum(self):
		if self.k < 0 :
			l = - self.k - 1
		else:
			l = self.k
		return l
	#
	def actual_nodes(self):
		l = self.angular_momentum()
		nodx = self.n-l-1
		return nodx
#----------------#
class BOUND_EMIN_EMAX:
	"""
	defines the minimum and maximum energy of the bound orbitals
	input:
	 n = principle quantum number
	 kappa: orbital QN
	 rv = -r*V(r) is potential energy
	 r  = radial grid points
	 e_initial = initial estimate of the bound energy
	 
	 classes used:
	  Thomas_Fermi (rV)
	"""
	def __init__(self,n,kappa,rv,r,e_initial):
		self.n  = n
		self.k  = kappa
		self.rv = rv
		self.r  = r
		self.E  = e_initial
	#
	def e_min_and_e_max(self):
		l = PHYSICS(self.n,self.k).angular_momentum()
		const1 = 0.5*l*(l+1)  
		# 1- determine E_min
		emin = 1.0
		np = len(self.r)
		for i in range(1,np):
			r2 = self.r[i]*self.r[i]
			term1 = (self.rv[i]/self.r[i])+(const1/r2)
			emin= min(emin,term1)
		# to avoid large value
		if emin < -c_au*c_au:
			emin = -c_au*c_au 
		EMIN = emin 
		#------------------#
		# 2- determine E_max
		last_np = np-1
		r2 = self.r[last_np]*self.r[last_np]
		term2 = ( self.rv[last_np]/self.r[last_np] )+(const1/r2)
		emax = min( term2 , 0.0)  
		EMAX = emax
		#------------------#
		# confinement of the initial estimate of energy
		E_av = self.E
		if self.E > EMAX or self.E < EMIN:
			E_av = 0.5*(EMAX+EMIN)
		#
		return EMIN,EMAX,E_av
#----------------#
class CONDN1:
	"""
	check the last term conditions 
	 condition1 = CONDN1
	 condition2 = CONDN2
	"""
	def __init__(self,an,bn,phi_x1,psi_x1,eta_error):
		self.a_n 	= an
		self.b_n 	= bn
		self.phi_x1 = phi_x1
		self.psi_x1 = psi_x1
		self.eta 	= eta_error
	#
	def cond1(self):
		mod_a = abs(self.a_n)
		mod_b = abs(self.b_n)
		term_left = max(mod_a,mod_b)
		#
		mod_px1 = abs(self.phi_x1)
		mod_qx1 = abs(self.psi_x1)
		term_right = self.eta *( max(mod_px1,mod_qx1) )
		if term_left < term_right:
			return True
		else:
			return False	
#
#----------------#
class CONDN2:
	"""
	check the last term conditions 
	condition1 = CONDN1
	condition2  = CONDN2
	"""
	def __init__(self,rj,kappa,ucr,theta,H,phi_x1,d_phi_x1,\
				psi_x1,d_psi_x1,eta_error):
		self.r_j	= rj
		self.k 		= kappa
		self.ucr 	= ucr
		self.theta 	= theta
		self.h 		= H
		self.px1 	= phi_x1
		self.d_px1 	= d_phi_x1
		self.qx1 	= psi_x1
		self.d_qx1 	= d_psi_x1
		self.eta 	= eta_error
	#
	def cond2(self):
		left1 = (self.r_j*self.d_px1)+(self.k*self.h*self.px1)+\
				(self.ucr*self.h*self.qx1)
		mod_left1 = abs(left1)
		left2 = (self.r_j*self.d_qx1)-(self.k*self.h*self.qx1)-\
				(self.theta*self.h*self.px1)
		mod_left2 = abs(left2)
		term_left = max(mod_left1 , mod_left2)
		#
		mod_px1 = abs(self.px1)
		mod_qx1 = abs(self.qx1)
		term_right = self.eta*(max(mod_px1,mod_qx1))
		if term_left < term_right:
			return True
		else:
			return False
#----------------#
# class 4
class KAPPA_SIGN:
	"""
	to get the sign of the kappa
	if k = -1,-2,-3, ...  return -ve
	else:  return +ve
	"""
	def __init__(self,kappa):
		self.k = kappa
	#
	def kappa_sgn(self):
		if self.k > 0:
			return "+ve"
		elif self.k == 0:
			print("error in kappa value kappa = ",self.k)
			return "error in kappa value"
		else:
			return "-ve"
#----------------#
class INWARD_OUTWARD_DATA:
	"""
	collect all data of inward and outward results in order to 
	facilitate all the processes
	input:
	r = grid points
	pr_out = outward solution of P(r)
	qr_out = outward solution of Q(r)
	pr_in = inward solution of P(r)
	qr_in = inward solution of Q(r)
	P_inw,Q_inw = are the matching points of inward and outward solution 
	"""
	# Getting the outward and the inward integerations
	# Used for new corrections
	def __init__(self,r,pr_out,qr_out,pr_in,qr_in,P_inw,Q_inw):
		self.R = r
		self.Pr_outward = pr_out
		self.Qr_outward = qr_out
		self.Pr_inward = pr_in
		self.Qr_inward = qr_in
		self.P_inw = P_inw
		self.Q_inw = Q_inw
	#
	def r_pq_inward_outward(self ):
		# outward
		r_out = []
		p_out = []
		q_out = []
		#
		item = -1
		for i in range(len(self.Pr_outward)):
			r_out.append(self.R[i])
			p_out.append(self.Pr_outward[i])
			q_out.append(self.Qr_outward[i])
			item+=1
		# inward
		r_in = []
		p_in = []
		q_in = []
		for j in range(len(self.Pr_inward)):
			k = item+1+j
			r_in.append(self.R[k])
			p_in.append(self.Pr_inward[j])
			q_in.append(self.Qr_inward[j])
		# all data for outward
		rpq_out = [r_out,p_out,q_out]
		# all data for inward
		rpq_in = [r_in,p_in,q_in]
		return rpq_out,rpq_in
#----------------#
class RESCALING_INWARD:
	"""
	rescaling the inward solution with the outward
	input:
	IOTP = turning point
	Pr_outward, Qr_outward = the outward solution of P,Q
	Pr_inward, Qr_inward = the inward solution of P,Q
	P_out = outward last P(r) value at IOTP
	P_inw = inward last P(r) value at IOTP
	"""
	def __init__(self,IOTPx,Pr_outwardx, Qr_outwardx,P_outx,Pr_inwardx,\
			Qr_inwardx,P_inwx):
		self.IOTPx = IOTPx
		self.Pr_outwardx = Pr_outwardx
		self.Qr_outwardx = Qr_outwardx
		self.P_outx = P_outx
		self.Pr_inwardx = Pr_inwardx
		self.Qr_inwardx = Qr_inwardx
		self.P_inwx = P_inwx
	#
	def pr_qr_matched(self):
		Pr_matched = []
		Qr_matched = []
		# outward
		for i in range( len(self.Pr_outwardx) ):
			Pr_matched.append(self.Pr_outwardx[i])
			Qr_matched.append(self.Qr_outwardx[i])
		# inward
		fact = self.P_outx/self.P_inwx
		for k in range(len(self.Pr_inwardx)):
			pr_i =  self.Pr_inwardx[k]*fact
			qr_i =  self.Qr_inwardx[k]*fact
			Pr_matched.append(pr_i)
			Qr_matched.append(qr_i)
		# finished
		return Pr_matched,Qr_matched
#----------------#
