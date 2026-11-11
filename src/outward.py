#!/usr/bin/python3


from pq_zero import PQ_ZERO
from convergence import CONVEGE,CONVEGE2
from pq_middle import PQ_MIDDLE
from physics import CONDN1,CONDN2,KAPPA_SIGN
from constants import outward_cycles,print_error_messages

# class 1
class TURNING_POINT:
	"""
	define the turning point of each orbital
	such that
	
	V(r_j) > E - l(l+1)/2r^2_j     if j>= t
	
	input:
	l = angular momentum
    r = grid points
    rv = rV(r)
    np = node points
    E = initial energy 
	"""
	def __init__(self, l, r, rv, np, E_init):
		self.l 		= l
		self.r 	    = r
		self.rv 	= rv
		self.np     = np
		self.E      = E_init
	#
	def turn_point(self):
		fl = 0.5 * self.l * ( self.l + 1 )   
		#
		idx = []
		for s in range(self.np):
			j = self.np - s - 1
			if j == 0:
				pass
			else:
				r2 = self.r[j]*self.r[j]
				term1 = fl/r2 
				value = self.E-(self.rv[j]/self.r[j])-(fl/r2) 
				if value > 0.0 :
					idx.append(j)
		#-------#
		# modified if the grid points are good
		if len(idx)>0:
			turn_point = idx[0]
		else:
			turn_point = self.np-1
		# 
		if turn_point > self.np-1:
			print(" ERROR: error increase the grid")
			return 0
		else:
			return turn_point
#-----------------------#
# class 2
class PQ_OUTWARD:
	"""
	the main class for outward solution
	# main outward solution

	N_r = actual nodes
	A,B,C,D = spline coefficients
	R = grid points
	RV = r V(r)
	NP = total node points
	E = initial estimated energy
	EMIN = minimum energy used
	EMAX = maximum energy used
	L = orbital angular momentum
	K = kappa
	SL = speed of light
	EPS = estimated error
	prt = print parameter takes 
	   0 for initial loop ,1 for main loop
	
	"""
	def __init__(self,nr,A,B,C,D,r,rv,np,E_init,emin,emax,l,kappa,\
			eta_error,prt):
		self.N_r = nr
		self.A  = A
		self.B  = B
		self.C  = C
		self.D  = D
		self.R  = r
		self.RV = rv
		self.NP = np
		self.E = E_init
		self.EMIN = emin
		self.EMAX = emax
		self.L = l
		self.K = kappa
		self.eta = eta_error
		self.prt = prt
	#
	def pq_outward(self):
		# maximum iteration is outward_cycles = 1000
		# is determined in constants file
		scf = 0
		#
		nodesr = -1   # number of nodes to enter while loop
		#
		avoid_zero = True 
		countzz = 0
		#
		while nodesr != self.N_r and scf < outward_cycles:
			# 
			r_t = TURNING_POINT(self.L,self.R,self.RV,self.NP,self.E)
			IOTP = r_t.turn_point()
			#
			En= self.E
			nodesr = 0
			#  ****  OUTWARD SOLUTION.  
			# initial trerms of P and Q are 0.0 
			Pr_outward = [0.0]
			Qr_outward = [0.0]
			# 
			stab = []  # initialize stab in case IOTP==0
			#
			#print(" IOTP = ",IOTP)
			#print(" --> ",len(self.R))
			# avoid that IOTP is the last point in the grid
			if IOTP == len(self.R)-1:
				# avoid this and get error message
				if print_error_messages == True:
					print(" Error: increase the grid/potential/loops")
				break
			# main outward calculations
			else: 
				for k in range(IOTP):
					#print(" k count = ",k)
					# first grid point
					if k==0:   
						i = k
						j = k+1
						P_ra = Pr_outward[i]
						Q_ra = Qr_outward[i]
						#
						pqz = PQ_ZERO(i,j,self.R,self.A,self.B,\
							self.C,self.D,IOTP,self.K,P_ra,\
							Q_ra,self.eta,self.E,-1,-1)
						P_rb,Q_rb,converg,stab = pqz.pq_zero()
						# no node count is used for the first point
						# check convergence
						if converg == 0:
							ran = self.R[i]
							rbn = self.R[j]
							#
							converg_fn = CONVEGE(i,j,ran,rbn,P_ra,Q_ra,\
							self.A,self.B,self.C,self.D,self.R,IOTP,self.K,\
							self.eta,self.E)
							P_rb,Q_rb,convergence = converg_fn.get_converg()
							#
							if convergence == 0:
								print(" STOP() CHECK: (get_converg) in\
								 PQ_OUTWARD fn.")
								print(" for zero point ")
								print(" ra = ",ran," rb = ",rbn)
								break
						# normalization
						p_norm = P_rb / abs(P_rb)
						q_norm = Q_rb / abs(P_rb)
						#
						Pr_outward.append(p_norm)
						Qr_outward.append(q_norm)
						#
					# the grid points rather than the first point
					else:
						#
						initial = k   #  = i
						final = k+1     #  = j
						P_ra = Pr_outward[initial]
						Q_ra = Qr_outward[initial]
						#
						pq = PQ_MIDDLE(initial,final,self.R,self.A,self.B,\
							self.C,self.D,IOTP,self.K,P_ra,\
							Q_ra,self.eta,self.E)
						P_rb,Q_rb,converg = pq.pq_middle()
						# check convergence
						if converg == 0:
							ran = self.R[initial]
							rbn = self.R[final]
							#
							converg_fn = CONVEGE(initial,final,ran,rbn,\
								P_ra,Q_ra,self.A,self.B,self.C,self.D,\
								self.R,IOTP,self.K,self.eta,self.E)
							P_rb,Q_rb,convergence = converg_fn.get_converg()
							#
							if convergence == 0:
								#print(" -> try CONVERG2")
								# 
								ran = self.R[initial]
								rbn = self.R[final]
								#
								converg_fn = CONVEGE2(initial,final,ran,\
								rbn,P_ra,Q_ra,self.A,self.B,self.C,self.D,\
									self.R,IOTP,self.K,self.eta,self.E)
								P_rb,Q_rb,convergence =\
									   converg_fn.get_converg()
								if convergence == 0:
									# NOW we use CONVEGE2
									print(" STOP() CHECK: (CONVEGE2) in\
									 PQ_OUTWARD fn.")
									print(" for middle points with")
									print(" ra = ",ran," rb = ",rbn)
									break
								#
						# count the node points 
						P_prev = Pr_outward[k]
						P_now = P_rb
						if P_prev*P_now < 0.0: 
							# it means if the sign is changed
							nodesr+=1
						#
						Pr_outward.append(P_rb)
						Qr_outward.append(Q_rb)
				#------------------------------#   
				scf+=1
				#print("\n")
				if self.prt == 0 :
					print(" cycle number ",scf)
				else:
					pass
				#
				print(" initial energy E = ",self.E)
				if nodesr == self.N_r:
					print(" the orbital has the correct count of nodes ",\
							nodesr )
				else:
					print(" our orbital have ",nodesr,\
					" while the correct nodes are :",self.N_r)
				#
				# Change the energy    
				# C  ****  TOO MANY NODES.
				if nodesr > self.N_r:
					#
					self.EMAX = self.E
					self.E = 0.5*(self.EMIN+self.E)

				#C  ****  TOO FEW NODES.
				# increase energy 
				elif nodesr < self.N_r:
					self.EMIN = self.E
					self.E = 0.5*(self.E+self.EMAX)
				else:
					self.E = self.E
				#
				print(" the next new energy is E = ",self.E)
				if self.prt == 0 :
					print("\n")
				else:
					pass
		if nodesr < 0:
			return IOTP,self.E,Pr_outward,Qr_outward,P_out,Q_out,stab,nodesr
		else:
			# last outward values
			P_out = Pr_outward[-1]
			Q_out = Qr_outward[-1]
			return IOTP,self.E,Pr_outward,Qr_outward,P_out,Q_out,stab,nodesr
#-----------------------#
#
class PZ:
	"""
	determine the first term in the series of the PQ_ZERO function 
	as in GRASP programe 
	see: 
	"Foundations of the Relativistic Theory of Atomic and Molecular
	 Structure" by Ian P. Grant, Harry M. Quiney
	 https://doi.org/10.1016/S0065-2199(08)60105-0
	"""
	def __init__(self,normx,a_node,rbx,sx):
		self.normx   = normx
		self.a_node  = a_node
		self.rbx 	 = rbx
		self.sx 	 = sx
	#
	def pz_parameter(self):
		factor = 1.0 / ( self.normx**(1.0/2.0) )
		pzz = (self.a_node / (self.rbx**self.sx) ) * factor
		return pzz

#-----------------------#
