#!/usr/bin/python3

from physics import KAPPA_SIGN,CONDN1,CONDN2
from constants import c_au,end_n


class PQ_MIDDLE:
	"""
	equation 12
	
	input:
	i = initial index of ra point = 0
    j = final index of rb point   = 1
    value of i,j |>> i == 0 and j == 1 <<|
     
    r = grid points
    A,B,C,D are the spline coefficients
    r_turn = end point of the grid (turning point)
    kappa = kappa of the orbital
    c_au = speed of light
    P_ri = values of P_ra at i points
    Q_ri = values of Q_ra at i points
    eta_error = tolerence of the error
    E_init = initial estimate of energy
	"""
	def __init__(self,i,j,r,A,B,C,D,r_turn,kappa,P_ri,Q_ri,eta_error,E_init):
		self.i 		= i
		self.j 		= j
		self.r 		= r
		self.A 		= A
		self.B 		= B
		self.C 		= C
		self.D 		= D
		self.IOTP 	= r_turn
		self.K 		= kappa
		self.P_ra 	= P_ri
		self.Q_ra 	= Q_ri 
		self.eta 	= eta_error
		self.E 		= E_init
	#
	def pq_middle(self):
		# check the output convergence 
		# if convergence occurs = 1  if not = 0
		convergence = 1 
		#
		a = self.i
		b = self.j
		ra = self.r[a]
		rb = self.r[b]
		h = rb - ra
		#
		h2 = h*h
		h3 = h2 * h
		#
		# nu values
		index = self.i
		if index >= len(self.A):
			index = len(self.A) - 1
		
		d0 = self.A[index]
		d1 = self.B[index]
		d2 = self.C[index]
		d3 = self.D[index]
		#
		ra2 = ra * ra
		ra3 = ra2 * ra
		d1E = d1 - self.E
		#
		u0 =  (d0 + (d1E*ra) + (d2*ra2) + (d3* ra3)) / c_au
		u1 =  ((d1E + (2*d2*ra) + (3*d3*ra2) )*h) / c_au
		u2 =  (( d2 + (3*d3 *ra) )*h2 )/ c_au
		u3 =  (d3 *h3) / c_au
		#
		ut = u0 + u1 + u2 + u3
		#
		# parameters used for all
		u1ch = u1 - (2*c_au*h)
		kappa_sign = KAPPA_SIGN(self.K).kappa_sgn()
		ucr = ut-(2*c_au*rb)
		# for middle region
		# i.e j > 1 , i != 0
		if self.i > 0 and self.j > 1: 
			# the middle region for calculation
			# with ra > 0
			# parameters used in common
			u0cra = u0-(2*c_au*ra)
			const1 = h/ra
			#
			a_coeff = []
			b_coeff = []
			#
			# total number of series 
			end = False
			#
			# a0 , b0
			a0 = self.P_ra
			a_coeff.append(a0) 
			b0 = self.Q_ra
			b_coeff.append(b0)
			#
			# the higher coefficients
			# a1 , b1
			a1 = -(const1/1.0)*( ((self.K+1-1)*a0) +(u0cra*b0) ) 
			a_coeff.append(a1) 
			b1 =(const1/1.0)*((u0*a0)+((self.K-1+1)*b0))
			b_coeff.append(b1)
			#
			# a2 , b2
			a2 = -(const1/2.0)*( ((self.K+2-1)*a1) +(u0cra*b1)+(u1ch*b0) ) 
			a_coeff.append(a2) 
			b2 =(const1/2.0)*((u0*a1)+((self.K-2+1)*b1)+(u1*a0))
			b_coeff.append(b2)
			#
			# a3 , b3
			a3 = -(const1/3.0)*( ((self.K+3-1)*a2)+(u0cra*b2)+(u1ch*b1)+(u2*b0) ) 
			a_coeff.append(a3) 
			b3 =(const1/3.0)*((u0*a2)+((self.K-3+1)*b2)+(u1*a1)+(u2*a0))
			b_coeff.append(b3)
			#
			# a4 , b4
			a4 = -(const1/4.0)*( ((self.K+4-1)*a3) +(u0cra*b3)+(u1ch*b2)+(u2*b1)+(u3*b0) ) 
			a_coeff.append(a4) 
			b4 =(const1/4.0)*((u0*a3)+((self.K-4+1)*b3)+(u1*a2)+(u2*a1)+(u3*a0))
			b_coeff.append(b4)
			#
			# the P(x)=\phi(x) and Q(x)=\psi(x) at r=rb or x=1 are
			Px1 = sum(a_coeff)
			Qx1 = sum(b_coeff)
			# first derivatives of P(x) and Q(x) at r=rb or x=1 are
			d_Px1 = (0)*a_coeff[0]+(1)*a_coeff[1]+(2)*a_coeff[2]+\
					(3)*a_coeff[3]+(4)*a_coeff[4]
			d_Qx1 = (0)*b_coeff[0]+(1)*b_coeff[1]+(2)*b_coeff[2]+\
					(3)*b_coeff[3]+(4)*b_coeff[4]
			#
			# the higher coefficients with checking the conditions        
			idx = 5
			while idx< end_n and end==False:   
				n = idx
				n1 = n-1
				n2 = n-2
				n3 = n-3
				n4 = n-4
				#
				a_n = -(const1/n)*( ((self.K+n-1)*a_coeff[n1]) +\
								(u0cra*b_coeff[n1])+\
								(u1ch*b_coeff[n2])+(u2*b_coeff[n3])+\
								(u3*b_coeff[n4]) ) 
				a_coeff.append(a_n)
				b_n = (const1/n)*( (u0*a_coeff[n1])+\
								((self.K-n+1)*b_coeff[n1])+\
								(u1*a_coeff[n2])+(u2*a_coeff[n3])+\
								(u3*a_coeff[n4]) )
				b_coeff.append(b_n)
				#
				Px1 = Px1+a_n
				d_Px1 = d_Px1+((n)*a_n)
				#
				Qx1 = Qx1+b_n
				d_Qx1 = d_Qx1+((n)*b_n)
				#
				# check conditions
				cond1 = CONDN1(a_n,b_n,Px1,Qx1,self.eta).cond1()
				cond2 = CONDN2(rb,self.K,ucr,ut,h,Px1,d_Px1,\
							Qx1,d_Qx1,self.eta).cond2()
				#
				idx+=1
				if cond1 == True and cond2 == True:
					P_r2 = Px1
					Q_r2 = Qx1
					#
					# now end the loop
					end = True
			#
			if idx == end_n:
				convergence=0
			#
			# then get the P and Q for r=rb
			if convergence == 1: # convergence occurs
				return P_r2,Q_r2,convergence
			if convergence == 0: # convergence NOT occurs
				return 0,0,convergence
		# If it is the initial or zero grid points
		else:
			print(" ERROR: (pq_middle) It is not the middle region but (0-ra) or (0,r_j)")
			return 0,0,0

#----------------#
class PQ_MIDDLE2:
	"""
	Only for not converging cases
	
	equation 12
	
	input:
	i = initial index of ra point = 0
    j = final index of rb point   = 1
    value of i,j |>> i == 0 and j == 1 <<|
     
    r = grid points
    A,B,C,D are the spline coefficients
    r_turn = end point of the grid (turning point)
    kappa = kappa of the orbital
    c_au = speed of light
    P_ri = values of P_ra at i points
    Q_ri = values of Q_ra at i points
    eta_error = tolerence of the error
    E_init = initial estimate of energy
	"""
	def __init__(self,i,j,r,A,B,C,D,r_turn,kappa,P_ri,Q_ri,eta_error,E_init):
		self.i 		= i
		self.j 		= j
		self.r 		= r
		self.A 		= A
		self.B 		= B
		self.C 		= C
		self.D 		= D
		self.IOTP 	= r_turn
		self.K 		= kappa
		self.P_ra 	= P_ri
		self.Q_ra 	= Q_ri 
		self.eta 	= eta_error
		self.E 		= E_init
	#
	def pq_middle(self):
		# check the output convergence 
		# if convergence occurs = 1  if not = 0
		convergence = 1 
		#
		a = self.i
		b = self.j
		ra = self.r[a]
		rb = self.r[b]
		h = rb - ra
		#
		h2 = h*h
		h3 = h2 * h
		#
		# nu values
		index = self.i
		if index >= len(self.A):
			index = len(self.A) - 1
		
		d0 = self.A[index]
		d1 = self.B[index]
		d2 = self.C[index]
		d3 = self.D[index]
		#
		ra2 = ra * ra
		ra3 = ra2 * ra
		d1E = d1 - self.E
		#
		u0 =  (d0 + (d1E*ra) + (d2*ra2) + (d3* ra3)) / c_au
		u1 =  ((d1E + (2*d2*ra) + (3*d3*ra2) )*h) / c_au
		u2 =  (( d2 + (3*d3 *ra) )*h2 )/ c_au
		u3 =  (d3 *h3) / c_au
		#
		ut = u0 + u1 + u2 + u3
		#
		# parameters used for all
		u1ch = u1 - (2*c_au*h)
		kappa_sign = KAPPA_SIGN(self.K).kappa_sgn()
		ucr = ut-(2*c_au*rb)
		# for middle region
		# i.e j > 1 , i != 0
		if self.i > 0 and self.j > 1: 
			# the middle region for calculation
			# with ra > 0
			# parameters used in common
			u0cra = u0-(2*c_au*ra)
			const1 = h/ra
			#
			a_coeff = []
			b_coeff = []
			#
			# total number of series
			end = False
			#
			# a0 , b0
			a0 = self.P_ra
			a_coeff.append(a0) 
			b0 = self.Q_ra
			b_coeff.append(b0)
			#
			# the higher coefficients
			#
			# a1 , b1
			a1 = -(const1/1.0)*( ((self.K+1-1)*a0) +(u0cra*b0) ) 
			a_coeff.append(a1) 
			b1 =(const1/1.0)*((u0*a0)+((self.K-1+1)*b0))
			b_coeff.append(b1)
			#
			# a2 , b2
			a2 = -(const1/2.0)*( ((self.K+2-1)*a1) +(u0cra*b1)+(u1ch*b0) ) 
			a_coeff.append(a2) 
			b2 =(const1/2.0)*((u0*a1)+((self.K-2+1)*b1)+(u1*a0))
			b_coeff.append(b2)
			#
			# a3 , b3
			a3 = -(const1/3.0)*( ((self.K+3-1)*a2)+(u0cra*b2)+(u1ch*b1)+(u2*b0) ) 
			a_coeff.append(a3) 
			b3 =(const1/3.0)*((u0*a2)+((self.K-3+1)*b2)+(u1*a1)+(u2*a0))
			b_coeff.append(b3)
			#
			# a4 , b4
			a4 = -(const1/4.0)*( ((self.K+4-1)*a3) +(u0cra*b3)+(u1ch*b2)+(u2*b1)+(u3*b0) ) 
			a_coeff.append(a4) 
			b4 =(const1/4.0)*((u0*a3)+((self.K-4+1)*b3)+(u1*a2)+(u2*a1)+(u3*a0))
			b_coeff.append(b4)
			#
			# the P(x)=\phi(x) and Q(x)=\psi(x) at r=rb or x=1 are
			Px1 = sum(a_coeff)
			Qx1 = sum(b_coeff)
			# first derivatives of P(x) and Q(x) at r=rb or x=1 are
			d_Px1 = (0)*a_coeff[0]+(1)*a_coeff[1]+(2)*a_coeff[2]+\
					(3)*a_coeff[3]+(4)*a_coeff[4]
			d_Qx1 = (0)*b_coeff[0]+(1)*b_coeff[1]+(2)*b_coeff[2]+\
					(3)*b_coeff[3]+(4)*b_coeff[4]
			#
			# the higher coefficients with checking the conditions        
			idx = 5
			while idx< end_n and end==False:   
				n = idx
				n1 = n-1
				n2 = n-2
				n3 = n-3
				n4 = n-4
				#
				a_n = -(const1/n)*( ((self.K+n-1)*a_coeff[n1]) +\
								(u0cra*b_coeff[n1])+\
								(u1ch*b_coeff[n2])+(u2*b_coeff[n3])+\
								(u3*b_coeff[n4]) ) 
				a_coeff.append(a_n)
				b_n = (const1/n)*( (u0*a_coeff[n1])+\
								((self.K-n+1)*b_coeff[n1])+\
								(u1*a_coeff[n2])+(u2*a_coeff[n3])+\
								(u3*a_coeff[n4]) )
				b_coeff.append(b_n)
				#
				Px1 = Px1+a_n
				d_Px1 = d_Px1+((n)*a_n)
				#
				Qx1 = Qx1+b_n
				d_Qx1 = d_Qx1+((n)*b_n)
				#
				# check conditions
				cond1 = CONDN1(a_n,b_n,Px1,Qx1,self.eta).cond1()
				cond2 = CONDN2(rb,self.K,ucr,ut,h,Px1,d_Px1,\
							Qx1,d_Qx1,self.eta).cond2()
				#
				idx+=1
				# change condition than the other function
				if cond1 == True or cond2 == True:
					P_r2 = Px1
					Q_r2 = Qx1
					#
					# now end the loop
					end = True
			#
			if idx == end_n:
				convergence=0
			#
			# then get the P and Q for r=rb
			if convergence == 1: # convergence occurs
				return P_r2,Q_r2,convergence
			if convergence == 0: # convergence NOT occurs
				return 0,0,convergence
		# If it is the initial or zero grid points
		else:
			print(" ERROR: (pq_middle) It is not the middle region but (0-ra) or (0,r_j)")
			return 0,0,0
#----------------#
