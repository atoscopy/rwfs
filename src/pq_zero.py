#!/usr/bin/python3

from physics import CONDN1, CONDN2, KAPPA_SIGN
from constants import c_au,half,end_n

class PQ_ZERO:
	"""
	a solution for the initial term series 
	
	# input
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
    ri,rj = are the values of r_i and r_j
            ri,rj == -1  if no need for GET_CONVERG 
            function or no need for grid subdivision
            any other values for use of GET_CONVERG subdivision fn.
	"""
	# (i,j,R,A,B,C,D,IOTP,K,SL,P_ra,Q_ra,EPS,E,ri,rj)
	def __init__(self,i,j,r,A,B,C,D,r_turn,kappa,P_ri,Q_ri,eta_error,E_init,ri,rj):
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
		self.ri 	= ri
		self.rj 	= rj
	#
	def pq_zero(self):
		# check the output convergence 
		# if occurs = 1  if not = 0
		convergence = 1 
		#
		a = self.i
		b = self.j
		if self.ri == self.rj == -1:
			ra = self.r[a]
			rb = self.r[b]
		else:
			ra = self.ri
			rb = self.rj
		#
		h = rb - ra
		h2 = h*h
		h3 = h2 * h
		#
		# nu values
		#index = self.i
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
		u1ch = u1 - 2*c_au*h
		kappa_sign = KAPPA_SIGN(self.K).kappa_sgn()
		ucr = ut-(2*c_au*rb)
		#
		# ONLY the first term for ra = 0 and rb 
		# i.e j = 1 , i = 0 so h = rb 
		if self.i == 0 and self.j == 1: 
			#
			a_coeff = []
			b_coeff = []
			#
			end = False
			#
			# ======= 1 ======
			# case i 
			#
			if u0 != 0:
				#
				k2 = self.K * self.K
				u02 = u0*u0
				#
				s = (k2-u02)**half
				t = 0
				#
				# a0 , b0
				a0 = 1.0
				a_coeff.append(a0)   
				b0 = -(s+self.K)/u0
				b_coeff.append(b0)
				#
				# the higher coefficients
				#
				# a1 , b1
				A1 = u1*a_coeff[0]     # = M1
				B1 = u1ch *b_coeff[0]  # = R1
				a1 =((-u0*A1)-((s+1-self.K)*B1))/((2*s)+1)
				a_coeff.append(a1) 
				b1 =(((s+1+self.K)*A1)-(u0*B1))/((2*s)+1)
				b_coeff.append(b1)
				#
				# a2 , b2
				A2 = (u1*a_coeff[1])+(u2*a_coeff[0])
				B2 = (u1ch*b_coeff[1])+(u2*b_coeff[0])
				a2 = (-(u0*A2)-((s+2-self.K)*B2))/(2*((2*s)+2))
				a_coeff.append(a2)
				b2 = (((s+2+self.K)*A2)-(u0*B2))/(2*((2*s)+2))
				b_coeff.append(b2)
				#
				# a3 , b3
				A3 = (u1*a_coeff[2])+(u2*a_coeff[1])+(u3*a_coeff[0])
				B3 = (u1ch*b_coeff[2])+(u2*b_coeff[1])+(u3*b_coeff[0])
				#
				a3 = (-(u0*A3)-((s+3-self.K)*B3))/(3*((2*s)+3))
				a_coeff.append(a3)
				b3 = (((s+3+self.K)*A3)-(u0*B3))/(3*((2*s)+3))
				b_coeff.append(b3)
				#
				# the P(x)=\phi(x) and Q(x)=psi(x) at r=rb or x=1 are
				Px1 = sum(a_coeff)
				Qx1 = sum(b_coeff)
				#
				# first derivatives of P(x)=phi(x) and Q(x)=psi(x) 
				# at r=rb or x=1 are
				d_Px1 = s*a_coeff[0]+(1+s)*a_coeff[1]+(2+s)*a_coeff[2]+\
						(3+s)*a_coeff[3]
				d_Qx1 = s*b_coeff[0]+(1+s)*b_coeff[1]+(2+s)*b_coeff[2]+\
						(3+s)*b_coeff[3]
				#
				# the higher coefficients with checking the conditions
				idx = 4
				while idx < end_n and end==False:   
					n = idx
					n1 = n-1
					n2 = n-2
					n3 = n-3
					A_n = u1*a_coeff[n1]+u2*a_coeff[n2]+u3*a_coeff[n3]
					B_n = u1ch*b_coeff[n1]+u2*b_coeff[n2]+u3*b_coeff[n3]
					#
					a_n = (-(u0*A_n)-((s+n-self.K)*B_n))/(n*((2*s)+n))
					a_coeff.append(a_n)
					b_n = (((s+n+self.K)*A_n)-(u0*B_n))/(n*((2*s)+n))
					b_coeff.append(b_n)
					#
					Px1 = Px1+a_n
					d_Px1 = d_Px1+((n+s)*a_n)
					#
					Qx1 = Qx1+b_n
					d_Qx1 = d_Qx1+((n+s)*b_n)
					#
					# check conditions
					#
					cond1 = CONDN1(a_n,b_n,Px1,Qx1,self.eta).cond1()
					cond2 = CONDN2(rb,self.K,ucr,ut,h,Px1,d_Px1,\
							Qx1,d_Qx1,self.eta).cond2()
					idx+=1
					if cond1 == True and cond2 == True:
						Q_r2 = Qx1
						P_r2 = Px1
						# now end the while loop
						end = True
				#
				if idx == end_n: # convergence not occur
					convergence=0
				#
			#
			# ======= 2 ======
			# case ii    (eqns. 28-30)
			#
			if u0==0 and kappa_sign=="-ve":
				#
				mod_kappa = abs(self.K)
				s = mod_kappa
				t = 1
				#
				# a0 , b0
				a0 = 1.0
				a_coeff.append(a0)   
				b0 = u1/((2*mod_kappa)+ 1)
				b_coeff.append(b0)
				# the higher coefficients
				#
				# a1 , b1
				a1 = 0.0
				a_coeff.append(a1) 
				b1 =((u1*a1)+(u2*a0))/((2*mod_kappa)+1+1)
				b_coeff.append(b1)
				#
				# a2 , b2
				a2 =-(u1ch*b0)/2
				a_coeff.append(a2) 
				b2 =((u1*a2)+(u2*a1)+(u3*a0) )/((2*mod_kappa)+2+1)
				b_coeff.append(b2)
				#
				# a3 , b3
				a3 =(-(u1ch*b1)-(u2*b0))/3
				a_coeff.append(a3) 
				b3 =((u1*a3)+(u2*a2)+(u3*a1))/((2*mod_kappa)+3+1)
				b_coeff.append(b3)
				#
				# a4 , b4
				a4 =(-(u1ch*b2)-(u2*b1)-(u3*b0))/4
				a_coeff.append(a4) 
				b4 =((u1*a4)+(u2*a3)+(u3*a2) )/((2*mod_kappa)+4+1)
				b_coeff.append(b4)
				#
				# the P(x)=phi(x) and Q(x)=psi(x) at r=rb or x=1 are
				Px1 = sum(a_coeff)
				Qx1 = sum(b_coeff)
				#
				# first derivatives of P(x) and Q(x) at r=rb or x=1 are
				d_Px1 = s*a_coeff[0]+(1+s)*a_coeff[1]+(2+s)*a_coeff[2]+\
						(3+s)*a_coeff[3]+(4+s)*a_coeff[4]
				d_Qx1 = (s+t)*b_coeff[0]+(1+s+t)*b_coeff[1]+(2+s+t)*b_coeff[2]+\
						(3+s+t)*b_coeff[3]+(4+s+t)*b_coeff[4]
				# the higher coefficients with checking the conditions
				idx = 5
				while idx< end_n and end==False:  
					n = idx
					n1 = n-1
					n2 = n-2
					n3 = n-3
					n4 = n-4
					a_n = (-(u1ch*b_coeff[n2])-(u2*b_coeff[n3])-(u3*b_coeff[n4]))/(n)
					a_coeff.append(a_n)
					b_n = ((u1*a_coeff[n])+(u2*a_coeff[n1])+(u3*a_coeff[n2]))/(2*mod_kappa+n+1)
					b_coeff.append(b_n)
					#
					Px1 = Px1+a_n
					d_Px1 = d_Px1+((n+s)*a_n)
					#
					Qx1 = Qx1+b_n
					d_Qx1 = d_Qx1+((n+s+t)*b_n)
					#
					# check conditions
					cond1 = CONDN1(a_n,b_n,Px1,Qx1,self.eta).cond1()
					cond2 = CONDN2(rb,self.K,ucr,ut,h,Px1,d_Px1,\
							Qx1,d_Qx1,self.eta).cond2()
					#
					idx+=1
					if cond1 == True and cond2 == True:
						Q_r2 = Qx1
						P_r2 = Px1
						# now end the loop
						end = True
				if idx == end_n:
					convergence=0
				#
			#
			# ======= 3 ======
			# case iii  
			#
			if u0==0 and kappa_sign=="+ve":
				#
				mod_kappa = abs(self.K)
				s = mod_kappa+1
				t = -1
				#
				# a0 , b0
				# add the condition of b0 = +1 or -1 so that
				# a0 is positive
				if u1ch > 0.0:
					b0 = -1.0
				else:
					b0 = 1.0
				b_coeff.append(b0)
				a0 = -u1ch*b0 /((2*mod_kappa)+1)
				a_coeff.append(a0) 
				# the higher coefficients
				#
				# a1 , b1
				b1 =0.0
				b_coeff.append(b1)
				a1 = (-(u1ch*b1)-(u2*b0))/((2*mod_kappa)+1+1)
				a_coeff.append(a1) 
				#
				# a2 , b2
				b2 =(u1*a0)/(2)
				b_coeff.append(b2)
				a2 =(-(u1ch*b2)-(u2*b1)-(u3*b0))/((2*mod_kappa)+2+1)
				a_coeff.append(a2) 
				#
				# a3 , b3
				b3 =((u1*a1)+(u2*a0))/(3)
				b_coeff.append(b3)
				a3 =(-(u1ch*b3)-(u2*b2)-(u3*b1))/((2*mod_kappa)+3+1)
				a_coeff.append(a3) 
				#
				# a4 , b4
				b4 = ((u1*a2)+(u2*a1)+(u3*a0))/(4)
				b_coeff.append(b4)
				a4 = (-(u1ch*b4)-(u2*b3)-(u3*b2))/((2*mod_kappa)+4+1)
				a_coeff.append(a4) 
				#
				# the P(x)=\phi(x) and Q(x)=\psi(x) at r=rb or x=1 are
				Px1 = sum(a_coeff)
				Qx1 = sum(b_coeff)
				# first derivatives of P(x) and Q(x) at r=rb or x=1 are
				d_Px1 = s*a_coeff[0]+(1+s)*a_coeff[1]+(2+s)*a_coeff[2]+\
						(3+s)*a_coeff[3]+(4+s)*a_coeff[4]
				d_Qx1 = (s+t)*b_coeff[0]+(1+s+t)*b_coeff[1]+(2+s+t)*b_coeff[2]+\
						(3+s+t)*b_coeff[3]+(4+s+t)*b_coeff[4]
				#
				# the higher coefficients with checking the conditions
				idx = 5
				while idx< end_n and end==False:   
					n = idx
					n1 = n-1
					n2 = n-2
					n3 = n-3
					n4 = n-4
					b_n = ((u1*a_coeff[n2])+(u2*a_coeff[n3])+(u3*a_coeff[n4]))/(n)
					b_coeff.append(b_n)
					a_n = (-(u1ch*b_coeff[n])-(u2*b_coeff[n1])-(u3*b_coeff[n2]))/\
						  ((2*mod_kappa)+n+1)
					a_coeff.append(a_n)
					#
					Px1 = Px1+a_n
					d_Px1 = d_Px1+((n+s)*a_n)
					#
					Qx1 = Qx1+b_n
					d_Qx1 = d_Qx1+((n+s+t)*b_n)
					#
					# check conditions
					cond1 = CONDN1(a_n,b_n,Px1,Qx1,self.eta).cond1()
					cond2 = CONDN2(rb,self.K,ucr,ut,h,Px1,d_Px1,\
							Qx1,d_Qx1,self.eta).cond2()
					#
					idx+=1
					if cond1 == True and cond2 == True:
						Q_r2 = Qx1
						P_r2 = Px1
						# now end the loop
						end = True
				if idx == end_n:
					convergence=0
				#
			#-------------------#
			# OUTPUT
			# collect all data for the zero point to get PZ coefficient
			zero_out = [s,t,a_coeff,b_coeff]  # stab
			# then get the P and Q for r=rb
			if convergence == 1: # convergence occurs
				return P_r2,Q_r2,convergence,zero_out
			if convergence == 0: # convergence NOT occurs
				return 0,0,convergence,zero_out
		# higher terms for ra != 0 and rb  i.e j > 1 , i != 0 
		else:
			print(" ERROR: It is not the first point (0-ra) or (0,r_j)")
			return 0,0,0,0
