#!/usr/bin/python3


from pq_zero import PQ_ZERO
from pq_middle import PQ_MIDDLE,PQ_MIDDLE2
from constants import c_au,loops,print_error_messages


class CONVEGE:
	"""
	check the convergence of the power series solution
	
		# the ra and rb interval should be divided into several points 
		# to get convergence

		# 1- Input
		# ra , P(ra),Q(ra) 
		# rb 
		# i,j index

		# 2- Output
		# numbers of zero points in between
		# P(rb),Q(rb)
		
		where A,B,C,D = spline coefficients
	"""
	def __init__(self,i,j,ri,rj,P_ri,Q_ri,A,B,C,D,r,r_turn,kappa,eta_error,E_init):
		self.i 		= i
		self.j 		= j
		self.ra 	= ri
		self.rb 	= rj
		self.P_ra 	= P_ri
		self.Q_ra 	= Q_ri 
		self.A 		= A
		self.B 		= B
		self.C 		= C
		self.D 		= D
		self.r 		= r
		self.IOTP 	= r_turn
		self.K 		= kappa
		self.eta 	= eta_error
		self.E 		= E_init
	#
	
	def get_converg(self):
		# force the convergence to occur step by step
		tot_convergence = 0
		#
		P_out = [self.P_ra]
		Q_out = [self.Q_ra]
		#
		h = self.rb - self.ra
		h_new = h
		# 
		# spline coefficients
		#
		ra_x     = []
		rb_x     = []
		h_x      = []
		con_verg = []
		#
		ra_n = self.ra
		rb_end = self.rb
		#
		# solve for the new division
		#
		end = 0
		#
		h_new = h_new/2.0 
		nsteps = 1
		#
		rx = -1
		ry = -1
		#
		#loops = 600
		while nsteps < loops and end == 0 :
			#
			ra_n = ra_n
			rb_n = ra_n + h_new
			#
			P_ra = P_out[-1]
			Q_ra = Q_out[-1]
			#
			rx = ra_n
			ry = rb_n
			#
			if self.i == 0: # zero first step
				pq = PQ_ZERO(self.i,self.j,self.r,self.A,self.B,self.C,\
						self.D,self.IOTP,self.K,self.P_ra,\
						self.Q_ra,self.eta,self.E,rx,ry)
				P1,Q1,conv,_stab = pq.pq_zero()
			else:
				pq = PQ_MIDDLE(self.i,self.j,self.r,self.A,self.B,\
						self.C,self.D,self.IOTP,self.K,self.P_ra,\
						self.Q_ra,self.eta,self.E)
				P1,Q1,conv = pq.pq_middle()
			#
			if conv == 0:
				h_new = h_new/2.0 
				#
				ra_n = ra_n
				rb_n = ra_n + h_new
			#
			# check end or not
			#if conv == 1 and rb_n == rb_end:
			if conv == 1 and abs(rb_n - rb_end) <= 1e-14 * max(abs(rb_end), 1.0):
				end = 1
				tot_convergence = 1
			#
			if conv == 1: # get the other side  
				# save data first
				ra_x.append(ra_n)
				rb_x.append(rb_n)
				h_x.append(h_new)
				con_verg.append( conv )
				P_out.append(P1)
				Q_out.append(Q1)
				#
				ra_n = rb_n
				rb_n = rb_end
				#
				h_new = rb_n - ra_n			
			#
			nsteps+=1
			#
		#check the output
		zeros = 0
		if end == 1 :
			for ix in range(len(P_out)-1):
				prev = ix
				now  = ix + 1
				if P_out[now]*P_out[prev] < 0.0:
					# it means if the sign is changed
					zeros+=1
			# end
			if zeros > 0:
				if print_error_messages == True:
					print(" The number of nodes are ",zeros,\
						  " between points ri:",self.ra," rj:",self.rb)
					print(" The grid should be denser")
			#
			P_r2 = P_out[-1]
			Q_r2 = Q_out[-1]
			#
			return P_r2,Q_r2,tot_convergence
		#
		if end ==0 and nsteps == loops:
			if print_error_messages == True:
				#print(" NOTE: (get_converg) The grid should be denser")
				#print(" Or use number of loops > ",loops)
				#
				#print(" use different grid type or ")
				#print(" denser grid or use number of loops > ",loops)
				print(" use different grid/potential/loops > ",loops)
			return 0,0,tot_convergence

#
#

class CONVEGE2:
	"""
	check the convergence of the power series solution
	
		# the ra and rb interval should be divided into several points 
		# to get convergence

		# 1- Input
		# ra , P(ra),Q(ra) 
		# rb 
		# i,j index

		# 2- Output
		# numbers of zero points in between
		# P(rb),Q(rb)
		
		where A,B,C,D = spline coefficients
	"""
	def __init__(self,i,j,ri,rj,P_ri,Q_ri,A,B,C,D,r,r_turn,kappa,eta_error,E_init):
		self.i 		= i
		self.j 		= j
		self.ra 	= ri
		self.rb 	= rj
		self.P_ra 	= P_ri
		self.Q_ra 	= Q_ri 
		self.A 		= A
		self.B 		= B
		self.C 		= C
		self.D 		= D
		self.r 		= r
		self.IOTP 	= r_turn
		self.K 		= kappa
		self.eta 	= eta_error
		self.E 		= E_init
	#
	
	def get_converg(self):
		# force the convergence to occur step by step
		tot_convergence = 0
		#
		P_out = [self.P_ra]
		Q_out = [self.Q_ra]
		#
		h = self.rb - self.ra
		h_new = h
		# 
		# spline coefficients
		#
		ra_x     = []
		rb_x     = []
		h_x      = []
		con_verg = []
		#
		ra_n = self.ra
		rb_end = self.rb
		#
		# solve for the new division
		#
		end = 0
		#
		h_new = h_new/2.0 
		nsteps = 1
		#
		rx = -1
		ry = -1
		#
		while nsteps < loops and end == 0 :
			#
			ra_n = ra_n
			rb_n = ra_n + h_new
			#
			P_ra = P_out[-1]
			Q_ra = Q_out[-1]
			#
			rx = ra_n
			ry = rb_n
			#
			if self.i == 0: # zero first step
				pq = PQ_ZERO(self.i,self.j,self.r,self.A,self.B,self.C,\
						self.D,self.IOTP,self.K,self.P_ra,\
						self.Q_ra,self.eta,self.E,rx,ry)
				P1,Q1,conv,_stab = pq.pq_zero()
			else:
				pq = PQ_MIDDLE2(self.i,self.j,self.r,self.A,self.B,\
						self.C,self.D,self.IOTP,self.K,self.P_ra,\
						self.Q_ra,self.eta,self.E)
				P1,Q1,conv = pq.pq_middle()
			#
			if conv == 0:
				h_new = h_new/2.0 
				#
				ra_n = ra_n
				rb_n = ra_n + h_new
			#
			# check end or not
			if conv == 1 and rb_n == rb_end:
				end = 1
				tot_convergence = 1
			#
			if conv == 1: # get the other side  
				# save data first
				ra_x.append(ra_n)
				rb_x.append(rb_n)
				h_x.append(h_new)
				con_verg.append( conv )
				P_out.append(P1)
				Q_out.append(Q1)
				#
				ra_n = rb_n
				rb_n = rb_end
				#
				h_new = rb_n - ra_n			
			#
			nsteps+=1
			#
		#check the output
		zeros = 0
		if end == 1 :
			for ix in range(len(P_out)-1):
				prev = ix
				now  = ix + 1
				if P_out[now]*P_out[prev] < 0.0:
					# it means if the sign is changed
					zeros+=1
			# end
			if zeros > 0:
				if print_error_messages == True:
					print(" The number of nodes are ",zeros,\
						  " between points ri:",self.ra," rj:",self.rb)
					print(" The grid should be denser")
			#
			P_r2 = P_out[-1]
			Q_r2 = Q_out[-1]
			#
			return P_r2,Q_r2,tot_convergence
		#
		if end ==0 and nsteps == loops:
			if print_error_messages == True:
				print(" NOTE: (get_converg) The grid should be denser")
				print(" Or use number of loops > ",loops)
			return 0,0,tot_convergence
