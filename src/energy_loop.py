#!/usr/bin/python3

from outward import PQ_OUTWARD
from inward import R_INF,PQ_INWARD,ENERGY_CORRECTION
from physics import RESCALING_INWARD
from spline import NORM_INTEGERAL,PQ_NORMALIZED,EXPECTATION_R,\
					EXPECTATION_R_general,EXPECTATION_R_GENERAL_fx
from constants import c_au


class ENERGY:
	"""
	E_min = lower limit of Energy
	E_max = upper limit of energy
	E_val = value of the energy from outward/inward calculation
	Delta_E = correction of the energy
	"""
	def __init__(self,E_min, E_max,E_val, Delta_E):
		self.E_min = E_min
		self.E_max = E_max
		self.E_val = E_val
		self.Delta_E = Delta_E
	#
	def energy(self):
		E_new = self.E_val + self.Delta_E
		#
		EMIN_n = self.E_min
		EMAX_n = self.E_max
		#---===========----#
		if self.Delta_E > 0.0: # +ve DE
			EMAX = EMAX_n
			EMIN = self.E_val
			E    = E_new
		
		elif self.Delta_E < 0.0:  # -ve DE
			EMAX = self.E_val
			EMIN = EMIN_n
			E    = E_new
		else:
			EMAX = EMAX_n
			EMIN = EMIN_n
			E    = E_new
		## 
		if EMIN < -c_au*c_au:
			emin = -c_au*c_au 
			EMIN = emin 
		if EMAX > 0.0:
			EMAX = self.E_max
		if E > EMAX or E < EMIN:
			E = 0.5*(EMAX+EMIN)
		return E,EMIN,EMAX
#----------------#
class ENERGY_LOOP:
	"""
	This is the main loop cycle to get the orbital
	"""
	def __init__(self,energy,emin,emax,DE,eta_error,N_r,A,B,C,D,r,rv,\
	            np,l,kappa):
		self.E = energy
		self.EMIN = emin
		self.EMAX = emax
		self.DE = DE
		self.eta = eta_error
		self.N_rx = N_r
		self.A  = A
		self.B  = B
		self.C  = C
		self.D  = D
		self.R  = r
		self.RV = rv
		self.NP = np
		self.L = l
		self.K = kappa
	#
	def energy_loop(self):
		item = 1
		#
		DELL = abs(self.eta)
		#
		scf = 100  # always correct just number greater that error value
		#
		initial_energy = self.E
		Eva            = self.E
		#
		EMIN_loop = self.EMIN
		EMAX_loop = self.EMAX
		#
		E1 = self.E
		#
		stop = 0
		#
		output = []
		#
		DE = self.DE
		#
		# to stop the while if no change in consestency
		consestency_before = 0
		#print(" before loop scf=",scf,"DELL=",DELL," DE=",DE)
		while  scf > DELL and stop == 0:
			#
			print(" ===================== ")
			print(" iteration number ",item)
			eg = ENERGY(EMIN_loop,EMAX_loop,Eva,DE)
			Ex,EMIN_loop,EMAX_loop = eg.energy()
			# main loop
			# outward
			out_pq = PQ_OUTWARD(self.N_rx,self.A,self.B,self.C,self.D,\
					self.R,self.RV,self.NP,Ex,EMIN_loop,EMAX_loop,\
					self.L,self.K,self.eta,1)
			IOTP,Eva,Pr_outward,Qr_outward,P_out,Q_out,stab_ze,_=\
			         out_pq.pq_outward()
			#
			# Inward
			rf = R_INF(self.L,Eva,self.R,self.A,self.B,self.C,self.D,\
					IOTP,self.K)
			r_inf,P_r0,Q_r0 = rf.r_infinity()
			#
			in_pq = PQ_INWARD(r_inf,P_r0,Q_r0,self.R,self.A,self.B,\
					self.C,self.D,IOTP,self.K,self.eta,Eva)
			Pr_inward,Qr_inward,P_inw,Q_inw =in_pq.pq_inward()
			# combine two solutions
			rescale = RESCALING_INWARD(IOTP,Pr_outward,Qr_outward,\
					P_out,Pr_inward,Qr_inward,P_inw)
			Pr_matched,Qr_matched = rescale.pr_qr_matched()
			#
			nt = NORM_INTEGERAL(Pr_matched,Qr_matched,self.R)
			integt,R_new,Ax,Bx,Cx,Dx = nt.normalization_integeral()
			#
			pqnorm = PQ_NORMALIZED(Pr_matched,Qr_matched,integt)
			Pr_norm,Qr_norm = pqnorm.pq_normalized()
			# expectation value of orbital radius < r >
			# < r >
			r_expec  = EXPECTATION_R_GENERAL_fx(Pr_norm,Qr_norm,self.R)
			#
			expectation_r = r_expec.r_expectation()
			# <r^{2} >
			r_expec_2 = EXPECTATION_R_GENERAL_fx(Pr_norm,Qr_norm,\
						self.R,operator=lambda r: r**2)
			expectation_r2 = r_expec_2.r_expectation()
			# < 1/r >
			exp_inv_r = EXPECTATION_R_GENERAL_fx(Pr_norm,Qr_norm,\
						self.R,operator=lambda r: 1.0/r)
			expectation_inv_r = exp_inv_r.r_expectation()
			# correction energy
			den = ENERGY_CORRECTION(integt,P_out,Q_out,P_inw,Q_inw)
			DE = den.energy_correction()
			# reference consestency
			#
			E2 = Eva
			consestency =  abs(E2 - E1)
			#
			print(" DE (energy correction):",DE)
			print(" consestency =",consestency," eta (relative error):",DELL)
			#print(" E1",E1," E2=Eva", '{:<.30f}'.format(E2) )
			#print(" scf :",scf)
			#
			if consestency == consestency_before:
				stop = 1
			#
			if consestency == 0:
				Eva = E1
				stop = 1
			else:
				if consestency < DELL:
					scf = consestency
				if consestency > DELL:
					E1  = Eva
			#			
			# append and save all output needed
			outwer = ['{:<.30f}'.format(Eva),R_new,Pr_norm,Qr_norm,\
					  Pr_matched,Qr_matched,stab_ze,integt,\
					  expectation_r,expectation_r2,expectation_inv_r]
			output.append(outwer)
			#
			consestency_before = consestency
			item+=1
		#================#
		if stop == 1:
			if len(output) >= 2:
				return output[-2]
			else:
				return output[-1]
		else:
			return output[-1]
