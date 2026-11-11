#!/usr/bin/python3

from data_input import DATA

# data input for reading potential function

class POTENTIAL:
	"""
	This part to determine the potential used in the calculations
	
	
	"""
	def pot_fn(self):
		ln1 = " Choose which potential function used:\n"
		ln2 = " 1 - Thomas Fermi potential "
		ln3 = "   V(r) = -1/r [ Z phi_m - Z_{inf} (phi_m -1)]"
		ln4 = " 2 - Couloumb potential (Unscreened) "
		ln5 = "    V(r) = -Z/r"
		ln6 = " 3- Yukawa potential"
		ln7 = "    V(r) = -Z/r * e^{- alpha r}"
		ln8 = " 4- Multi Yukawa potential"
		ln9 = "    V(r) = -Z/r sum_{i=1}^{n_max} c_i e^{- alpha_i r}"
		ln10 = " 5- Screened potential" 
		ln11 = "    V(r) = -1/r [ 1+(Z-1) e^{-alpha r}] "
		ln12 = " 6- Screened-hydrogenic potential" 
		ln13 = "    V(r) = -Z_{eff}/r  = - (Z-alpha)/r "
		ln14 = " 7- Read from input file containing two columns of" 
		ln15 = "    r , v(r) in header and their values along the file" 
		##======##
		print(ln1);print(ln2);print(ln3);print(ln4)
		print(ln5);print(ln6);print(ln7);print(ln8)
		print(ln9);print(ln10);print(ln11);print(ln12)
		print(ln13);print(ln14);print(ln15)
		print(" 8- Uniformly charged sphere (finite nuclear size)")
		print("    rV(r) = -Z [ 3r/(2R) - r^3/(2R^3) ]  r<R ;  -Z  r>R")
		print("    R = sqrt(5/3) R_rms  (R_rms = rms nuclear radius)")
		print(" 9- Fermi charge distribution (finite nuclear size)")
		print("    rho(r) = rho_0 / ( 1 + exp[(r-c)/a] )")
		print("    a = t/(4 ln 3), t = 2.3 fm ; c fixed by R_rms")
		return
	# choose the potrential form
	def define_potential(self):
		#ref_list = [1,2,3,4,5,6]
		ref_list = [1,2,3,4,5,6,7,8,9]
		self.pot_fn()
		print("\n")
		print(" Enter the index/number of potential used:")
		pot_num = input()
		# check correction 
		d1 = pot_num.strip()
		#print(d1)
		#
		# check if it is number as integer or not 
		out = DATA().check_number_integer(d1)
		while out == False:
			m1 = " Error: Enter the index/number of potential used:"
			m2 = " Enter 1 or 2 or 3 or ...."
			print(m1);print(m2)
			pot_num = input()
			d1 = pot_num.strip()
			# check if it is number as integer or not 
			out = DATA().check_number_integer(d1)
			if out== True and int(d1) not in ref_list:
				out = False
		#
		"""
		while d1 == "" or int(d1) not in ref_list:
			m1 = " Error: Enter the index/number of potential used:"
			m2 = " Enter 1 or 2 or 3 or ...."
			print(m1);print(m2)
			pot_num = input()
			d1 = pot_num.strip()
		"""
		return int(d1)
		
		
	# define function for Yukawa alpha parameter
	def alpha(self):
		d1 = " Enter the alpha screening parameter"
		print(d1)
		alpha_p = input()
		# check if it is number or not 
		out = DATA().check_number(alpha_p)
		while out == False:
			d2 = " Error: Enter the alpha screening parameter"
			print(d2)
			alpha_p = input()
			# check if it is number or not 
			out = DATA().check_number(alpha_p)
		return float(alpha_p)
	#
	# define  n, c, alpha for multi Yukawa parameters
	def multi_yukawa_nmax(self):
		d1 = " Enter the maximum series number n_max as integer:"
		print(d1)
		nmax = input()
		# check if it is number as integer or not 
		out = DATA().check_number_integer(nmax)
		while out == False:
			d2 = " Error: Enter the maximum series number n_max as integer:"
			print(d2)
			nmax = input()
			# check if it is number as integer or not 
			out = DATA().check_number_integer(nmax)
		return int(nmax)
		
	def multi_yukawa_parameters(self):
		nmax = self.multi_yukawa_nmax()
		terms = []
		p1 = " Enter c"
		p2 = " Enter alpha"
		for t in range(1, int(nmax)+1):
			p11 = p1+str(t)
			print(p11)
			cp = input()
			# check if it is number or not 
			out = DATA().check_number(cp)
			while out == False:
				d2 = " Error: "+p11
				print(d2)
				cp = input()
				# check if it is number or not 
				out = DATA().check_number(cp)
			# 
			p21 = p2+str(t)
			print(p21)
			ap = input()
			# check if it is number or not 
			out1 = DATA().check_number(ap)
			while out1 == False:
				d3 = " Error: "+p21
				print(d3)
				ap = input()
				# check if it is number or not 
				out1 = DATA().check_number(ap)
			#
			
			mm = (float(cp.strip()) ,float(ap.strip()) )
			terms.append(mm)
		#
		return nmax,terms
	#
	#
	## read input potential file 
	def potent_file(self):
		import os
		p1 = " Enter the potential file with data r, V(r)"
		print(p1)
		pot_file =  input()
		#
		pot_file = pot_file.strip()
		#
		while pot_file== "" or not os.path.isfile(pot_file):
			d3 = " Error: "+p1
			print(d3)
			pot_file =  input()
			pot_file = pot_file.strip()
		#
		return pot_file
	#
	# ---- finite nuclear size (options 8 and 9) ----
	# rms nuclear charge radius in fm, either given directly or
	# estimated from the mass number A
	def nuclear_rms_radius(self):
		from potential import rms_radius_from_mass_number
		p1 = " Enter the nuclear mass number A to estimate the rms radius"
		p2 = " (R_rms = 0.836 A^(1/3) + 0.570 fm), or press Enter to"
		p3 = " give the rms radius R_rms (fm) directly:"
		print(p1)
		print(p2)
		print(p3)
		a_in = input().strip()
		while a_in != "" and (DATA().check_number_integer(a_in) == False \
				or int(a_in) <= 0):
			print(" Error: A must be a positive integer, or press Enter")
			print(p3)
			a_in = input().strip()
		if a_in != "":
			rrms = rms_radius_from_mass_number(int(a_in))
			print(" The estimated rms radius is ",rrms," fm")
			return rrms
		d1 = " Enter the rms nuclear charge radius R_rms in fm:"
		print(d1)
		rr = input().strip()
		while DATA().check_number(rr) == False or float(rr) <= 0.0:
			print(" Error: R_rms must be a positive number (fm)")
			print(d1)
			rr = input().strip()
		return float(rr)
	#
	# 10%-90% surface thickness of the Fermi distribution
	def fermi_skin_thickness(self):
		from constants import fermi_skin_t_fm
		d1 = " Enter the 10%-90% skin thickness t in fm"
		d2 = " (press Enter for the default "+str(fermi_skin_t_fm)+" fm):"
		print(d1)
		print(d2)
		tt = input().strip()
		while tt != "" and (DATA().check_number(tt) == False or \
				float(tt) <= 0.0):
			print(" Error: t must be a positive number (fm)")
			print(d2)
			tt = input().strip()
		if tt == "":
			return float(fermi_skin_t_fm)
		return float(tt)

###########################
