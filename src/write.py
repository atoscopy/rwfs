#!/usr/bin/python3




class WRITE_DATA:
	def __init__(self,file_name,n,kappa,E,nnp,r,pr,qr,zz,nn):
		self.f = file_name
		self.n = n
		self.k = kappa
		self.E = E
		self.np = nnp
		self.r = r
		self.pr = pr
		self.qr = qr
		self.zz = zz  # atomic number
		self.nn = nn  # number of electrons
	#
	def write(self):
		f = open(self.f,"w")
		f.write("\n")   # new empty line
		# print atomic number and electron number
		lk1 = "atomic number (Z)        = "+'{:<8d}'.format(self.zz)
		lk2 = "number of electrons (Ne) = "+'{:<8d}'.format(self.nn)
		f.write(lk1+"\n")
		f.write(lk2+"\n")
		f.write("\n")   # new empty line
		#
		ln0 = "       n     kappa       Energy(a.u)         NP"
		f.write(ln0+"\n")
		line1 = '{:>8d}'.format(self.n)+'{:>8d}'.format(self.k)+\
		'{:>25.16E}'.format(self.E)+'{:>8d}'.format(self.np)
		f.write(line1+"\n")
		f.write("\n")   # new empty line
		lne1  = "          r                       P(r)"
		lne2 = "                     Q(r)"
		f.write(lne1+lne2+"\n")
		for i in range(self.np):
			 lnx = '{:<25.16E}'.format(self.r[i])+\
			 '{:<25.16E}'.format(self.pr[i])+'{:<25.16E}'.format(self.qr[i])
			 f.write(lnx+"\n")
		f.write("\n")   # new empty line
		#
		f.close()
#----------------#
class WRITE_DATA_MULTI:
	"""
	write MULTIPLE orbitals (radial wave functions) into a single file.
	the atomic number / number of electrons header is written once at
	the top of the file, followed by one self-contained block per
	orbital (n, kappa, energy, NP, then the r,P(r),Q(r) table), in the
	same layout WRITE_DATA uses for a single orbital.
	input:
	 file_name     = output file name
	 zz            = atomic number
	 nn            = number of electrons
	 orbitals_data = list of tuples, one per orbital, each of the form
	                 (n, kappa, E, np, r, pr, qr)
	"""
	def __init__(self,file_name,zz,nn,orbitals_data):
		self.f = file_name
		self.zz = zz  # atomic number
		self.nn = nn  # number of electrons
		self.orbitals_data = orbitals_data
	#
	def write(self):
		f = open(self.f,"w")
		f.write("\n")   # new empty line
		# print atomic number and electron number (once for the file)
		lk1 = "atomic number (Z)        = "+'{:<8d}'.format(self.zz)
		lk2 = "number of electrons (Ne) = "+'{:<8d}'.format(self.nn)
		f.write(lk1+"\n")
		f.write(lk2+"\n")
		lk3 = "number of orbitals       = "+'{:<8d}'.format(len(self.orbitals_data))
		f.write(lk3+"\n")
		f.write("\n")   # new empty line
		#
		for orb in self.orbitals_data:
			n,kappa,E,nnp,r,pr,qr = orb
			#
			ln0 = "       n     kappa       Energy(a.u)         NP"
			f.write(ln0+"\n")
			line1 = '{:>8d}'.format(n)+'{:>8d}'.format(kappa)+\
			'{:>25.16E}'.format(E)+'{:>8d}'.format(nnp)
			f.write(line1+"\n")
			f.write("\n")   # new empty line
			lne1  = "          r                       P(r)"
			lne2 = "                     Q(r)"
			f.write(lne1+lne2+"\n")
			for i in range(nnp):
				 lnx = '{:<25.16E}'.format(r[i])+\
				 '{:<25.16E}'.format(pr[i])+'{:<25.16E}'.format(qr[i])
				 f.write(lnx+"\n")
			f.write("\n")   # new empty line separating orbital blocks
		#
		f.close()
#----------------#
class WRITE_DATA_GRASP:
	def __init__(self,file_name,n,kappa,E,nnp,pz,r,pr,qr):
		self.f = file_name
		self.n = n
		self.k = kappa
		self.E = E
		self.np = nnp
		self.pz = pz
		self.r = r
		self.pr = pr
		self.qr = qr
	#
	def REFINE(self,x,index):
		"""
		index = 0   for energy eigenvalue
		index = 1 for P(r), Q(r) or r values
		"""
		x2 = x.split('E')
		power = x2[1]
		left = x2[0].split('.')
		# check if it is positive or negative
		if '-' in left[0]:
			if index == 0: # for eigenvalues
				my_left = ' 0.'
			else:
				my_left = '-0.'
		else:
			my_left = ' 0.'
		my_num = str( abs(int(left[0])) )+ left[1]
		power = int(x2[1]) + 1
		if -10 < power < 10:
			if '-' in str(power):
				f_power = '-0'+ str(abs(power))
			else:
				f_power = '+0'+str(abs(power))
		else:
			f_power = str(power)
		out = my_left+my_num +'E'+f_power
		# modify values of zero
		out_grasp = my_left+my_num +'D'+f_power
		# zero values
		if abs(float(out)) == 0 : 
			return "  "+my_left+my_num +'D+00' 
		# other values than zero
		else:
			return "  "+out_grasp
	#
	def write_grasp(self):
		f = open(self.f,"w")
		f.write("G92RWF"+"\n")   # new empty line
		line1 = '{:>8d}'.format(self.n)+'{:>8d}'.format(self.k)+\
		self.REFINE( '{:>25.15E}'.format(self.E), 0)+\
		'{:>8d}'.format(self.np)
		f.write(line1+"\n")
		# pz
		lne1 = self.REFINE('{:>25.15E}'.format(self.pz),1)
		#lne1 = '{:>25.16E}'.format(self.pz)
		f.write(lne1+"\n")
		# P(r)
		for i in range(self.np):
			 lnx1 = self.REFINE( '{:>25.15E}'.format(self.pr[i]),1)
			 f.write(lnx1+"\n")
		# Q(r)
		for i in range(self.np):
			 lnx2 = self.REFINE( '{:>25.15E}'.format(self.qr[i]),1)
			 f.write(lnx2+"\n")
		# r
		for i in range(self.np):
			 lnx3 = self.REFINE( '{:>25.15E}'.format(self.r[i]),1)
			 f.write(lnx3+"\n")
		#
		f.close()
#----------------#
class WRITE_DATA_GRASP_MULTI:
	"""
	write MULTIPLE orbitals into a single file in GRASP rwfn format.
	the "G92RWF" magic header is written ONCE at the top of the file,
	followed by one record block per orbital -- exactly the layout
	real multi-orbital GRASP .rwf files use, and readable by GRASP
	tools that expect several orbitals concatenated in one file.
	each orbital's block is: header line (n,kappa,E,np), pz, then the
	P(r), Q(r), r arrays in turn -- identical to WRITE_DATA_GRASP's
	single-orbital layout, just repeated per orbital.
	input:
	 file_name     = output file name
	 orbitals_data = list of tuples, one per orbital, each of the form
	                 (n, kappa, E, np, pz, r, pr, qr)
	"""
	def __init__(self,file_name,orbitals_data):
		self.f = file_name
		self.orbitals_data = orbitals_data
	#
	def REFINE(self,x,index):
		"""
		index = 0   for energy eigenvalue
		index = 1 for P(r), Q(r) or r values
		"""
		x2 = x.split('E')
		power = x2[1]
		left = x2[0].split('.')
		# check if it is positive or negative
		if '-' in left[0]:
			if index == 0: # for eigenvalues
				my_left = ' 0.'
			else:
				my_left = '-0.'
		else:
			my_left = ' 0.'
		my_num = str( abs(int(left[0])) )+ left[1]
		power = int(x2[1]) + 1
		if -10 < power < 10:
			if '-' in str(power):
				f_power = '-0'+ str(abs(power))
			else:
				f_power = '+0'+str(abs(power))
		else:
			f_power = str(power)
		out = my_left+my_num +'E'+f_power
		# modify values of zero
		out_grasp = my_left+my_num +'D'+f_power
		# zero values
		if abs(float(out)) == 0 :
			return "  "+my_left+my_num +'D+00'
		# other values than zero
		else:
			return "  "+out_grasp
	#
	def write_grasp(self):
		f = open(self.f,"w")
		f.write("G92RWF"+"\n")   # header written ONCE for the whole file
		#
		for orb in self.orbitals_data:
			n,kappa,E,nnp,pz,r,pr,qr = orb
			#
			line1 = '{:>8d}'.format(n)+'{:>8d}'.format(kappa)+\
			self.REFINE( '{:>25.15E}'.format(E), 0)+\
			'{:>8d}'.format(nnp)
			f.write(line1+"\n")
			# pz
			lne1 = self.REFINE('{:>25.15E}'.format(pz),1)
			f.write(lne1+"\n")
			# P(r)
			for i in range(nnp):
				 lnx1 = self.REFINE( '{:>25.15E}'.format(pr[i]),1)
				 f.write(lnx1+"\n")
			# Q(r)
			for i in range(nnp):
				 lnx2 = self.REFINE( '{:>25.15E}'.format(qr[i]),1)
				 f.write(lnx2+"\n")
			# r
			for i in range(nnp):
				 lnx3 = self.REFINE( '{:>25.15E}'.format(r[i]),1)
				 f.write(lnx3+"\n")
		#
		f.close()
#----------------#
class WRITE_SUMMARY:
	"""
	write a summary table of ALL calculated orbitals in this run:
	orbital label, n, kappa, energy (a.u. and eV), <r>, <r^2>, and NP
	(the number of radial grid points in that orbital's final
	solution). Also records the atomic data (Z, Ne) and the potential
	type/parameters used for the whole run, since these are shared by
	every orbital in the sequence.
	input:
	 file_name     = output file name
	 zz            = atomic number
	 nn            = number of electrons
	 potential_label = descriptive string of the potential type/
	                   parameters used (e.g. "Coulomb potential
	                   (unscreened)", "Yukawa potential (alpha = 0.5)")
	 summary_data  = list of tuples, one per orbital, each of the form
	                 (label, n, kappa, E_au, E_eV, r_exp, r2_exp, npx)
	"""
	def __init__(self,file_name,zz,nn,potential_label,potential_constants,\
			summary_data,not_converged=None):
		self.f = file_name
		self.zz = zz
		self.nn = nn
		self.potential_label = potential_label
		self.potential_constants = potential_constants
		self.summary_data = summary_data
		self.not_converged = not_converged if not_converged is not None else []
	#
	def not_converged_lines(self):
		lines = []
		if len(self.not_converged) == 0:
			lines.append("  (none -- all requested orbitals converged)")
		else:
			for item in self.not_converged:
				label,n,kappa,reason = item
				lines.append("  "+'{:<10}'.format(label)+\
					"(n = "+str(n)+", kappa = "+str(kappa)+")  -- "+reason)
		return lines
	#
	def potential_constants_lines(self):
		lines = []
		if len(self.potential_constants) == 0:
			lines.append("  (no additional free parameters)")
		else:
			for name,value in self.potential_constants:
				lines.append("  "+'{:<10}'.format(name)+" = "+str(value))
		return lines
	#
	def header_line(self):
		return '{:<10}{:>6}{:>8}{:>20}{:>16}{:>16}{:>16}{:>6}{:>8}'.format(\
			"Orbital","n","kappa","E (a.u.)","E (eV)","<r>","<r^2>","Pot","NP")
	#
	def row_line(self,row):
		label,n,kappa,E_au,E_eV,r_exp,r2_exp,pot_idx,npx = row
		return '{:<10}{:>6d}{:>8d}{:>20.10E}{:>16.6E}{:>16.8E}{:>16.8E}{:>6d}{:>8d}'.\
			format(label,n,kappa,E_au,E_eV,r_exp,r2_exp,pot_idx,npx)
	#
	def write(self):
		f = open(self.f,"w")
		f.write("\n")
		lk1 = "atomic number (Z)         = "+'{:<8d}'.format(self.zz)
		lk2 = "number of electrons (Ne)  = "+'{:<8d}'.format(self.nn)
		lk3 = "potential used            = "+self.potential_label
		f.write(lk1+"\n")
		f.write(lk2+"\n")
		f.write(lk3+"\n")
		f.write("potential constants:\n")
		for line in self.potential_constants_lines():
			f.write(line+"\n")
		f.write("\n")
		ln0 = self.header_line()
		f.write(ln0+"\n")
		f.write("-"*len(ln0)+"\n")
		for row in self.summary_data:
			f.write(self.row_line(row)+"\n")
		f.write("\n")
		f.write("Orbitals NOT converged (excluded from the table above):\n")
		for line in self.not_converged_lines():
			f.write(line+"\n")
		f.close()
	#
	def print_summary(self):
		lk1 = "atomic number (Z)         = "+'{:<8d}'.format(self.zz)
		lk2 = "number of electrons (Ne)  = "+'{:<8d}'.format(self.nn)
		lk3 = "potential used            = "+self.potential_label
		print(lk1)
		print(lk2)
		print(lk3)
		print("potential constants:")
		for line in self.potential_constants_lines():
			print(line)
		print("")
		ln0 = self.header_line()
		print(ln0)
		print("-"*len(ln0))
		for row in self.summary_data:
			print(self.row_line(row))
		print("")
		print("Orbitals NOT converged (excluded from the table above):")
		for line in self.not_converged_lines():
			print(line)
#----------------#
class WRITE_POTENTIAL_TF:
	"""
	NOt used
	this version  for Thomas-Fermi only
	"""
	def __init__(self,file_name,r,rv,phi):
		self.f = file_name
		self.r = r
		self.rv = rv
		self.phi = phi
	#
	def write_potential(self):
		# determine the half of the data
		np = len(self.r)
		if np%2 == 0: # even
			half_index = int(np/2)
			count = half_index
		else:
			half_index = int( (np+1)/2)
			count = half_index
		# open file
		f = open(self.f,"w")
		f.write("\n")   # new empty line
		f.write("Thomas-Fermi potential"+"\n\n")   
		lne1 = "    ========= r =================== r * V(r) ===============  phi(r) ======"
		f.write(lne1+lne1+"\n")
		for i in range(half_index-1):
			lnx1 = '{:>25.15E}'.format(self.r[i])+\
				'{:>25.15E}'.format(self.rv[i])+\
				'{:>25.15E}'.format(self.phi[i])
			lnx2 = '{:>25.15E}'.format(self.r[count])+\
				'{:>25.15E}'.format(self.rv[count])+\
				'{:>25.15E}'.format(self.phi[count])
			f.write(lnx1+lnx2+"\n")
			count+=1
		# write the last line
		if np%2 == 0: # even
			lnx3 = '{:>25.15E}'.format(self.r[half_index-1])+\
				'{:>25.15E}'.format(self.rv[half_index-1])+\
				'{:>25.15E}'.format(self.phi[half_index-1])
			lnx4 = '{:>25.15E}'.format(self.r[count-1])+\
				'{:>25.15E}'.format(self.rv[count-1])+\
				'{:>25.15E}'.format(self.phi[count-1])
			f.write(lnx3+lnx4+"\n")
		else: # odd
			lnx5 = '{:>25.15E}'.format(self.r[half_index-1])+\
				'{:>25.15E}'.format(self.rv[half_index-1])+\
				'{:>25.15E}'.format(self.phi[half_index-1])
			f.write(lnx5+"\n")
		#
		f.close()

#----------------#
class WRITE_POTENTIAL:
	"""
	general form for all potential types
	"""
	def __init__(self,file_name,r,rv):
		self.f = file_name
		self.r = r
		self.rv = rv
	#
	def write_potential(self):
		# determine the half of the data
		np = len(self.r)
		# open file
		f = open(self.f,"w")
		f.write("\n")   # new empty line
		#f.write("Thomas-Fermi potential"+"\n\n")  
		f.write("Potential function used"+"\n\n")  
		lne1 = "    ========= r ===================    V(r) ==============="
		f.write(lne1+"\n")
		for i in range(np):
			if float(self.r[i]) > 0.0:
				poten = float(self.rv[i]) / float(self.r[i])
				lnx1 = '{:>25.15E}'.format(self.r[i])+\
					'{:>25.15E}'.format(poten)
				f.write(lnx1+"\n")
		#
		f.close()
#----------------#
