#!/usr/bin/python3

from constants import  max_Z
from physics import PHYSICS

# data input


class ATOMIC:
	"""
	define the atomic number used
	"""
	#
	def atomic_range_Z(self):
		z_max = max_Z       
		atomic = [""]  
		for i in range(z_max+1):
			z = str(i+1)
			atomic.append(z)
		return atomic
	
	#
	def atomic_number(self):
		d1 = " Enter the atomic number (Z) of the atom:"
		print(d1)
		atomic_z = input() 
		# atomic range
		atomic = self.atomic_range_Z()
		# check correction 
		while atomic_z.strip() not in atomic or int(atomic_z) == 0:
			print(" Error: please enter correct atomic number")
			print(d1)
			atomic_z = input() 
		# 
		return int(atomic_z)
	#
	def number_electrons(self,zz):
		# zz is the atomic number
		t0 = " Enter the number of electrons in the atom/ion:"
		t1 = " Positive or negative ions are all allowed."
		t2 = " If the atom is neutral ( press Enter )"
		t3 = " i.e., the atomic number = the number of electrons"
		print(t0)
		print(t1)
		print(t2)
		print(t3)
		num = input()
		# all atomic possible range
		atomic = self.atomic_range_Z()
		#
		while num.strip() not in atomic:
			print(" Error: please enter correct number of electrons")
			print(t0)
			num = input()
		if str(num).strip() == "":
			num = int(zz)
		if int(num) > int(zz):
			er = " this is negative ion"
			print(er)
		if int(num) < int(zz):
			er1 = " this is positive ion"
			print(er1)
		return int(num)
#----------------#
# define the orbital used
class ORBITAL:
	"""
	this class determines the orbital used in our calculations
	"""
	# data used in all
	shell = ['s','p','d','f','g','h','i','k','l','m','n']
	
	# 1-
	def orbital_type(self):
		px1 = " Enter the orbital you want either as spectroscopic notation"
		px2 = " as (1s, 2p-, 2p, ...) or as (n,kappa) values without brackets."
		px3 = " Press ( 1 ) for spectroscopic notation or press ( 2 ) for"
		px4 = " n,kappa  notation: "
		print(px1)
		print(px2)
		print(px3)
		print(px4)
		orbital = input() 
		while orbital.strip() == "" or orbital.strip() not in ["1","2"]:
			print(" Error: please enter correct number either 1 or 2")
			orbital = input() 
		return int(orbital)
	#
	#===============#
	# aided functions for spectroscopic notations
	# define the orbital
	def read_datax(self,orbitalx):
		if orbitalx.find(",")== -1:
			return "spectroscopic"
		else:
			return "nkappa"
	#
	def loop_spectros(self,data,shell):
		#
		# case1: if shell are in shell list
		#
		out = [] # all output
		#
		shell_out = []
		for s in shell:
			if s in data:
				shell_out.append(s)
		if len(shell_out) == 1:
			out.append(True)
		else:
			out.append(False)
		# add if space is found
		m1 = data.strip()
		m2 = m1.find(" ")
		if m2 == -1:
			out.append(True)
		else:
			out.append(False)
		# check if false is found
		k = True
		for i in out:
			if i == False:
				k = False
		# case 2
		# if n ,shell not compatible
		error = 1
		if k == True:
			sh = shell_out[0]
			ns = data.split(sh)
			n = int(ns[0])
			# check if n, shell are correct 
			if n <= len(shell):
				for tt in range(n):
					if shell[tt] == sh.strip():
						error = 0
			else:
				for rr in range(len(shell)):
					if shell[rr] == sh.strip():
						error = 0
			# 
			if error == 0:
				return True
			else:
				return False
		else:
			return False
	#
	def nk_extract2(self,orbtx):
		if orbtx.find(",") == -1:
			return False
		else:
			return True
	#
	def loop_nkappa2(self,data):
		if self.nk_extract2(data) == True:
			return True
		else:
			return False
	#
	def angular_l(self,orbital_shell,shell):
		# print l values for each shell
		l_out = -1
		l = 0
		for t in range(len(shell)):
			if shell[t] == orbital_shell:
				l_out =  l
			l+=1
		return l_out
	#
	def kappa2(self,l,outsign):
		# return kappa values for each l
		if outsign == "-":
			ka = l
		else:
			ka = -(l+1)
		return ka
	# 2 
	def orbital_spectroscopic(self,shell):
		px1 = " Enter the orbital you want in spectroscopic notation"
		px2 = " as 1s, 2p-, 2p, ... , 5g-, 5g, ... etc."
		print(px1)
		print(px2)
		orbital = input() 
		x = orbital
		# x is the input 
		datam = self.read_datax(x)
		# check if 
		check = 1
		while self.loop_spectros(x,shell) == False or datam !=\
			"spectroscopic" or check  == 0:
			print(" ERROR: not physical spectroscopic input: ",x)
			print(px1)
			print(px2)
			orbital = input() 
			#
			x = orbital
			if self.loop_nkappa2(x) == True:
				check  = 0
		#
		m1 = x.strip()
		#
		shell_out = []
		for s in shell:
			if s in m1:
				shell_out.append(s)
		if len(shell_out) == 1:
			pass
		else:
			print( " ERROR: check the orbital ",m1)
		#
		sh = shell_out[0]
		m2 = m1.split(sh)
		n = int( m2[0]) 
		sh_sign  = m2[1]
		#
		if sh_sign=='-':
			outsign = '-'
		else:
			outsign = '+'
		#
		L = self.angular_l(sh,shell)
		#
		#
		kx = self.kappa2(L,outsign)
		# # return n , kappa values
		return int(n) , int(kx)   
	#===============#
	# aid functions for n, kappa notation
	#
	def kappax(self,n):
		# return kappa values for each n
		kappa = []
		for i in range(-n,n,1):
			if i == 0:
				pass
			else:
				kappa.append(i)
		return kappa
	#
	def nkapppa_check(self,n,kappa):
		# return True if correct n,kappa and False if not
		k = self.kappax(n)
		if kappa in k:
			return True 
		else:
			return False
	#
	def nk_extract(self,orbtx):
		m = orbtx.split(",")
		n  = int( m[0] )
		kappa = int( m[1] )
		return n,kappa
	#
	def loop_nkappa(self,data):
		nn , kk = self.nk_extract(data)
		# check if data is correct 
		if self.nkapppa_check(nn,kk) == True:
			return True
		else:
			return False
	#		
	# 3
	def orbital_nk(self,shell):
		px1 = " Enter the orbital you want in (n,kappa) notation "
		px2 = " without brackets as 1,-1  or 3,-3  ...... etc."
		print(px1)
		print(px2)
		orbital = input() 
		x = orbital
		# x is the input 
		datam = self.read_datax(x)
		# check if 
		check = 1
		#
		while self.loop_nkappa2(x)==False or datam != "nkappa" or\
		 check == 0 or self.loop_nkappa(x) == False:
			print(" ERROR: not physical n,kappa input: ",x)
			print(px1)
			print(px2)
			orbital = input() 
			#
			x = orbital
			if self.loop_spectros(x,shell) == True:
				check  = 0
		#
		nn , kk = self.nk_extract(x)
		return int(nn) , int(kk) # return n , kappa values
	#===============#
	def read(self):
		shell = ['s','p','d','f','g','h','i','k','l','m','n']
		nkappa_or_spectroscopic = self.orbital_type()
		if nkappa_or_spectroscopic == 1:
			n,kappa = self.orbital_spectroscopic(shell)
			return n,kappa 
		else:
			n,kappa = self.orbital_nk(shell)
			return n,kappa 
	#===============#
	# parse a single spectroscopic-notation token (e.g. "1s","2p-","3d-")
	# returns (n,kappa) or None if the token is not a valid orbital
	def parse_spectroscopic_token(self,token,shell):
		tok = token.strip()
		if tok == "":
			return None
		if self.loop_spectros(tok,shell) == False:
			return None
		#
		shell_out = []
		for s in shell:
			if s in tok:
				shell_out.append(s)
		if len(shell_out) != 1:
			return None
		#
		sh = shell_out[0]
		parts = tok.split(sh)
		if len(parts) != 2:
			return None
		try:
			n = int(parts[0])
		except ValueError:
			return None
		sh_sign = parts[1]
		outsign = '-' if sh_sign == '-' else '+'
		#
		L = self.angular_l(sh,shell)
		if L == -1:
			return None
		kx = self.kappa2(L,outsign)
		return int(n),int(kx)
	#
	# parse a single (n,kappa)-notation token, with or without
	# surrounding parentheses, e.g. "(1,-1)","1,-1","2,-3"
	# returns (n,kappa) or None if the token is not a valid orbital
	def parse_nkappa_token(self,token):
		tok = token.strip().strip("()")
		if "," not in tok:
			return None
		parts = tok.split(",")
		if len(parts) != 2:
			return None
		try:
			n = int(parts[0].strip())
			kappa = int(parts[1].strip())
		except ValueError:
			return None
		if self.nkapppa_check(n,kappa) == False:
			return None
		return n,kappa
	#
	# parse a whole line of space-separated orbital tokens, in either
	# notation. returns (orbitals, bad_tokens)
	def parse_orbitals_line(self,line,notation,shell):
		tokens = line.split()
		orbitals = []
		bad_tokens = []
		for tok in tokens:
			if notation == 1:
				res = self.parse_spectroscopic_token(tok,shell)
			else:
				res = self.parse_nkappa_token(tok)
			if res is None:
				bad_tokens.append(tok)
			else:
				orbitals.append(res)
		return orbitals,bad_tokens
	#
	#shell_letters = ['s','p','d','f','g','h','i','k','l','m','n']
	def orbital_label(self,nx,kappax):
		lx = PHYSICS(nx,kappax).angular_momentum()
		sign = '-' if kappax > 0 else ''
		return str(nx)+self.shell[lx]+sign
	#
	# read a sequence (list) of orbitals in one line, all sharing the
	# same notation, which is chosen once at the start
	def read_multiple(self):
		shell = self.shell
		#
		# 1- choose the notation once for the whole sequence
		notation = self.orbital_type()
		#
		if notation == 1:
			px1 = " Enter ALL the orbitals you want, separated by spaces,"
			px2 = " in spectroscopic notation, for example:"
			px3 = " 1s 2s 2p- 2p 3d- 4f"
		else:
			px1 = " Enter ALL the orbitals you want, separated by spaces,"
			px2 = " in (n,kappa) notation (parentheses optional, no"
			px3 = " space inside a pair), for example: (1,-1) (2,-1) 2,1 3,-2"
		#
		print(px1)
		print(px2)
		print(px3)
		line = input()
		orbitals,bad_tokens = self.parse_orbitals_line(line,notation,shell)
		#
		while len(bad_tokens) > 0 or len(orbitals) == 0:
			if len(bad_tokens) > 0:
				print(" ERROR: could not parse the following orbital(s): ",\
					  " ".join(bad_tokens))
			else:
				print(" ERROR: please enter at least one orbital")
			print(px1)
			print(px2)
			print(px3)
			line = input()
			orbitals,bad_tokens = self.parse_orbitals_line(line,notation,shell)
		#
		
		# summary
		#print(" ----------------------------------------")
		print(" The following ",len(orbitals)," orbital(s) will be calculated:")
		for idx in range(len(orbitals)):
			nn,kk = orbitals[idx]
			label = self.orbital_label(nn,kk)
			print("  ",idx+1,") n =",nn,", kappa =",kk,", spec = ",label)
		#print(" ----------------------------------------")
		return orbitals
# END
