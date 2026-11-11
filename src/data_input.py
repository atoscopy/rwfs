#!/usr/bin/python3

from constants import eta, max_Z


class ATOM:
	#====#
	def atomic_number(self,d):
		atomic_z = input(d) 
		# check correction 
		while atomic_z.strip() == "" or int(atomic_z)<=0 or\
		 int(atomic_z) > max_Z :
			dx1 = " ERROR: please enter correct atomic number"
			print(dx1)
			atomic_z = input(d) 
		if int(atomic_z) <= 0:
			er1 = " ERROR: (atomic_number) the atomic number"
			er2 = " must be >= 1 (like Hydrogen)"
			print(er1+er2)
		if int(atomic_z) > max_Z:
			erx1 = " ERROR: (atomic_number) the atomic number"
			erx2 = "must be <= 109 (like Meitnerium)"
			print(erx1+erx2)
		# 
		return int(atomic_z)
	#====#
	def number_electrons(self,zz,p):
		cont = False
		# zz is the atomic number
		t1 = " If the atom is neutral ( press Enter <--: )"
		t2 = " i.e. the atomic number = the number of electrons"
		print(t1)
		print(t2)
		#p = " Enter the number of electrons in the atom/ion:"
		print(p)
		num = input()
		if str(num).lstrip() == "":
			num = int(zz)
			cont = True
		else:
			if int(num) > int(zz):
				er = " ERROR: (number_electrons) the ion should be +ve ion"
				print(er)
			if int(num) < 0:
				er2 = " ERROR: (number_electrons) the number of electrons"
				er3 = "   >= 1 (Hydrogen case)"
				print(er2)
				print(er3)
		while  int(num)> int(zz) or int(num) < 0 or cont == False :
			print(" Error: please enter correct number of electrons")
			print(p)
			num = input()
			if str(num).lstrip() == "" or int(num)== int(zz):
				num = int(zz)
				cont = True
		return int(num)
#-----------------------#
class DEFAULT:
	"""
	In case of changing the parameters used
	"""
	#
	def yes_or_no(self,x_list,y_list):
		for i in x_list:
			print(i)
		yn = input()
		while yn not in ["y","Y","n","N"]:
			for k in y_list:
				print(k)
			yn = input()
		return yn
	#
	def change_defaults(self):
		t1 = " Do you want to change the default parameters (y/n):"
		x_list = [t1]
		t2 = " ...the answer should be <y> or <n>"
		y_list = [t2,t1]
		include_x = self.yes_or_no(x_list,y_list)
		return include_x
	#
	def change_relative_error(self):
		t1 = " Do you want to change the relative error (y/n):"
		x_list = [t1]
		t2 = " ...the answer should be <y> or <n>"
		y_list = [t2,t1]
		include_x = self.yes_or_no(x_list,y_list)
		return include_x
	#
	def change_grid(self):
		x1 = " Do you want to change the grid parameters?"
		x_list = [x1]
		t2 = " ...the answer should be <y> or <n>"
		y_list = [t2,x1]
		include_x = self.yes_or_no(x_list,y_list)
		return include_x
	#
	def print_TF_potential(self):
		x1 = " Do you want to print the potential function ( V(r) )?"
		#x2 = " and the Thomas-Fermi function (Phi)?"
		#x_list = [x1,x2]
		x_list = [x1]
		t2 = " ...the answer should be <y> or <n>"
		#y_list = [t2,x1,x2]
		y_list = [t2,x1]
		include_x = self.yes_or_no(x_list,y_list)
		return include_x
	#
	def grasp_output(self):
		x1 = " Do you want to print the output orbital in GRASP format?"
		x_list = [x1]
		t2 = " ...the answer should be <y> or <n>"
		y_list = [t2,x1]
		include_x = self.yes_or_no(x_list,y_list)
		return include_x
#-----------------------#
class DATA:
	"""
	This class to control all the non default output
	"""
	# ===============#
	def check_number(self,user_number):
		try:
			data = float(user_number)
			return True
		except ValueError:
			return False
	#
	def check_number_integer(self,user_number):
		try:
			data = int(user_number)
			return True
		except ValueError:
			return False
	# ===============#
	# eta error
	def eta_error(self):
		ch_eta = DEFAULT().change_relative_error()
		#
		if ch_eta in ["y","Y"]:
			p1 = " The relative error used value is  "
			p2 = " enter the new relative error (eta) you want"
			print(p1,eta)
			print(p2)
			new_eta = input()
			while self.check_number(new_eta) == False:
				p3 = " Error: not correct input"
				print(p3)
				print(p2)
				new_eta = input()
			#
			return float(new_eta)
		else:
			return eta
	#
	def grid_type(self):
		print(" There are two kinds of grids")
		print("\n")
		print(" type (1) the default")
		print(" defined by r_i = r_{i-1} x d_r")
		print(" i = 1,2,3,....,np")
		print(" with r_0 = 0.0")
		print("\n")
		print(" type (2) (grasp)")
		print(" defined by r_i = r_0 (e^{t_i} - 1)")
		print(" where r_0 = T/Z, t_i = h(i-1), i = 1,2,3,...,np")
		#
		print("\n")
		print(" print 1 for type (1) or 2 for type (2)")
		g_type = input() 
		while g_type.strip() == "" or g_type not in ["1","2"]:
			print(" Error: please enter correct number either 1 or 2")
			g_type = input() 
		return int(g_type)
	#
	# grid parameters
	def non_default_grid(self):
		ch_grid = DEFAULT().change_grid()
		if ch_grid in ["y","Y"]:
			tp = self.grid_type()
			if tp == 1:
				p1 = " The grid points used are"
				p2 = " np = 3000"
				p3 = " r1 = 1.0E-8"
				p4 = " dr = 1.012"
				p5 = " enter the new grid points you want"
				p6 = " enter np: "
				print(p1)
				print(p2)
				print(p3)
				print(p4)
				print(p5)
				print(p6)
				new_np = input()
				while self.check_number_integer(new_np) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p5)
					print(p6)
					new_np = input()
				#
				p7 = " Enter r1: "
				print(p7)
				new_r1 = input()
				while self.check_number(new_r1) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p7)
					new_r1 = input()
				#
				p8 = " Enter dr: "
				print(p8)
				new_dr = input()
				while self.check_number(new_dr) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p8)
					new_dr = input()
				## 1 for type 1
				return [1,new_np, new_r1, new_dr]
			if tp == 2:
				p1 = " The default grasp grid points used are"
				p2 = " np = 590    (N)"
				p3 = " T  = 1.0E-5 (H)"
				p4 = " H  = 0.05   (RNT)"
				p5 = " enter the new grid points you want"
				p6 = " enter np: "
				print(p1)
				print(p2)
				print(p3)
				print(p4)
				print(p5)
				print(p6)
				new_np = input()
				while self.check_number_integer(new_np) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p5)
					print(p6)
					new_np = input()
				#
				p7 = " Enter T: "
				print(p7)
				new_r1 = input()
				while self.check_number(new_r1) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p7)
					new_r1 = input()
				#
				p8 = " Enter H: "
				print(p8)
				new_dr = input()
				while self.check_number(new_dr) == False:
					x1 = " Error: not correct input"
					print(x1)
					print(p8)
					new_dr = input()
				## 2 for type 2
				return [2,new_np, new_r1, new_dr]
		else:
			return [0]
	#
	def TF_potential(self):
		prt_tf = DEFAULT().print_TF_potential()
		if prt_tf in ["y","Y"]:
			p1 = " the potential data r, V(r) are saved in file name"
			p2 = " 'potential' if you want another file name, enter"
			p3 = " the file name: (press enter to save the default)"
			print(p1)
			print(p2)
			print(p3)
			pot_file =  input()
			#
			pot_file = pot_file.strip()
			#
			if pot_file== "":
				pot_name = "potential"
			else:
				pot_name = pot_file
		else:
			pot_name =""
		return pot_name
	#
	def prt_grasp_format(self):
		grasp_output = DEFAULT().grasp_output()
		if grasp_output in ["y","Y"]:
			p1 = " to save the output as in form used in GRASP package"
			p2 = " the default file name is 'grasp.rwf', enter the "
			p3 = " file name: (press enter to save the default)"
			print(p1)
			print(p2)
			print(p3)
			grasp_file =  input()
			#
			grasp_file = grasp_file.strip()
			if grasp_file== "":
				grasp_name = "grasp.rwf"
			else:
				grasp_name = grasp_file
		else:
			grasp_name = ""
		return grasp_name
