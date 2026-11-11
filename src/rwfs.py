#!/usr/bin/python3


from grid import GRID,GRID_90,GRID_C
from potential import THOMAS_FERMI,COULOMB,YUKAWA,MULTI_YUKAWA,\
				SCREENED,SCREENED_HYDROGENIC,READ_pote_file_rv,\
				UNIFORM_SPHERE,FERMI_NUCLEUS
from spline import SPLINE,NORM_INTEGERAL,INTEG_VALUE,PQ_NORMALIZED,\
			Cubic_Spline_Interpolation_Function
from physics import PHYSICS,BOUND_EMIN_EMAX,INWARD_OUTWARD_DATA,\
		RESCALING_INWARD
from outward import KAPPA_SIGN,PQ_OUTWARD,PZ
from inward import R_INF,PQ_INWARD,ENERGY_CORRECTION
from energy_loop import ENERGY_LOOP
from constants import E_init,eta,Hartree_to_eV
from data_input import ATOM,DATA,DEFAULT
from data_input2 import ORBITAL,ATOMIC
from data_input3 import POTENTIAL
from write import WRITE_DATA,WRITE_DATA_MULTI,WRITE_DATA_GRASP,\
			WRITE_DATA_GRASP_MULTI,WRITE_SUMMARY,WRITE_POTENTIAL
#-----#

#----------------#
# ATOMIC DATA
Z = ATOMIC().atomic_number()
N_ele = ATOMIC().number_electrons(Z)
print(" The number of electrons in atom/ion is ",N_ele)
#
dbo = ORBITAL()
# read a SEQUENCE of orbitals to be calculated in this run
orbitals = dbo.read_multiple()
#
# keep the initial estimate of the energy as given in constants.py so
# that every orbital in the sequence starts from the same default
# estimate, rather than inheriting the converged energy of a previous
# orbital in the list
E_init_default = E_init
#
# spectroscopic-style label for an (n,kappa) pair, used to build
# distinct output file names for each orbital in the sequence
shell_letters = ['s','p','d','f','g','h','i','k','l','m','n']
def orbital_label(nx,kappax):
	lx = PHYSICS(nx,kappax).angular_momentum()
	sign = '-' if kappax > 0 else ''
	return str(nx)+shell_letters[lx]+sign
#----------------#
#check default values
use_default = DEFAULT().change_defaults()
if use_default in ["y","Y"]:
	# 1
	eta = DATA().eta_error()
	# 2
	# grid type
	new_grid = DATA().non_default_grid()
	if len(new_grid) == 1: # default values
		rg = GRID_90(3000,1.0E-8,1.012)
		r = rg.grid_90()
	if len(new_grid) == 4:  # non default values
		if new_grid[0]==1 :  # type 1
			rg = GRID_90(int(new_grid[1]),float(new_grid[2]),\
			float(new_grid[3]))
			r = rg.grid_90()
		if new_grid[0]==2 :  # type 2
			rg = GRID(int(new_grid[1]),float(new_grid[2]),\
			float(new_grid[3]),float(Z) )
			r = rg.grid()
	# 3
	prt_potential = DATA().TF_potential()
	# 4
	grasp_prt = DATA().prt_grasp_format()
else:
	# grid
	rg = GRID_90(3000,1.0E-8,1.012)
	r = rg.grid_90()
	# potential and grasp output will be removed
	prt_potential = ""
	grasp_prt = ""
# grid points
np = len(r)
########################
# define the potential used for the atom/ion
poten1 = POTENTIAL()
pot_form = poten1.define_potential()
# now define the new parameters 
if pot_form == 1:
	TF = THOMAS_FERMI(r,Z,N_ele)
	phi = TF.phi()
	rv = TF.potential()
if  pot_form == 2:
	rv = COULOMB(r,Z).potential()
if  pot_form == 3:
	alpha_r = POTENTIAL().alpha()
	rv = YUKAWA(r,Z,alpha_r).potential()
if  pot_form == 4:
	nmax,terms = POTENTIAL().multi_yukawa_parameters()
	rv = MULTI_YUKAWA(r,Z,nmax,terms).potential()
if  pot_form == 5:
	alpha_s = POTENTIAL().alpha()
	rv = SCREENED(r,Z,alpha_s).potential()
if  pot_form == 6:
	alpha_6 = POTENTIAL().alpha()
	rv = SCREENED_HYDROGENIC(r,Z,alpha_6).potential()
if  pot_form == 7:
	file_name_input = POTENTIAL().potent_file()
	r_new,v_new = READ_pote_file_rv(file_name_input).read()
	if len(r_new) < 2:
		mx = " ERROR: potential file "
		mx2 = "has fewer than 2 data points"
		raise ValueError(mx+file_name_input+mx2)
	# now interpolate on the default grid 
	# avoid r = 0.0
	if r_new[0] == 0.0:
		r_new = r_new[1:]
		v_new = v_new[1:]
	# create new r * V(r)
	rv_new = []
	for t in range(len(r_new)):
		p1 = r_new[t]*v_new[t]
		rv_new .append(p1)
	# define the spline coeffients for the potential rV(r)
	a1,b1,c1,d1 = SPLINE(r_new,rv_new).spline()
	interpolate = Cubic_Spline_Interpolation_Function(r_new,\
					a1,b1,c1,d1)
	# define the last term in the default r 
	rmax = max(r_new)
	# effective asymptotic charge implied by the file's own tail
	Z_eff_tail = -rv_new[-1]     # since rV(r) -> -Z_eff at large r

	rv = []
	for rr in r:
		if rr <= rmax:
			rv.append(interpolate(rr))
		else:
			# analytic Coulomb tail beyond the tabulated data
			rv.append(-Z_eff_tail)
	#
	# now define np again
	np = len(r)
if  pot_form == 8:
	# uniformly charged sphere, finite nuclear size
	R_rms_nuc = POTENTIAL().nuclear_rms_radius()
	nuc_model = UNIFORM_SPHERE(r,Z,R_rms_nuc)
	rv = nuc_model.potential()
if  pot_form == 9:
	# Fermi charge distribution, finite nuclear size
	R_rms_nuc = POTENTIAL().nuclear_rms_radius()
	t_skin = POTENTIAL().fermi_skin_thickness()
	nuc_model = FERMI_NUCLEUS(r,Z,R_rms_nuc,t_skin)
	rv = nuc_model.potential()

# build a human-readable label describing the potential type/parameters
# used for this run, shared by every orbital and recorded in the summary
potential_names = {
	1: "Thomas-Fermi potential",
	2: "Coulomb potential (unscreened)",
	3: "Yukawa potential",
	4: "Multi Yukawa potential",
	5: "Screened potential",
	6: "Screened-hydrogenic potential",
	7: "Potential read from input file",
	8: "Uniformly charged sphere (finite nucleus)",
	9: "Fermi charge distribution (finite nucleus)",
}
potential_label = potential_names.get(pot_form,"Unknown potential")
if pot_form == 3:
	potential_label += " (alpha = "+str(alpha_r)+")"
if pot_form == 4:
	potential_label += " (n_max = "+str(nmax)+")"
if pot_form == 5:
	potential_label += " (alpha = "+str(alpha_s)+")"
if pot_form == 6:
	potential_label += " (alpha = "+str(alpha_6)+")"
if pot_form == 7:
	potential_label += " (file = "+file_name_input+")"
if pot_form == 8:
	potential_label += " (R_rms = "+str(R_rms_nuc)+" fm)"
if pot_form == 9:
	potential_label += " (R_rms = "+str(R_rms_nuc)+" fm, t = "+\
					   str(t_skin)+" fm)"

# full list of ALL free constants/parameters used to define the
# potential, for the detailed summary (not just the short label above)
potential_constants = []
if pot_form == 1:
	potential_constants.append( ("Z",Z) )
	potential_constants.append( ("N_e",N_ele) )
if pot_form == 2:
	potential_constants.append( ("Z",Z) )
if pot_form == 3:
	potential_constants.append( ("Z",Z) )
	potential_constants.append( ("alpha",alpha_r) )
if pot_form == 4:
	potential_constants.append( ("Z",Z) )
	potential_constants.append( ("n_max",nmax) )
	for t_idx in range(len(terms)):
		c_t,a_t = terms[t_idx]
		potential_constants.append( ("c_"+str(t_idx+1),c_t) )
		potential_constants.append( ("alpha_"+str(t_idx+1),a_t) )
if pot_form == 5:
	potential_constants.append( ("Z",Z) )
	potential_constants.append( ("alpha",alpha_s) )
if pot_form == 6:
	potential_constants.append( ("Z",Z) )
	potential_constants.append( ("alpha",alpha_6) )
if pot_form == 7:
	potential_constants.append( ("file",file_name_input) )
	potential_constants.append( ("Z_eff_tail",Z_eff_tail) )
if pot_form == 8 or pot_form == 9:
	potential_constants.append( ("Z",Z) )
	for nm,val in nuc_model.constants():
		potential_constants.append( (nm,val) )

# define the spline coeffients for the potential r.V(r)
A,B,C,D = SPLINE(r,rv).spline()
#
# print the potential once, since it only depends on Z/N_ele/pot_form,
# not on the individual orbitals
if prt_potential == "":
	pass
else:
	potential_file_name = prt_potential
	pot = WRITE_POTENTIAL(potential_file_name,r,rv)
	pot.write_potential()

########################
# loop over the requested sequence of orbitals
########################
n_total = len(orbitals)
#
# accumulate every orbital's results here; all orbitals are written
# together into ONE combined file after the loop finishes, rather
# than one file per orbital
radial_data_all = []
grasp_data_all = []
summary_data = []
not_converged_data = []
#
for orb_idx in range(n_total):
	n_qn,kappa_qn = orbitals[orb_idx]
	n     = n_qn
	kappa = kappa_qn
	#
	label = orbital_label(n_qn,kappa_qn)
	print("\n")
	print(" ================================================")
	print(" Orbital ",orb_idx+1,"/",n_total,": ",label,\
		  " (n =",n_qn,", kappa =",kappa_qn,")")
	print(" ================================================")
	#
	ph = PHYSICS(n,kappa)
	N_0 = ph.actual_nodes()
	l   = ph.angular_momentum()
	#
	print(" Starting calculations ...")
	print("\n")
	print(" the orbital defined by n =",n,",kappa =",kappa)
	print(" has ",N_0," node/s.")
	#
	# energy limits for calculation
	# always start from the default initial estimate given in
	# constants.py, not from a previous orbital's converged energy
	energy = BOUND_EMIN_EMAX(n,kappa,rv,r,E_init_default)
	EMIN,EMAX,E = energy.e_min_and_e_max()
	print(" EMIN = ",EMIN," EMAX = ",EMAX)
	print(" initial E used = ",E)
	print("\n")
	# start the outward integration
	E_init_orb = E
	pqout = PQ_OUTWARD(N_0,A,B,C,D,r,rv,np,E_init_orb,EMIN,EMAX,l,kappa,eta,0)
	IOTP,E,Pr_outward,Qr_outward,P_out,Q_out,stab_zero,cal_nodes =\
	   pqout.pq_outward()
	### check if the outward solution gave the final node count
	# continue the rest else pass
	if cal_nodes == N_0:
		# start the inward integration
		rf = R_INF(l,E,r,A,B,C,D,IOTP,kappa)
		r_inf,P_r0,Q_r0 = rf.r_infinity()
		prw = PQ_INWARD(r_inf,P_r0,Q_r0,r,A,B,C,D,IOTP,kappa,eta,E)
		Pr_inward,Qr_inward,P_inw,Q_inw =prw.pq_inward()
		# combine the outward and inward data
		all_data = INWARD_OUTWARD_DATA(r,Pr_outward,Qr_outward,\
								Pr_inward,Qr_inward,P_inw,Q_inw)
		RPQ_outward,RPQ_inward = all_data.r_pq_inward_outward()
		# rescaling inward with outward
		rescalen = RESCALING_INWARD(IOTP,Pr_outward,Qr_outward,P_out,Pr_inward,\
				   Qr_inward,P_inw)
		Pr_matched,Qr_matched = rescalen.pr_qr_matched()
		P_combined = Pr_outward+Pr_inward
		#
		nt = NORM_INTEGERAL(Pr_matched,Qr_matched,r)
		integ,R_new,A1,B1,C1,D1 = nt.normalization_integeral()
		#
		in_iv = INTEG_VALUE(RPQ_outward[1],RPQ_outward[2],RPQ_outward[0])
		outward_integ,A_out,B_out,C_out,D_out = in_iv.integeral_value()
		#
		out_iv = INTEG_VALUE(RPQ_inward[1],RPQ_inward[2],RPQ_inward[0])
		inward_integ,A_in,B_in,C_in,D_in = out_iv.integeral_value()
		# normalization
		pqnorm = PQ_NORMALIZED(Pr_matched,Qr_matched,integ)
		Pr_norm,Qr_norm = pqnorm.pq_normalized()
		#
		# If IOTP is the last point the output is 0.0
		# so we should stop here
		if integ != 0.0:
			#print(integ,P_out,Q_out,P_inw,Q_inw)
			den = ENERGY_CORRECTION(integ,P_out,Q_out,P_inw,Q_inw)
			DE = den.energy_correction()
			#
			egl = ENERGY_LOOP(E,EMIN,EMAX,DE,eta,N_0,A,B,C,D,r,rv,np,l,kappa)
			Eva,r_grid,Pr_norm_final,Qr_norm_final,Pr_match_f,Qr_match_f,\
			stab_final,integ_norm,expectation_r,expect_r2,\
				expect_inv_r = egl.energy_loop()
			###############################################
			##   Now check conversion occur 
			#################################################
			# end main calculations for this orbital
			print("\n")
			print(" Final Output orbital data (",label,"):")
			print(" < r >  = ",expectation_r)
			print(" < r^{2} >  = ",expect_r2)
			print(" < 1/r >  = ",expect_inv_r)
			# print the final energy in electron volts unit
			energy_eV = Hartree_to_eV*float(Eva)
			print(" The orbital energy = ",Eva," a.u.")
			print(" The orbital energy in eV unit is = ",energy_eV," eV.")
			########################
			# collect the data for this orbital; nothing is written to disk
			# yet -- ALL orbitals are combined into one file after the loop
			Ex = float(Eva)
			nnpx = len(r_grid)
			rx = r_grid
			prx = Pr_norm_final
			qrx  = Qr_norm_final
			#
			radial_data_all.append( (n_qn,kappa_qn,Ex,nnpx,rx,prx,qrx) )
			#
			summary_data.append( (label,n_qn,kappa_qn,Ex,energy_eV,\
								  expectation_r,expect_r2,pot_form,nnpx) )
			#
			if grasp_prt != "":
				# calculating the pz parameter
				zs = stab_final[0]
				a_0 = stab_final[2][0]
				rb = r_grid[1]
				#
				pzfn = PZ(integ_norm,a_0,rb,zs)
				pz = pzfn.pz_parameter()
				#
				grasp_data_all.append( (n_qn,kappa_qn,Ex,nnpx,pz,rx,prx,qrx) )
			#
			print("\n")
			print(" Orbital ",label," calculated (",orb_idx+1,"/",n_total,")")
		else:
			print("\n")
			print(" Orbital ",label," NOT converged")
			not_converged_data.append( (label,n_qn,kappa_qn,\
				"normalization integral = 0 (IOTP at last grid point "+\
				"or grid does not extend far enough)") )
	# if the outward gave error number of nodes 
	else:
		print("\n")
		print(" Orbital ",label," NOT converged")
		not_converged_data.append( (label,n_qn,kappa_qn,\
			"outward solution did not reach the correct node count "+\
			"(got "+str(cal_nodes)+", expected "+str(N_0)+")") )
#
########################
# write ALL orbitals into a single combined radial-wave-function file
########################
radial_file = "radial.rwf"
WRITE_DATA_MULTI(radial_file,Z,N_ele,radial_data_all).write()
print(" All ",n_total," orbital(s) saved to ",radial_file)
#
# write ALL orbitals into a single combined GRASP-format file
if grasp_prt == "":
	pass
else:
	grasp_file_output = grasp_prt
	WRITE_DATA_GRASP_MULTI(grasp_file_output,grasp_data_all).write_grasp()
	print(" All ",n_total," orbital(s) saved to ",grasp_file_output,\
		  " (GRASP format)")
#
print("\n")
print(" Finished ",n_total," orbital(s) calculated.")
#
########################
# summary of the whole run: orbital, energy, <r>, <r^2>, NP
########################
print("\n")
print(" ================================================")
print(" SUMMARY OF CALCULATION")
print(" ================================================")
summ = WRITE_SUMMARY("summary.dat",Z,N_ele,potential_label,\
		potential_constants,summary_data,not_converged_data)
summ.print_summary()
summ.write()
print(" summary saved to  summary.dat")


# END
