# Description

import sys	# needed? interfering???
from subprocess  import Popen

from fort_dp     import fort_dp

pcpv_location = '../PCPV/'

class Simmo(object):
	"""docstring for Simmo"""
	def __init__(self, solar_cell, light, other_para, 
		max_num_BMs, var_BM_min, max_order_PWs, p):
		self.solar_cell    = solar_cell
		self.light         = light
		self.other_para    = other_para
		self.max_num_BMs   = max_num_BMs
		self.var_BM_min    = var_BM_min
		self.max_order_PWs = max_order_PWs
		self.p             = p


	def inclusion_n(self):
		return self.solar_cell.inclusion.n(self.light.Lambda)

	def background_n(self):
		return self.solar_cell.background.n(self.light.Lambda)

	def superstrate_n(self):
		return self.solar_cell.superstrate.n(self.light.Lambda)

	def substrate_n(self):
		return self.solar_cell.substrate.n(self.light.Lambda)
		
	# def run(self):
	def fortran_command_str(self):
		"""Return a string that runs the simmo, to execute in a shell"""
		inclusion_n   = self.inclusion_n()
		background_n  = self.background_n()
		superstrate_n = self.superstrate_n()
		substrate_n   = self.substrate_n()
		if self.light.Lambda > self.var_BM_min:
			max_n         = self.solar_cell.inclusion.max_n()
			nval_tmp      = self.max_num_BMs*inclusion_n.real/max_n
			nval 		  = round(nval_tmp)
			ordre_ls      = round(self.max_order_PWs*nval/self.max_num_BMs)
		else:
			nval 		  = self.max_num_BMs
			ordre_ls      = self.max_order_PWs			


		command_to_run  = "%(pcpv_location)spcpv.exe %(parallel)d %(Lambda)s %(nval)d \
		%(ordre_ls)d %(d_in_nm)d %(debug)d %(mesh_file)s %(mesh_format)d \
		%(nb_typ_el)d %(re_n_eff_super)s %(im_n_eff_super)s \
		%(re_n_eff_sub)s %(im_n_eff_sub)s %(re_n_eff_bkg)s %(im_n_eff_bkg)s \
		%(re_n_eff_inc)s %(im_n_eff_inc)s %(has_substrate)d \
		%(theta)s %(phi)s %(h)s %(lx)s %(ly)s %(tol)s %(E_H_field)d \
		%(i_cond)d %(itermax)d %(pol)d %(traLambda)d %(PropModes)d \
		%(PrintSolution)d %(PrintSupModes)d %(PrintOmega)d %(PrintAll)d \
		%(Checks)d %(q_average)d %(plot_real)d %(plot_imag)d %(plot_abs)d \
		%(incident)d %(whatincident)d %(outincident)d %(Loss)d" % {
			'pcpv_location' : pcpv_location,
			'parallel' 	    : self.p,
			'nval'          : nval,
		    'ordre_ls'      : ordre_ls,
		    'h'             : fort_dp(self.solar_cell.height),
		    'lx'            : fort_dp(self.solar_cell.lx),
		    'ly'            : fort_dp(self.solar_cell.ly),
		    'd_in_nm'       : self.solar_cell.period,
		    'has_substrate' : self.solar_cell.has_substrate,
		    'mesh_file'     : self.solar_cell.mesh_file,
		    'mesh_format'   : self.solar_cell.mesh_format,
		    'nb_typ_el'     : self.solar_cell.nb_typ_el,
		    'Loss'          : self.solar_cell.loss,
		    're_n_eff_super': fort_dp(superstrate_n.real),
		    'im_n_eff_super': fort_dp(superstrate_n.imag),
		    're_n_eff_sub'  : fort_dp(substrate_n.real),
		    'im_n_eff_sub'  : fort_dp(substrate_n.imag),
		    're_n_eff_bkg'  : fort_dp(background_n.real),
		    'im_n_eff_bkg'  : fort_dp(background_n.imag),
		    're_n_eff_inc'  : fort_dp(inclusion_n.real),
		    'im_n_eff_inc'  : fort_dp(inclusion_n.imag),
			'Lambda'	    : fort_dp(self.light.Lambda),
		    'theta'         : fort_dp(self.light.theta),
		    'phi'           : fort_dp(self.light.phi),
		    'pol'           : self.light.pol,
		    'tol'           : fort_dp(self.other_para.tol),
		    'debug'         : self.other_para.debug,
		    'E_H_field'     : self.other_para.E_H_field,
		    'i_cond'        : self.other_para.i_cond,
		    'itermax'       : self.other_para.itermax,
		    'traLambda'     : self.other_para.traLambda,
		    'PropModes'     : self.other_para.PropModes,
		    'PrintSolution' : self.other_para.PrintSolution,
		    'PrintSupModes' : self.other_para.PrintSupModes,
		    'PrintOmega'    : self.other_para.PrintOmega,
		    'PrintAll'      : self.other_para.PrintAll,
		    'Checks'        : self.other_para.Checks,
		    'q_average'     : self.other_para.q_average,
		    'plot_real'     : self.other_para.plot_real,
		    'plot_imag'     : self.other_para.plot_imag,
		    'plot_abs'      : self.other_para.plot_abs,
		    'incident'      : self.other_para.incident,
		    'whatincident'  : self.other_para.what4incident,
		    'outincident'   : self.other_para.out4incident,
		}

		return command_to_run
