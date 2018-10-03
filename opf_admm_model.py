############################################################################
# This file builds the opf_admm_model class that represents a subproblem 
# ADMM algorithm parameters should be defined in admmoption
# Package Pypower 5.1.3 is used in this application
#############################################################################

from numpy import flatnonzero as find
from numpy import pi, Inf, ix_, ones, zeros, arange, append, unique, array, \
	reshape, vstack, transpose, dot, conj, exp, r_, maximum, minimum, absolute
from scipy.sparse import lil_matrix, vstack, hstack, csr_matrix as sparse


from math import ceil

from pypower.idx_brch import F_BUS, T_BUS, RATE_A
from pypower.idx_bus import BUS_AREA, BUS_TYPE, REF, VA, VMAX, VMIN
from pypower.idx_gen import PMAX, PMIN, QMAX, QMIN, PG, QG, GEN_BUS
from pypower.idx_cost import COST
from pypower.makeSbus import makeSbus
from pypower.dSbus_dV import dSbus_dV
from pypower.d2Sbus_dV2 import d2Sbus_dV2
from pypower.pips import pips
from copy import deepcopy



from pypower.dIbr_dV import dIbr_dV
from pypower.dSbr_dV import dSbr_dV
from pypower.dAbr_dV import dAbr_dV

from pypower.d2Sbus_dV2 import d2Sbus_dV2
from pypower.dSbr_dV import dSbr_dV
from pypower.dIbr_dV import dIbr_dV
from pypower.d2AIbr_dV2 import d2AIbr_dV2
from pypower.d2ASbr_dV2 import d2ASbr_dV2


from pypower.makeYbus import makeYbus

import numpy
import json

class opf_admm_model(object):
	"""This class encapsulates the local opf problem and implements admm steps
	including x-update(solving local opf), send data to neighbors, receive data
	from neighbors, z-update and y-update

	@author Junyao Guo
	"""

	def __init__(self):
		# initialize all the fields
		self.region = {
			'bus': None,
			'gen': None,
			'gencost': None,
			'branch':None,
			'Ybus': None,
			'Yf':None,
			'Yt':None,
			'nb': 0,
			'ng': 0,
			'nl':0,
			'ntl': 0,
			'nnbor': 0,
			'nwait': 0, #number of neighbors to wait for until next local update
			'Varefs': 0,
			'baseMVA':100
		}
        # local variables
		self.var = {
			'Vm': None, 'Va': None, 'Pg': None, 'Qg': None,
			'zmd': None, 'zms': None, 'zad': None, 'zas': None,
			'ymd': None, 'yms': None, 'yad': None, 'yas': None,
			'rho': None, 'AVmd': None, 'AVms': None, 'AVad': None, 'AVas': None
		}
        # get the indices
		self.idx = {
			'rbus': {'int': None, 'ext': None}, #local and global index of buses physically in the region
			'vbus': None, #global index of buses considered in the subproblem including copies of boundary buses 
			'gen': None, 'tl': [], 'fbus': [], 'tbus': [], #local indices of from and to buses of tielines
			'mapping': None, #idx mapping: key:global idx, value: local idx
			# idx for local variable x
			'var': {'iVa': None, 'iVm': None, 'iPg': None, 'iQg': None}
		}

		self.mtx = {
			'Adiff': None, 'Asum': None  # A matrix used in ADMM
		}

		# problem formulation for the pypower.pips solver
		self.pb = {
			'x0': None, 'A': None, 'l': None, 'u': None,
			'xmin': None, 'xmax': None, 'opt': None, 'solution': None,
			'f_fcn': None, 'gh_fcn': None, 'hess_fcn': None, 'converge': False,
			'primalgap': [Inf], 'dualgap': [Inf]
		}

		self.ID = None
		self.nbor = {}
		self.pipes = None
		self.admmopt = admmoption()
		self.recvmsg = {}
		self.gapAll = None 

    # configure the physical constraints 
	def config(self, ID, partition, bus, gen, gencost,branch, Ybus,Yf,Yt, genBus, tieline, pipes, na, baseMVA):
		self.region['baseMVA'] = baseMVA
		self.ID = ID
		self.gapAll = [10**8] * na
		self.idx['rbus']['ext'] = find(bus[:, BUS_AREA] == ID)
		self.idx['vbus'] = deepcopy(self.idx['rbus']['ext'])
		self.idx['gen'] = find(bus[ix_(genBus.astype(int), [BUS_AREA])] == ID).tolist()
		self.region['Varefs'] = bus[bus[:, BUS_TYPE] == REF, VA] * (pi / 180)
		print('Configuring worker %d ...' % (ID,))
		# configure neighbors
		for row in tieline:
			idx = find(row[2:4] == ID)
			if idx.tolist() != []:
				nid = row[3 - idx[0]].tolist()
				if nid not in self.nbor:
					self.nbor[nid] = neighbor()
				fb = row[idx].tolist()      #this is global idx
				tb = row[1 - idx].tolist()  #this is global idx
				lidx = row[-1].tolist()
				self.region['ntl'] += 1
				self.nbor[nid].config(fb[0], tb[0], lidx)
				self.idx['vbus'] = append(self.idx['vbus'], tb) #include neighboring buses to own region
		nbors = self.nbor.keys()
		self.region['nnbor'] = len(nbors)
		self.region['nwait'] = ceil(self.region['nnbor'] * self.admmopt.nwaitPercent)
		counter = 0

		# configure indices		
		self.idx['vbus'] = unique(self.idx['vbus']).tolist()
		self.idx['rbus']['ext'] = self.idx['rbus']['ext'].tolist()
		self.idx['mapping'] = {key: value for (value, key) in list(enumerate(self.idx['vbus']))}
		self.idx['rbus']['int'] = [self.idx['mapping'][i] for i in self.idx['rbus']['ext']]

		# set the local indices of tielines
		for nbor in nbors:
			tielines = self.nbor[nbor].tlidx['ext']
			idx = [counter + i for i in range(len(tielines))]
			self.nbor[nbor].tlidx['int'] += idx 
			counter = idx[-1] + 1
			self.idx['tl'] += tielines
			for t in range(len(tielines)):
				self.idx['fbus'] += [self.idx['mapping'][self.nbor[nbor].fbus[t]]]
				self.idx['tbus'] += [self.idx['mapping'][self.nbor[nbor].tbus[t]]]

		# get local buses, generator and costs
		self.region['bus'] = bus[self.idx['vbus'], :]
		self.region['gen'] = gen[self.idx['gen'], :] 
		self.region['gencost'] = gencost[self.idx['gen'], :] 
		self.region['nb'] = self.region['bus'].shape[0]
		self.region['ng'] = self.region['gen'].shape[0]
		self.region['Ybus'] = Ybus[:,self.idx['vbus']].tocsr()[self.idx['vbus'],:]
		self.region['Yf'] = Yf[:,self.idx['vbus']].tocsr()[self.idx['vbus'],:]
		self.region['Yt'] = Yt[:,self.idx['vbus']].tocsr()[self.idx['vbus'],:]
		


		# fetch the indices for variables
		st, nb, ng = 0, self.region['nb'], self.region['ng']
		iv = self.idx['var']
		iv['iVa'] = range(st, st + nb) 
		iv['iVm'] = range(st + nb, st + 2 * nb) 
		iv['iPg'] = range(st + 2 * nb, st + 2 * nb + ng) 
		iv['iQg'] = range(st + 2 * nb + ng, st + 2 * nb + 2 * ng) 
		
		
		
		
		
		tl1 = find(bus[branch[:, F_BUS].astype(int), BUS_AREA] == ID)
		tl2 = find(bus[branch[:, T_BUS].astype(int), BUS_AREA] == ID)
		
		tlall = list(set(tl1) | set(tl2)) 
		
		self.region['branch'] = branch[ix_(tlall, )]
		

		
		self.region['nl'] = self.region['branch'].shape[0]

		self.region['bus'][:,0] = [i for i in range(self.region['nb'])]
		for i in range(self.region['ng']):
			self.region['gen'][i,0] = self.idx['mapping'][self.region['gen'][i,0]]
			
		for i in range(self.region['nl']):
			self.region['branch'][i,0] = self.idx['mapping'][self.region['branch'][i,0]]
		for i in range(self.region['nl']):
			self.region['branch'][i,1] = self.idx['mapping'][self.region['branch'][i,1]]
		


		# configure communication pipes
		if pipes:
			self.pipes = pipes[ID] 

		# configure A matrix used in ADMM that indicates tie lines
		row = append(arange(self.region['ntl']), arange(self.region['ntl']))
		col = array(self.idx['fbus'] + self.idx['tbus'])
		data_diff = append(ones(self.region['ntl']) * self.admmopt.beta_diff, \
			ones(self.region['ntl']) * self.admmopt.beta_diff * -1)
		data_sum = append(ones(self.region['ntl']) * self.admmopt.beta_sum, \
			ones(self.region['ntl']) * self.admmopt.beta_sum)
		self.mtx['Adiff'] = sparse((data_diff, (row, col)), \
			shape = (self.region['ntl'], self.region['nb']))
		self.mtx['Asum'] = sparse((data_sum, (row, col)), \
			shape = (self.region['ntl'], self.region['nb']))

	# ----- initialize variables and set upper and lower bounds for variable -----
	def var_init(self):
		self.var['zmd'] = zeros(self.region['ntl'])
		self.var['zms'] = ones(self.region['ntl'])
		self.var['zad'] = zeros(self.region['ntl'])
		self.var['zas'] = zeros(self.region['ntl'])
		self.var['ymd'] = zeros(self.region['ntl'])
		self.var['yms'] = zeros(self.region['ntl'])
		self.var['yad'] = zeros(self.region['ntl'])
		self.var['yas'] = zeros(self.region['ntl'])
		self.var['rho'] = self.admmopt.rho_0
		if self.admmopt.init == 'flat':
			self.var['Va'] = self.region['Varefs'][0] * ones(self.region['nb'])
			self.var['Vm'] = (self.region['bus'][:, VMAX] + \
				self.region['bus'][:, VMIN]) / 2
			self.var['Pg'] = (self.region['gen'][:, PMAX] + \
				self.region['gen'][:, PMIN]) / 2
			self.var['Qg'] = (self.region['gen'][:, QMAX] + \
				self.region['gen'][:, QMIN]) / 2
		Adiff = self.mtx['Adiff']
		Asum = self.mtx['Asum']
		self.var['AVmd'] = Adiff * self.var['Vm']
		self.var['AVms'] = Asum * self.var['Vm']
		self.var['AVad'] = Adiff * self.var['Va']
		self.var['AVas'] = Asum * self.var['Va']

		Vamin, Vamax = -pi * ones(self.region['nb']), pi * ones(self.region['nb'])
		refbus = find(self.region['bus'][:, BUS_TYPE] == REF)
		if refbus.tolist() != []: # set the reference bus upper = lower bound
			Vamin[refbus], Vamax[refbus] = self.region['Varefs'], self.region['Varefs']
		Vmmin = self.region['bus'][:, VMIN]
		Vmmax = self.region['bus'][:, VMAX]
		Pgmin = self.region['gen'][:, PMIN]
		Pgmax = self.region['gen'][:, PMAX]
		Qgmin = self.region['gen'][:, QMIN]
		Qgmax = self.region['gen'][:, QMAX]
		self.pb['xmin'] = r_[Vamin, Vmmin, Pgmin, Qgmin]
		self.pb['xmax'] = r_[Vamax, Vmmax, Pgmax, Qgmax]
		self.pb['xmin'][self.pb['xmin'] == -Inf] = -1e10 #replace Inf with numerical proxies
		self.pb['xmax'][self.pb['xmax'] == Inf] = 1e10

	# ----- local opf solver for x-update using pypower.pips solver ------
	def pipslopf_solver(self):
		pb = self.pb 
		pb['x0'] = r_[self.var['Va'], self.var['Vm'], self.var['Pg'], self.var['Qg']]
		f_fcn = lambda x: self.admmopf_costfcn(x)
		gh_fcn = lambda x: self.admmopf_consfcn(x)
		
		hess_fcn = lambda x, lmbda, cost_mult: self.admmopf_hessfcn(x, lmbda, cost_mult)

		pb['solution'] = pips(f_fcn, pb['x0'], pb['A'], pb['l'], pb['u'], \
			pb['xmin'], pb['xmax'], gh_fcn, hess_fcn)
		# print("Subproblem %d is solved" % (self.ID,))
		return pb['solution']

	# ------ update local variables x according to self.pb['solution'] ------
	def update_x(self):
		iVa = self.idx['var']['iVa']
		iVm = self.idx['var']['iVm']
		iPg = self.idx['var']['iPg']
		iQg = self.idx['var']['iQg']
		Adiff = self.mtx['Adiff']
		Asum = self.mtx['Asum']
		if self.pb['solution']['x'].any():
			x = self.pb['solution']['x']
			self.var['Va'] = x[iVa]
			self.var['Vm'] = x[iVm]
			self.var['Pg'] = x[iPg]
			self.var['Qg'] = x[iQg]
			self.var['AVmd'] = Adiff * self.var['Vm']
			self.var['AVms'] = Asum * self.var['Vm']
			self.var['AVad'] = Adiff * self.var['Va']
			self.var['AVas'] = Asum * self.var['Va']
	
	# ----- evaluate cost function and gradient ----------
	def admmopf_costfcn(self, x = None):
		if x is None:
			x = self.pb['x0']
		nx = len(x)
		iv = self.idx['var']
		var = self.var
		# idx ranges
		iVa = iv['iVa']
		iVm = iv['iVm']
		iPg = iv['iPg']
		rho = self.var['rho']
		Adiff = self.mtx['Adiff']
		Asum = self.mtx['Asum']
		gencost = self.region['gencost']
		pg = x[iPg]
		vm = x[iVm]
		va = x[iVa]
		gc = dot(gencost[:,COST], pg ** 2) + dot(gencost[:,COST + 1], pg) + sum(gencost[:,COST + 2])
		mmd = Adiff * vm - var['zmd']
		mad = Adiff * va - var['zad']
		mms = Asum * vm - var['zms']
		mas = Asum * va - var['zas']
		# f is the augmented lagrangian 
		f = gc + var['ymd'].T.dot(mmd) + var['yad'].T.dot(mad) \
			+ var['yms'].T.dot(mms) + var['yas'].T.dot(mas) \
			+ rho / 2 * mmd.T.dot(mmd) + rho / 2 * mad.T.dot(mad) \
			+ rho / 2 * mms.T.dot(mms) + rho / 2 * mas.T.dot(mas) 
		# ----- evaluate cost gradient ------
		df = zeros(nx)
		df[iVa] = df[iVa] + sparse((var['yad'] + rho * mad).T).dot(Adiff) + \
			sparse((var['yas'] + rho * mas).T).dot(Asum)
		df[iVm] = df[iVm] + sparse((var['ymd'] + rho * mmd).T).dot(Adiff) + \
			sparse((var['yms'] + rho * mms).T).dot(Asum)
		df[iPg] = df[iPg] + 2 * pg.T * gencost[:, COST] + gencost[:, COST + 1]
		return f, df 

	# ---- evaluate nonlinear constraints and their gradients --------
	def admmopf_consfcn(self, x = None):
		if x is None:
			x = self.pb['x0']
		nx = len(x)
		nb = self.region['nb']
		ng = self.region['ng']
		iv = self.idx['var']
		bus = self.region['bus']
		gen = self.region['gen']
		branch = self.region['branch']
		Ybus = self.region['Ybus']
		Yf = self.region['Yf']
		Yt = self.region['Yt']
		baseMVA = self.region['baseMVA']
		
		ridx = self.idx['rbus']['int']
		# idx ranges
		iVa = iv['iVa']
		iVm = iv['iVm']
		iPg = iv['iPg']
		iQg = iv['iQg']
		# grab Pg and Qg
		gen[:, PG] = x[iPg]
		gen[:, QG] = x[iQg]
		# rebuid Sbus
		Sbus = makeSbus(1, bus, gen)
		# reconstruct V
		Va, Vm = x[iVa], x[iVm]
		V = Vm * exp(1j * Va)
		# evaluate power flow equations
		mis = V * conj(Ybus * V) - Sbus
		g = r_[mis.real, mis.imag]
		row = ridx +  [i + nb for i in ridx] 
		g = g[row]
		
		
		il = find((branch[:, RATE_A] != 0) & (branch[:, RATE_A] < 1e10))
		nl2 = len(il) 
		
		
		Ybus2, Yf2, Yt2 = makeYbus(baseMVA, bus, branch)
		Yf = Yf2[il,:]
		Yt = Yt2[il,:]
		
		if nl2 > 0:
			flow_max = (branch[il, RATE_A] / baseMVA)**2
			flow_max[flow_max == 0] = Inf
        
            ## compute branch power flows
            ## complex power injected at "from" bus (p.u.)
			Sf = V[ branch[il, F_BUS].astype(int) ] * conj(Yf * V)
            ## complex power injected at "to" bus (p.u.)
			St = V[ branch[il, T_BUS].astype(int) ] * conj(Yt * V)
			h = r_[ Sf.real**2 - flow_max,   ## branch P limits (from bus)
                        St.real**2 - flow_max ]  ## branch P limits (to bus)
		
		else: h = array([])
		
		
		
		# ---- evaluate constraint gradients -------
		# compute partials of injected bus powers
		dSbus_dVm, dSbus_dVa = dSbus_dV(Ybus, V) # w.r.t. V
		neg_Cg = sparse((-ones(ng), (gen[:, GEN_BUS], range(ng))), (nb, ng))
		

		
		
		# construct Jacobian of equality constraints (power flow) and transpose it
		dg = lil_matrix((2 * nb, nx))
		blank = sparse((nb, ng))
		dg = vstack([ \
			#P mismatch w.r.t Va, Vm, Pg, Qg
			hstack([dSbus_dVa.real, dSbus_dVm.real, neg_Cg, blank]),
			# Q mismatch w.r.t Va, Vm, Pg, Qg
			hstack([dSbus_dVa.imag, dSbus_dVm.imag, blank, neg_Cg])
			], "csr")
		dg = dg[row, :]
		dg = dg.T
		
		
		
		if nl2 > 0:
        ## compute partials of Flows w.r.t. V                  ## power
			dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft = \
                    dSbr_dV(branch[il, :], Yf, Yt, V)
			dFf_dVa = dFf_dVa.real
			dFf_dVm = dFf_dVm.real
			dFt_dVa = dFt_dVa.real
			dFt_dVm = dFt_dVm.real
			Ff = Ff.real
			Ft = Ft.real
			
			## squared magnitude of flow (of complex power or current, or real power)
			df_dVa, df_dVm, dt_dVa, dt_dVm = \
                dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft)
			## construct Jacobian of inequality constraints (branch limits)
			## and transpose it.
			dh = lil_matrix((2 * nl2, nx))
			dh[:, r_[iVa, iVm].T] = vstack([
                hstack([df_dVa, df_dVm]),    ## "from" flow limit
                hstack([dt_dVa, dt_dVm])     ## "to" flow limit
            ], "csr")
			dh = dh.T

        
		else: dh = None
		
		h = array([])
		dh = None
		
		return h, g, dh, dg

	def admmopf_hessfcn(self, x = None, lmbda = None , cost_mult=1.0):
		if x is None:
			x = self.pb['x0']
		if lmbda is None: # multiplier for power flow constraints
			lmbda = {'eqnonlin': 3000 * ones(2 * len(self.idx['rbus']['int']))}
			# lmbda = {'eqnonlin': 3000 * ones(2 * self.region['nb'])}
		nx = len(x)
		nb = self.region['nb']
		nl = self.region['nl']
		ng = self.region['ng']
		rho = self.var['rho']
		Ybus = self.region['Ybus']
		iv = self.idx['var']
		
		branch = self.region['branch']
		baseMVA = self.region['baseMVA']
		bus = self.region['bus']
		
		iVa = iv['iVa']
		iVm = iv['iVm']
		# reconstruct V
		Va, Vm = x[iVa], x[iVm]
		V = Vm * exp(1j * Va)
		# ---- evaluate Hessian of f ----------
		Adiff = self.mtx['Adiff']
		Asum = self.mtx['Asum']
		d2f_dVm2 = d2f_dVa2 = Adiff.T.dot(Adiff) * rho + Asum.T.dot(Asum) * rho 
		data = self.region['gencost'][:, COST] * 2
		i = arange(0, ng)
		d2f_dPgQg2 = sparse((data, (i, i)), (2 * ng, 2 * ng))
		blankv = sparse((nb, nb))
		d2f_dV2 = vstack([hstack([d2f_dVa2, blankv]), hstack([blankv, d2f_dVm2])], 'csr')
		blankvp = sparse((2 * nb, 2 * ng))
		d2f = vstack([hstack([d2f_dV2, blankvp]), hstack([blankvp.T, d2f_dPgQg2])], 'csr')
		
		# ---- evaluate Hessian of power flow constraints ------
		nlam = len(lmbda["eqnonlin"]) // 2
		lamP = zeros(nb)  #for non-included pf balances use 0 as multiplier
		lamQ = zeros(nb)
		lamP[self.idx['rbus']['int']] = lmbda["eqnonlin"][:nlam]
		lamQ[self.idx['rbus']['int']] = lmbda["eqnonlin"][nlam:nlam + nlam]
		# lamP = lmbda["eqnonlin"][:nlam]
		# lamQ = lmbda["eqnonlin"][nlam:nlam + nlam]
		Gpaa, Gpav, Gpva, Gpvv = d2Sbus_dV2(Ybus, V, lamP)
		Gqaa, Gqav, Gqva, Gqvv = d2Sbus_dV2(Ybus, V, lamQ)

		
		d2G = vstack([
				hstack([
					vstack([hstack([Gpaa, Gpav]),
							hstack([Gpva, Gpvv])]).real +
					vstack([hstack([Gqaa, Gqav]),
							hstack([Gqva, Gqvv])]).imag,
					sparse((2 * nb, 2 * ng))]),
				hstack([
					sparse((2 * ng, 2 * nb)),
					sparse((2 * ng, 2 * ng))
				])
			], 'csr')
		
		
		##----- evaluate Hessian of flow constraints -----
		"""
		il = find((branch[:, RATE_A] != 0) & (branch[:, RATE_A] < 1e10))
		nl2 = len(il) 
		
		
		Ybus2, Yf2, Yt2 = makeYbus(baseMVA, bus, branch)
		Yf = Yf2[il,:]
		Yt = Yt2[il,:]
		
		
		
		if il is None:
			il = arange(nl)            ## all lines have limits by default
		nl2 = len(il) ## number of constrained lines
		
		
		nmu = int(len(lmbda["ineqnonlin"]) / 2)
		muF = lmbda["ineqnonlin"][:nmu]
		muT = lmbda["ineqnonlin"][nmu:nmu + nmu]



		f = branch[il, F_BUS].astype(int)    ## list of "from" buses
		t = branch[il, T_BUS].astype(int)    ## list of "to" buses
		## connection matrix for line & from buses
		Cf = sparse((ones(nl2), (arange(nl2), f)), (nl2, nb))
		## connection matrix for line & to buses
		Ct = sparse((ones(nl2), (arange(nl2), t)), (nl2, nb))
		dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St = \
                dSbr_dV(branch[il,:], Yf, Yt, V)

		Hfaa, Hfav, Hfva, Hfvv = d2ASbr_dV2(dSf_dVa.real, dSf_dVm.real,
                                                Sf.real, Cf, Yf, V, muF)
		Htaa, Htav, Htva, Htvv = d2ASbr_dV2(dSt_dVa.real, dSt_dVm.real,
                                                St.real, Ct, Yt, V, muT)
		nxtra = nx - 2 * nb
		d2H = vstack([
            hstack([
                vstack([hstack([Hfaa, Hfav]),
                        hstack([Hfva, Hfvv])]) +
                vstack([hstack([Htaa, Htav]),
                        hstack([Htva, Htvv])]),
                sparse((2 * nb, nxtra))
            ]),
            hstack([
                sparse((nxtra, 2 * nb)),
                sparse((nxtra, nxtra))
            ])
        ], "csr")
		
		
		"""
		
		
		

		return d2f + d2G #+d2H

	# update z according to received message
	def update_z(self):
		srcs = self.nbor.keys() # all neighbors
		var = self.var
		zmd_old, zms_old, zad_old, zas_old = deepcopy(var['zmd']), deepcopy(var['zms']),\
			deepcopy(var['zad']), deepcopy(var['zas'])
		for k in srcs:
			tlidx = self.nbor[k].tlidx['int'] # internal index of tie lines to k
			if k not in self.recvmsg:  #if neighbor k has not arrived
				var['zmd'][tlidx] = var['AVmd'][tlidx]
				var['zms'][tlidx] = var['AVms'][tlidx]
				var['zad'][tlidx] = var['AVad'][tlidx]
				var['zas'][tlidx] = var['AVas'][tlidx]
			else:
				if self.recvmsg[k].tID != self.ID or \
					len(self.recvmsg[k].fields['ymd']) != len(tlidx):
					print('The received message from %d does not match region %d !' % (k, self.ID))
					continue
				else:
					nborvar = self.recvmsg[k].fields
					print(self.ID, nborvar['convergeTable'], nborvar['rho'])
					var['zmd'][tlidx] = (var['ymd'][tlidx] - nborvar['ymd'] + var['rho'] * var['AVmd'][tlidx] - \
						nborvar['rho'] * nborvar['AVmd']) / (var['rho'] + nborvar['rho'])
					var['zms'][tlidx] = (var['yms'][tlidx] + nborvar['yms'] + var['rho'] * var['AVms'][tlidx] + \
						nborvar['rho'] * nborvar['AVms']) / (var['rho'] + nborvar['rho'])
					var['zad'][tlidx] = (var['yad'][tlidx] - nborvar['yad'] + var['rho'] * var['AVad'][tlidx] - \
						nborvar['rho'] * nborvar['AVad']) / (var['rho'] + nborvar['rho'])
					var['zas'][tlidx] = (var['yas'][tlidx] + nborvar['yas'] + var['rho'] * var['AVas'][tlidx] + \
						nborvar['rho'] * nborvar['AVas']) / (var['rho'] + nborvar['rho'])
					# print('z updated at %d using messages received from %d !' % (self.ID, k))
		# update dual gap
		self.pb['dualgap'] += [max(absolute(r_[zmd_old - var['zmd'], zms_old - var['zms'], \
			zad_old - var['zad'], zas_old - var['zas']]))]

	# update y (multiplier associated with Ax-z = 0)
	def update_y(self):
		var = self.var
		rho = self.var['rho'] 
		N = len(var['ymd'])
		ymax, ymin = self.admmopt.ymax, self.admmopt.ymin
		var['ymd'] += rho * (var['AVmd'] - var['zmd'])
		var['yms'] += rho * (var['AVms'] - var['zms'])
		var['yad'] += rho * (var['AVad'] - var['zad'])
		var['yas'] += rho * (var['AVas'] - var['zas'])
		var['ymd'] = maximum(ymin * ones(N), minimum(var['ymd'], ymax * ones(N)))
		var['yms'] = maximum(ymin * ones(N), minimum(var['yms'], ymax * ones(N)))
		var['yad'] = maximum(ymin * ones(N), minimum(var['yad'], ymax * ones(N)))
		var['yas'] = maximum(ymin * ones(N), minimum(var['yas'], ymax * ones(N)))

	# update rho and primal gap locally
	def update_rho(self):
		var = self.var
		# calculate and update primal gap first
		primalgap_old = deepcopy(self.pb['primalgap'][-1])
		self.pb['primalgap'] += [max(absolute(r_[var['AVmd'] - var['zmd'], var['AVms'] - var['zms'], \
			var['AVad'] - var['zad'], var['AVas'] - var['zas']]))]
		# update rho if necessary
		if self.pb['primalgap'][-1] > self.admmopt.theta * primalgap_old:
			var['rho'] *= self.admmopt.tau
			var['rho'] = minimum(var['rho'], self.admmopt.rho_max)
		# update local converge table
		self.gapAll[self.ID - 1] = self.pb['primalgap'][-1]

	# use the maximum rho among neighbors for local update		
	def choose_max_rho(self):
		srcs = self.recvmsg.keys()
		for k in srcs:
			rho_nbor = self.recvmsg[k].fields['rho']
			self.var['rho'] = maximum(self.var['rho'], rho_nbor)

	def send(self):
		#dest = self.pipes.keys()
		dest = self.nbor.keys()
		#print('LALALALALALALA')
		#print(dest)
		#print(self.nbor.keys())
		for k in dest:
			# prepare the message to be sent to neighbor k
			msg = message() 
			#print(k)
			
			msg.config(self.ID, k, self.var, self.nbor[k].tlidx['int'], self.gapAll)
			self.pipes[k].send(msg)

			# print("Message sent from %d to %d" % (self.ID, k))
			# print(msg)

	def recv(self):
		twait = self.admmopt.pollWaitingtime
		dest = list(self.pipes.keys())
		#print(self.nbor.keys() ,"comp", dest)
		
		
		#print(dest)
		recvFlag = [0] * self.region['nnbor']
		arrived = 0 # number of arrived neighbors
		pollround = 0
		
		# keep receiving from nbor 1 to nbor K in round until nwait neighbors arrived
		while arrived < self.region['nwait'] and pollround < 5: 
			for i in range(len(dest)):
				k = dest[i]
				
				while self.pipes[k].poll(twait): #read from pipe until get the last message
					self.recvmsg[k] = self.pipes[k].recv()
					print("APO", self.recvmsg[k].fID , "PROS", self.recvmsg[k].tID )
					"""
					print("STO RECVVVVV %d\n" %self.recvmsg[k].fID)
					print(type(self.recvmsg[k].fields['AVad']))
					print(type(self.recvmsg[k].fields['yad']))
					print(type(self.recvmsg[k].fields['rho']))
					print(type(self.recvmsg[k].fields['convergeTable']))
					"""
					#print("ME PIPES", (self.recvmsg[k].fields['convergeTable']))
					recvFlag[i] = 1
					# print("Message received at %d from %d" % (self.ID, k))
			arrived = sum(recvFlag)		
			pollround += 1

	def converge(self):
		# first update local converge table using received converge tables
		if self.recvmsg is not None:
			for k in self.recvmsg:
				table = self.recvmsg[k].fields['convergeTable'] 
				self.gapAll = list(map(min, zip(self.gapAll, table)))
		# check if all local primal gaps < tolerance
		if (max(self.gapAll) < self.admmopt.convergetol and  self.pb['dualgap'][-1] < self.admmopt.convergetol):
			return True
		else:
			return False		

class neighbor(object):
	""" This class encapsulates all related info regarding a neighbor of 
	a region
	"""
	def __init__(self):
		self.fbus = []
		self.tbus = []
		self.tlidx = {'ext':[], 'int':[]}

	def config(self, fb, tb, lidx):
		self.fbus.append(fb)
		self.tbus.append(tb)
		self.tlidx['ext'].append(lidx)

##--------ADMM parameters specification -------------------------------------
class admmoption(object):
	""" This class defines all the parameters to use in admm """

	def __init__(self):
		self.beta_diff = 2    # weight for (Vi - Vj) of tie line
		self.beta_sum = 0.5   # weight for (Vi + Vj) of tie line
		self.rho_0 = 1*10**4  # initial value for penalty rho, at the same order of initial cost function
		self.rho_max = 10**16 # upper bound for penalty rho
		self.tau = 1.02    # multiplier for increasing rho
		self.theta = 0.99 # multiplier for determining whether to update rho
		self.init = 'flat'    # starting point
		self.ymax = 10**16    # maximum y
		self.ymin = -10**16   # minimum y
		self.pollWaitingtime = 0.001 # waiting time of receiving from one pipe
		self.nwaitPercent = 0.1  # waiting percentage of neighbors (0, 1]
		self.iterMaxlocal = 2000 # local maximum iteration
		self.convergetol = 10**(-3) # convergence criteria for maximum primal gap

class message(object):
	""" This class defines the message region i sends to/receives from j """

	def __init__(self):
		self.fID = 0 #source region ID
		self.tID = 0 #destination region ID
		self.fields = {
			'AVmd': None, 'AVms': None, 'AVad': None, 'AVas': None,
			'ymd': None, 'yms': None, 'yad': None, 'yas': None,
			'rho': None, 'convergeTable': None
		}
        
        

	def config(self, f, t, var, tlidx, gapAll): #AVall and var are local variables of f region
		self.fID = f 
		self.tID = t 
		self.fields['AVmd'] = var['AVmd'][tlidx]
		self.fields['AVms'] = var['AVms'][tlidx]
		self.fields['AVad'] = var['AVad'][tlidx]
		
		self.fields['AVas'] = var['AVas'][tlidx]
		self.fields['ymd'] = var['ymd'][tlidx]
		self.fields['yms'] = var['yms'][tlidx]
		self.fields['yad'] = var['yad'][tlidx]
		self.fields['yas'] = var['yas'][tlidx]
		
		self.fields['rho'] = var['rho']
		
		self.fields['convergeTable'] = gapAll
		"""
		print("LALAAL %d\n" %f)
		print(type(self.fields['AVad']))
		print(type(self.fields['yad']))
		print(type(self.fields['rho']))
		print(type(self.fields['convergeTable']))
		"""
		
	def toJSON(self):
		return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


    

class message2(object):
	""" This class defines the message region i sends to/receives from j """

	
	def __init__(self):
		self.fID = 0 #source region ID
		self.tID = 0 #destination region ID
		self.fields = {
			'Pg': None, 'Qg': None, 'Vm': None, 'Va': None,
			 'AVmd': None, 'AVms': None, 'AVad': None, 'AVas': None,
			'ymd': None, 'yms': None, 'yad': None, 'yas': None,
			'rho': None, 'convergeTable': None
		}
        
	def config(self, f, t, var, tlidx, gapAll): #AVall and var are local variables of f region
	
		self.fID = f  
		self.fields['Pg'] = var['Pg'].tolist()
		self.fields['Qg'] = var['Qg'].tolist()
		self.fields['Vm'] = var['Vm'].tolist()
		self.fields['Va'] = var['Va'].tolist()

		self.tID = t 
		self.fields['AVmd'] = var['AVmd'][tlidx].tolist()
		self.fields['AVms'] = var['AVms'][tlidx].tolist()
		self.fields['AVad'] = var['AVad'][tlidx].tolist()
		self.fields['AVas'] = var['AVas'][tlidx].tolist()
		self.fields['ymd'] = var['ymd'][tlidx].tolist()
		self.fields['yms'] = var['yms'][tlidx].tolist()
		self.fields['yad'] = var['yad'][tlidx].tolist()
		self.fields['yas'] = var['yas'][tlidx].tolist()
		self.fields['rho'] = var['rho']
		self.fields['convergeTable'] = gapAll
		
class message3(object):
	""" This class defines the message region i sends to/receives from j """

	def __init__(self):
		self.fID = 0 #source region ID
		self.tID = 0 #destination region ID
		self.fields = {
			'AVmd': None, 'AVms': None, 'AVad': None, 'AVas': None,
			'ymd': None, 'yms': None, 'yad': None, 'yas': None,
			'rho': None, 'convergeTable': None
		}
        
        

	def config(self, f, t, var, tlidx, gapAll): #AVall and var are local variables of f region
		self.fID = f 
		self.tID = t 
		self.fields['AVmd'] = var['AVmd'][tlidx].tolist()
		self.fields['AVms'] = var['AVms'][tlidx].tolist()
		self.fields['AVad'] = var['AVad'][tlidx].tolist()
		self.fields['AVas'] = var['AVas'][tlidx].tolist()
		self.fields['ymd'] = var['ymd'][tlidx].tolist()
		self.fields['yms'] = var['yms'][tlidx].tolist()
		self.fields['yad'] = var['yad'][tlidx].tolist()
		self.fields['yas'] = var['yas'][tlidx].tolist()
		self.fields['rho'] = var['rho']
		self.fields['convergeTable'] = gapAll
