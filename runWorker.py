#############################################################
# This functions implements the ADMM iterations at one worker
# This functions can be run in parallel using the Process Module
#############################################################
from sys import stdout
from eth_rpc_client import Client
import sys

import web3
import matplotlib.pyplot as plt
from ethjsonrpc import EthJsonRpc
import matplotlib.pyplot as plt

from numpy import \
    ones, zeros, r_, sort, exp, pi, diff, arange, min, \
    argmin, argmax, logical_or, real, imag, any


from numpy import flatnonzero as find

from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF
from pypower.idx_gen import GEN_BUS, PG, QG, QMAX, QMIN, GEN_STATUS, \
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN
from pypower.idx_brch import F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, \
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX

from pypower.isload import isload
from pypower.run_userfcn import run_userfcn
from pypower.ppoption import ppoption
from opf_admm_model import message2, message3
from opf_admm_model import message
from pypower.int2ext import int2ext
from pypower.ext2int import ext2int
from pypower.makeYbus import makeYbus
# from opf_admm_model import opf_admm_model
from time import time
from numpy import dot
from pypower.idx_cost import COST
from pypower.api import case9, ppoption, runpf, printpf, runopf, rundcopf
import time as tm

from numpy import conj, zeros, c_, shape, ix_


import json
import binascii
import numpy
import copy
import codecs

from json import JSONEncoder


class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

class MyDecoder(json.JSONDecoder):
    def default(self, obj):
        if isinstance(obj, int):
            return numpy.integer(obj)
        elif isinstance(obj, float):
            return numpy.floating(obj)
        #elif isinstance(obj, numpy.ndarray):
            #return obj.tolist()
        else:
            return super(MyDecoder, self).default(obj)



def runWorker(ID, s, output):

	client = EthJsonRpc('10.64.83.200', 39000+ID)
	obj =  web3.Web3(web3.HTTPProvider('http://10.64.83.200:3900'+str(ID)))
	print(client.eth_accounts()[0])
		#print(client1.eth_accounts()[1])
		


		
        
	print("Worker %d initialized successfully!" % (ID,))
	
	all_gas=[]
	nu = 0 #iteration count
	itermax = s.admmopt.iterMaxlocal #get maximum iteration
	flag = False
	new_gas=0
	
	
	
	start_balance = obj.eth.getBalance(obj.toChecksumAddress(client.eth_accounts()[0]))
	
	
	
	
	j=0
	s.admmopt.tau = 1
	last_block=obj.eth.blockNumber
	while (1):
		trans_by_block=obj.eth.getBlockTransactionCount(last_block - j)
		if(trans_by_block>0):
			#dest = list(s.nbor.keys())
			#if(ID==1):getNextValues(ID,s,client,obj,last_block - j, trans_by_block)
			
			#s.update_z()
			#s.choose_max_rho()
			break
			
		j=j+1
	
	
	prev_convegenceTable=[]
	P_history=[]
	Q_history=[]
	start_time2 = time()
	
	while nu <= itermax and not flag:
		if s.recvmsg: 
			s.update_z()
			s.choose_max_rho()
			
			#print(ID,"i have rho", s.var["rho"])
			
			
		start_time = time()
		result = s.pipslopf_solver()
		
		end_time = time()

		if result['eflag'] and nu % 20 == 0:
			print('Subproblem %d at iteration %d solved!' % (ID, nu) )
			# print("Time for solving subproblem %d: %ssecs to %ssecs" % (ID, start_time, end_time))

		s.update_x()
		P_history.append(sum(s.var['Pg']))
		Q_history.append(sum(s.var['Qg']))
		

		if s.recvmsg:  # not the initialization
			s.update_y()
			s.update_rho()
			
			prev_convegenceTable = s.gapAll[s.ID - 1] 

		# check convergence
		#if(nu>10):flag = 
		flag =s.converge()
		s.recvmsg = {}  # clear the received messages 
		
		#s.send()
		
		
		dest = s.nbor.keys()
		
		for k in dest:
			# prepare the message to be sent to neighbor k
			#tm.sleep(0.05)
			msg3 = message() 
			msg3.config(s.ID, k, s.var, s.nbor[k].tlidx['int'], s.gapAll)
			data2=json.dumps(msg3.__dict__, cls=MyEncoder)
			
			json_data =json.JSONEncoder().encode(data2)
			
			
			
			json_data ='0x'+json_data.encode("utf-8").hex()
			
			
			sampleDict = obj.txpool.content.pending
			
			new_gas= new_gas +18000000000* obj.eth.estimateGas({'to': obj.toChecksumAddress( '0xa8085d8331f16a5b76690e665d9f5eaaaa85ee1c') ,'data' :json_data })
			
			
			
			
			
			client.eth_sendTransaction(from_address= client.eth_accounts()[0],
                      to_address="0xa8085d8331f16a5b76690e665d9f5eaaaa85ee1c",
                      gas=  0xf4240,
                      gas_price=18000000000,
                      value= 1, # 2441406250
                      data = json_data
                      )
			
			
		
		txpooljson = obj.txpool.content.pending.__dict__
		
		recvmsg={}
		
		twait = s.admmopt.pollWaitingtime
		dest = list(s.nbor.keys())
		recvFlag = [0] * s.region['nnbor']
		arrived = 0 # number of arrived neighbors
		pollround = 0
		while arrived < s.region['nwait'] and pollround < 5: 
			#for i in range(len(dest)):
				#k = dest[i]
				#s.recvmsg[k] = recvmsg[k]
				#recvFlag[i] = 1
		
				for k in txpooljson.keys():
                    
					txpooljson2 = txpooljson[k].__dict__
					last_nonce = max(txpooljson2, key=int)
					#print(ID,", last nonce",last_nonce)
					last_nonce_dict = txpooljson2[last_nonce].__dict__
			
					hexjson = last_nonce_dict['input'][2:]
					jsonstring = codecs.decode(hexjson, "hex").decode('utf-8')
	
					jsonvalues = json.JSONDecoder().decode(jsonstring)
					values_dict = json.loads(jsonvalues, cls=MyDecoder)

			
					temp_msg = message()
					if( values_dict['fID'] in s.nbor.keys()):
						for last_nonce in range( int(max(txpooljson2, key=int)),0,-1):
							last_nonce_dict = txpooljson2[str(last_nonce)].__dict__
			
							hexjson = last_nonce_dict['input'][2:]
							jsonstring = codecs.decode(hexjson, "hex").decode('utf-8')
	
							jsonvalues = json.JSONDecoder().decode(jsonstring)
							values_dict = json.loads(jsonvalues, cls=MyDecoder)
							if( values_dict['tID'] == s.ID):break
			
			
			
					#print(ID,"last nonce=",last_nonce,"from",values_dict['fID'])
					temp_msg.fID = values_dict['fID']
					temp_msg.tID = values_dict['tID']
					temp_msg.fields['AVmd'] = numpy.asarray(values_dict['fields']['AVmd'])
					temp_msg.fields['AVms'] = numpy.asarray(values_dict['fields']['AVms'])
					temp_msg.fields['AVad'] = numpy.asarray(values_dict['fields']['AVad'])
					temp_msg.fields['AVas'] = numpy.asarray(values_dict['fields']['AVas'])
					temp_msg.fields['ymd'] = numpy.asarray(values_dict['fields']['ymd'])
					temp_msg.fields['yms'] = numpy.asarray(values_dict['fields']['yms'])
					temp_msg.fields['yad'] = numpy.asarray(values_dict['fields']['yad'])
					temp_msg.fields['yas'] = numpy.asarray(values_dict['fields']['yas'])
					temp_msg.fields['rho'] = values_dict['fields']['rho']
					temp_msg.fields['convergeTable'] = values_dict['fields']['convergeTable']
			
			
			
					if(temp_msg.tID == s.ID):
						recvmsg[temp_msg.fID] = temp_msg
						#recvFlag[i] = 1
					
					arrived = len(recvmsg)		
					pollround += 1
			

		
		
		s.recvmsg = copy.deepcopy(recvmsg)

		all_gas.append(new_gas)
		nu += 1

	# record results 
	print("Worker %d finished!" % (ID,))
	
	
	
	for k in dest:
		
		starting_point = message2() 
		starting_point.config(s.ID, k, s.var, s.nbor[k].tlidx['int'], s.gapAll)
		data2=json.dumps(starting_point.__dict__, cls=MyEncoder)
		json_data =json.JSONEncoder().encode(data2)
		json_data ='0x'+json_data.encode("utf-8").hex()
		sampleDict = obj.txpool.content.pending
	
	
	
	
		client.eth_sendTransaction(from_address= client.eth_accounts()[0],
                      to_address="0xa8085d8331f16a5b76690e665d9f5eaaaa85ee1c",
                      gas=  0xf4240,
                      gas_price=18000000000,
                      value= 1, # 2441406250
                      data = json_data
                      )
	
	
	
	print(starting_point.__dict__)
	
	x, f, info, lmbda, output2 = result["x"], result["f"], result["eflag"], result["lmbda"], result["output"]
	nx = len(x)
	nb = s.region['nb']
	ng = s.region['ng']
	iv = s.idx['var']
	bus = s.region['bus']
	gen = s.region['gen']
	branch = s.region['branch']
	baseMVA = s.region['baseMVA']
	Ybus = s.region['Ybus']
	ridx = s.idx['rbus']['int']
	# idx ranges
	iVa = iv['iVa']
	iVm = iv['iVm']
	iPg = iv['iPg']
	iQg = iv['iQg']
	# grab Pg and Qg
	gen[:, PG] = x[iPg] * s.region["baseMVA"]
	gen[:, QG] = x[iQg] * s.region["baseMVA"]
	bus[:, PD] = bus[:, PD] * s.region["baseMVA"]
	bus[:, QD] = bus[:, QD] * s.region["baseMVA"]
	# reconstruct V
	Va, Vm = x[iVa], x[iVm]
	V = Vm * exp(1j * Va)
	
	#print(V)
	nl = shape(branch)[0] ## number of branches
	bus[:, VA] = Va * 180 / pi
	bus[:, VM] = Vm
	
	if shape(branch)[1] < MU_ANGMAX + 1:
		branch = c_[branch, zeros((nl, MU_ANGMAX + 1 - shape(branch)[1]))]
	
	
	Ybus2, Yf, Yt = makeYbus(baseMVA, bus, branch)
	#print(Yf)
	## compute branch flows
	Sf = V[ branch[:, F_BUS].astype(int) ] * conj(Yf * V)  ## cplx pwr at "from" bus, p["u"].
	St = V[ branch[:, T_BUS].astype(int) ] * conj(Yt * V)  ## cplx pwr at "to" bus, p["u"].
	
	
	branch[:, PF] = Sf.real * baseMVA
	branch[:, QF] = Sf.imag * baseMVA
	branch[:, PT] = St.real * baseMVA
	branch[:, QT] = St.imag * baseMVA
	
	

	
	
	#gen[:, VG] = Vm[ gen[:, GEN_BUS].astype(int) ]
	
	nlam = len(lmbda["eqnonlin"]) // 2
	lamP = zeros(nb)  #for non-included pf balances use 0 as multiplier
	lamQ = zeros(nb)
	lamP[s.idx['rbus']['int']] = lmbda["eqnonlin"][:nlam] / s.region["baseMVA"]
	lamQ[s.idx['rbus']['int']] = lmbda["eqnonlin"][nlam:nlam + nlam] / s.region["baseMVA"]
	
	ong = find( (gen[:, GEN_STATUS] > 0) & ~isload(gen) )
	objValue1=0
	objValue2=0
	
	

    
	fd = stdout
	fd.write('\n REGION %d' % (s.ID) )
	
	fd.write('\nBus/Area  Voltage          Generation             Load        ')
	fd.write('  Lambda($/MVA-hr)')
	fd.write('\n  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)')
	fd.write('     P        Q   ')
	fd.write('\n----- ------- --------  --------  --------  --------  --------')
	fd.write('  -------  -------')
	for i in range(nb):
		for glob, loc in s.idx["mapping"].items():
			if loc == i:
				busid=glob+1                
				#bus[i,BUS_I]=glob+1
				pass
		fd.write('\n%5d/' %busid)
		fd.write('%d%7.3f%9.3f' %tuple( bus[i,[ BUS_AREA, VM, VA]]))
		if bus[i, BUS_TYPE] == REF:
			fd.write('*')
		else:
			fd.write(' ')
		g  = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) & ~isload(gen))
		ld = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) & isload(gen))
		if any(g + 1):
			fd.write('%9.2f%10.2f' % (sum(gen[g, PG]), sum(gen[g, QG])))
			objValue1 = objValue1 + (lamP[i] ) * (sum(gen[g, PG]))
		else:
			fd.write('      -         -  ')
		if logical_or(bus[i, PD], bus[i, QD]) | any(ld + 1):
			if any(ld + 1):
				fd.write('%10.2f*%9.2f*' % (bus[i, PD] - sum(gen[ld, PG]), bus[i, QD] - sum(gen[ld, QG])))
				objValue2 = objValue2 + (lamP[i]) * (bus[i, PD] - sum(gen[ld, PG]))
			else:
				fd.write('%10.2f%10.2f ' % tuple(bus[i, [PD, QD]]))
		else:
			fd.write('       -         -   ')
            
		fd.write('%9.3f' % lamP[i])
		
		if abs( lamQ[i] ) >  1e-4 :
			fd.write('%8.3f' % lamQ[i])
		else:
			fd.write('     -')
	fd.write('\n                        --------  --------  --------  --------')
	nzld = find((bus[:, PD] != 0.0) | (bus[:, QD] != 0.0))
	onld = find( (gen[:, GEN_STATUS] > 0) & isload(gen) )
	
	fd.write('\n               Total: %9.2f %9.2f %9.2f %9.2f' %
            (sum(gen[ong, PG]), sum(gen[ong, QG]),
             sum(bus[nzld, PD]) - sum(gen[onld, PG]),
             sum(bus[nzld, QD]) - sum(gen[onld, QG])))
	fd.write('\n')
	
	
	print("Local iteration of worker %d is %d" % (ID, nu))
	# calculate local generation cost
	gencost = s.region['gencost']
	pg = s.var['Pg']
	
	"""
	objValue21=0
	for i in range(ng):
		if(pg[i]>0):
			objValue1 = objValue1 - gencost[i,COST+1]* pg[i]
			print(gencost[i,COST+1]* pg[i])
		else:
			objValue21 = objValue21 + gencost[i,COST+1]* (-pg[i])
			print(gencost[i,COST+1]* (-pg[i]))
	objValue2 = objValue21 - objValue2
	objValue = objValue1 + objValue2
	"""
	
	objValue = dot(gencost[:,COST], pg ** 2) + dot(gencost[:,COST + 1], pg) + sum(gencost[:,COST + 2])
	print(objValue)
	

	varDual = {'ymd': s.var['ymd'], 'yad': s.var['yad'], 'yms': s.var['yms'],
		'yas': s.var['yas']}
	varPrimal = {'Vm': s.var['Vm'], 'Va': s.var['Va'], 
		'Pg': s.var['Pg'], 'Qg': s.var['Qg']}
	Result = {'objValue': objValue, 'varPrimal': varPrimal, 'varDual': varDual, 'localiter': nu, 
		'primalgap': s.pb['primalgap'], 'dualgap': s.pb['dualgap']}
	
	
	
	
	
	
	
	ng = s.region['ng']
	for i in range(ng):
		for glob, loc in s.idx["mapping"].items():
			if gen[i,GEN_BUS] == loc:
				gen[i,GEN_BUS]=glob+1
				break
				
				
				
	for i in range(nb):
		for glob, loc in s.idx["mapping"].items():
			if loc == i:
				
				bus[i,BUS_I]=glob+1
				
				
				
				
	nl = s.region['nl']
	for i in range(nl):
		for glob, loc in s.idx["mapping"].items():
			#print(glob, loc, branch[i, F_BUS])
			if branch[i, F_BUS] == loc:
				#print(branch[tl1,F_BUS])
				branch[i,F_BUS]=glob+1
				break
	for i in range(nl):
		for glob, loc in s.idx["mapping"].items():
			#print(glob, loc, branch[i, F_BUS])
			if branch[i, T_BUS] == loc:
				#print(branch[tl1,F_BUS])
				branch[i,T_BUS]=glob+1
				break
	
 
				 
				 
				 



	
	success=flag
	reg_nb = find(bus[:, BUS_AREA] == ID)
	
	lamP = lamP[ix_(reg_nb, )]
	
	lamP2= numpy.array(lamP)
	
	reg_bus = bus[ix_(reg_nb, )]
	print(reg_bus[:,BUS_I])
	reg_lam = numpy.array((reg_bus[:,BUS_I], lamP))
	
	print(reg_lam)
	
	
	ppc={}
	ppc["bus"], ppc["gen"], ppc["branch"] = reg_bus, gen, branch
	
	ids1=reg_bus[:,BUS_I]
	ids=numpy.array(ids1)
	#ppc["success"] = success
	#ppc=ext2int(ppc)
	P_history = [i * s.region["baseMVA"] for i in P_history]
	Q_history = [i * s.region["baseMVA"] for i in Q_history]
	
	#results = int2ext(ppc)
	
	"""
	fig = plt.figure()
	curfig = fig.add_subplot(1, 2, 1)
	curfig.plot(P_history, color = 'red', linewidth = 2.5, label = 'P')
	curfig = fig.add_subplot(1, 2, 2)
	curfig.plot(Q_history, color = 'blue', linewidth = 2.5, label = 'Q')
	curfig.set_yscale('log')
	curfig.legend(loc='upper right')
	"""
	
	#plt.show()
	
	"""
	
	print(all_gas)
	
	while( obj.txpool.inspect.pending.__dict__): pass
	end_time2 = time()
	final_balance = obj.eth.getBalance(obj.toChecksumAddress(client.eth_accounts()[0]))
	print("sasasa",ID, "    ",start_balance - final_balance , "Time ", end_time2 - start_time2)
	"""
	
	output.put((ID - 1, Result, ppc["bus"], ppc["gen"], ppc["branch"] , success, lamP2,ids, P_history, Q_history))
	
	"""
	plt.plot(all_gas, color = 'red', label = 'Estimated Gas Used')
	plt.hlines((start_balance - final_balance) ,xmin=0, xmax=nu, color='green', label = 'Actual Gas Used')
	plt.ylabel('Gas Used (wei)')
	plt.xlabel('Iteration')
	plt.title('Region %d' %(ID ))
	plt.legend()
	plt.show()
	"""
	
def getNextValues(ID,s,client,obj,num, count):
 		

 	dest = list(s.nbor.keys())
 	flag=0
 	
 	for k in dest:
 		j=0
 		temp_nonce=0
 		while(1):
 			j=j+1
 			if(j>count): break
 			trans={}
 			trans= obj.eth.getTransactionFromBlock(num, count-j).__dict__
 			#print(trans['from'])
 			if(trans['from'].lower() == client.eth_accounts()[0]):
 				hexjson = trans['input'][2:]
 				jsonstring = codecs.decode(hexjson, "hex").decode('utf-8')
 				jsonvalues = json.JSONDecoder().decode(jsonstring)
 				values_dict = json.loads(jsonvalues, cls=MyDecoder)
 				
 				if( (trans['nonce'] > temp_nonce) and values_dict['tID'] ==k ) :
 					
 					temp_nonce = trans['nonce']
 					print(trans['nonce'])
 					flag=1
 					my_count = count-j
    
 			#block = client.eth_getBlockByNumber(num, num-1)
 		if(flag):
 			last_nonce_dict = obj.eth.getTransactionFromBlock(num, my_count).__dict__
 			hexjson = last_nonce_dict['input'][2:]
 			jsonstring = codecs.decode(hexjson, "hex").decode('utf-8')
 			jsonvalues = json.JSONDecoder().decode(jsonstring)
 			values_dict = json.loads(jsonvalues, cls=MyDecoder)
 			print(values_dict)
 			tlidx = s.nbor[k].tlidx['int']
 			s.var['zmd'][tlidx] = numpy.asarray(values_dict['fields']['AVmd'])
 			s.var['zms'][tlidx] = numpy.asarray(values_dict['fields']['AVms'])
 			s.var['zad'][tlidx] = numpy.asarray(values_dict['fields']['AVad'])
 			s.var['zas'][tlidx] = numpy.asarray(values_dict['fields']['AVas'])
 		
 			#s.var['AVmd'][tlidx] = numpy.asarray(values_dict['fields']['AVmd'])
 			
 			length_key = len(values_dict['fields'])
 			if (length_key >10):
 				s.var['Pg']= numpy.asarray(values_dict['fields']['Pg'])
 				s.var['Qg'] = numpy.asarray(values_dict['fields']['Qg'])
 				s.var['Va'] = numpy.asarray(values_dict['fields']['Va'])
 				s.var['Vm'] = numpy.asarray(values_dict['fields']['Vm'])
 				#s.var['AVms'][tlidx] = numpy.asarray(values_dict['fields']['AVms'])
 				#s.var['AVmd'][tlidx] = numpy.asarray(values_dict['fields']['AVmd'])
 				#s.var['AVad'][tlidx] = numpy.asarray(values_dict['fields']['AVad'])
 				#s.var['AVas'][tlidx] = numpy.asarray(values_dict['fields']['AVas'])
 				s.var['rho'] =  values_dict['fields']['rho'] 
 				
 				
 				#s.pb['x0'] = r_[s.var['Va'], s.var['Vm'], s.var['Pg'], s.var['Qg']]
 			#s.var['AVms'][tlidx] = numpy.asarray(values_dict['fields']['AVms'])
 			#s.var['AVad'][tlidx] = numpy.asarray(values_dict['fields']['AVad'])
 			#s.var['AVas'][tlidx] = numpy.asarray(values_dict['fields']['AVas'])
 		
 				s.var['ymd'][tlidx] = numpy.asarray(values_dict['fields']['ymd'])
 				s.var['yms'][tlidx] = numpy.asarray(values_dict['fields']['yms'])
 				s.var['yad'][tlidx] = numpy.asarray(values_dict['fields']['yad'])
 				s.var['yas'][tlidx] = numpy.asarray(values_dict['fields']['yas'])
 		
 			#s.var['rho'] = values_dict['fields']['rho']

