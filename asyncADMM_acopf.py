
import sys
import time
import os
import numpy as np 
from numpy import flatnonzero as find
import matplotlib.pyplot as plt

from scipy.io import loadmat
from os.path import join
from multiprocessing import Process, Pipe, Queue
import queue
from numpy import flatnonzero as find

from pypower.loadcase import loadcase
from pypower.case14 import case14
from pypower.case39 import case39
from pypower.case30 import case30
from pypower.ext2int import ext2int
from pypower.makeYbus import makeYbus
from pypower.int2ext import int2ext
from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF
from pypower.idx_gen import GEN_STATUS, GEN_BUS, PMAX, PMIN, QMAX, QMIN, PG, QG
from pypower.idx_cost import COST
from pypower.idx_brch import F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, \
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX
from pypower.isload import isload

from findTieline import findTieline
from runWorker import runWorker
from opf_admm_model import opf_admm_model

from pypower.api import printpf, runopf, ppoption

from numpy import array, exp, pi,logical_or
from numpy import conj, zeros, c_, shape, ix_








os.system("taskset -p 0xff %d" % os.getpid())

#---------------- basic system configuration  ---------------------------------
#caseFile = join('/pypower', 'case14')
# load a test power system


ppc = loadcase('case39.py')
# convert to internal numbering, remove out-of-service stuff 
na = 3# number of regions

"""
allgens = ppc["gen"][:,GEN_BUS]
print(allgens)

for i in range(len(allgens)):
    busto = find( ppc["bus"][:,BUS_I].astype(int) == allgens[i].astype(int))
    ppc["bus"][busto,BUS_AREA] = 2
    print(ppc["bus"][busto,:])
"""

nb = ppc["bus"].shape[0]

for i in range(nb):
    ppc["bus"][i,PD] = 1 * ppc["bus"][i,PD]


#ppc["bus"][:,PD] =  ppc["bus"][:,PD]



ppc_cen = ppc
ppopt = ppoption(PF_ALG=2)
r = runopf(ppc_cen)

gen = r["gen"]



objCent = r["f"]
#gencost = ppc_cen["gencost"]

bus=r["bus"]
nb = r["bus"].shape[0]
gen=r["gen"]
gencost=ppc_cen["gencost"]
objValue1=0
objValue2=0

"""
for i in range(nb):
            
            
            g  = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) &
                        ~isload(gen))
            ld = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) &
                        isload(gen))
            if any(g + 1):
                
                objValue1 = objValue1 + (bus[i, LAM_P]) * (sum(gen[g, PG]))
            

            if logical_or(bus[i, PD], bus[i, QD]) | any(ld + 1):
                if any(ld + 1):
                    
                    objValue2 = objValue2 + (bus[i, LAM_P]) * (bus[i, PD] - sum(gen[ld, PG]))
                    
objValue21=0

for i in range(len(gen)):
        if (gen[i,PG] > 0 ):
            objValue1 = objValue1 - gencost[i,COST+1]*gen[i,PG]
        else:
            objValue21 = objValue21 + gencost[i,COST+1]* (-gen[i,PG])
            
objValue2 = objValue21 - objValue2
objCent = objValue1 + objValue2 

"""



print(objCent)    


true_prices = bus[:,LAM_P]
bus_true_prices = bus[:,BUS_I]


ppc = ext2int(ppc)
baseMVA, bus, gen, branch, gencost = \
		ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"], ppc["gencost"]
slack = find(bus[:, BUS_TYPE] == 3) # an array of indices
gen_active = find(gen[:, GEN_STATUS] == 1)
genBus = gen[gen_active, GEN_BUS]

## convert to p.u.

    
bus[:,[PD, QD]] /= baseMVA

branch2=branch


gen[:,[PMAX, PMIN, QMAX, QMIN]] /= baseMVA

gencost[:, COST] *= baseMVA ** 2
gencost[:, COST + 1] *= baseMVA




## problem dimensions
nb = bus.shape[0]          # number of buses
nl = branch.shape[0]       # number of branches
ng = len(gen_active)   # number of piece-wise linear costs

## build admittance matrices
Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)

##---------- partition the system ---------------------------------------------

# the input format for partition should be a nb*1 matlab matrix with elements denoting area number
#partitionName = 'OP14_2region.mat' 
# partitionName = 'OP118_w05_40region_sizediff.mat'
# partitionName = 'OP118_rd_8region_scale100000.mat'
#partitionDir = '/Users/junyaoguo/Dropbox/ADMM/Partition'
#partitionFile = join(partitionDir,partitionName)
#partition = loadmat(partitionName)
# partitionnumber = 0
#op = partition['OP14'][0]
#op[10:] = 3
# op[:] = 1

# op = partition['C'][partitionnumber]
#bus[:, BUS_AREA] = op
op = bus[:, BUS_AREA]

tieline, ntl = findTieline(bus, branch)

#pdb.set_trace()


##---------- create all the communication pipes --------------------------------
edges = np.vstack({tuple(sorted(row)) for row in tieline[:, 2:4]}) if tieline.any() else np.array([])

"""
pipes = {}
for edge in edges.tolist():
	fend, tend =  Pipe()
	if edge[0] not in pipes:
		pipes[edge[0]] = {}
	pipes[edge[0]][edge[1]] = fend
	if edge[1] not in pipes:
		pipes[edge[1]] = {}
	pipes[edge[1]][edge[0]] = tend


"""

##----subproblem configuration including local opf and communication pipes-----
problem = []
output = Queue()
for i in range(na):
	s = opf_admm_model()
	s.config(i + 1, op, bus, gen, gencost,branch, Ybus, Yf,Yt, genBus, tieline, None, na, baseMVA)
	s.var_init()
	problem.append(s)


##----- run each worker in parallel ---------
procs = []
for i in range(na):
	procs += [Process(target = runWorker, args = (i + 1, problem[i], output))]

start_time = time.time()
start_clock = time.clock()
#print("ndfsifsdfidsjhfifhdi")
#sys.exit()
for proc in procs:
	proc.start()


# TIMEOUT = 70
# bool_list = [True] * na
# start = time.time()
# while time.time() - start <= TIMEOUT:
# 	for i in range(na):
# 		bool_list[i] = procs[i].is_alive()

# 	#print(bool_list)

# 	if np.any(bool_list):  
# 		time.sleep(.1)  
# 	else:
# 		break
# else:
# 	print("Timed out, killing all processes")
# 	for p in procs:
# 		p.terminate()

liveprocs = list(procs)	
results = []

while liveprocs:
	try:
		while 1:
			results.append(output.get(False))
			
	except queue.Empty:
		pass

	time.sleep(0.5)
	if not output.empty():
		continue

	liveprocs = [p for p in liveprocs if p.is_alive()]

for proc in procs:
	proc.join()

# results = []
# for proc in procs:
# 	while proc.is_alive():
# 		proc.join(timeout = 30)
# 		while True:
# 			try: 
# 				results.append(output.get())
# 			except Empty:
# 				break

## ------------get results ---------------------------
ttime = time.time()
tclock = time.clock()
totaltime = ttime - start_time
clocktime = tclock - start_clock

results = sorted(results, key=lambda x: x[0])

objTotal = 0





for k in range(na):
    
	if k != results[k][0]:
		print('Error: Result of worker %d not returned!' % (k+1,))
		break
	objTotal += results[k][1]['objValue']

gap = abs(objTotal - objCent) / objCent * 100
print('The convergence time is %f' % (totaltime,))
print('The convergence clock time is %f' % (clocktime,))
print('The objective function value is %f' % (objTotal,))
print('The gap in objective function is %f %%' % (gap,))

## ------------ plots of convergence -----------------
fig = plt.figure()
ppc={}


diff_P = []
Gen_ID= []


for k in range(na):
	if k != results[k][0]:
		print('Error: Result of worker %d not returned!' % (k+1,))
		break
	pgap = results[k][1]['primalgap']
	dgap = results[k][1]['dualgap']
	bus = results[k][2]
	gen = results[k][3]
	branch = results[k][4]
	success = results[k][5]
	prc =results[k][6]
	ids =results[k][7]
	P_history= results[k][8]
	Q_history= results[k][9]
	
	
	P_corr=0
	Q_corr=0

	
	for l in range(len(gen)):
		
		reg_gen = find(r['gen'][:,GEN_BUS].astype(int) ==  gen[l,GEN_BUS].astype(int))
		
		diff_P.append( 100*abs(r['gen'][reg_gen,PG] - gen[l,PG]) / r['gen'][reg_gen,PG])
		Gen_ID.append(gen[l,GEN_BUS].astype(int))
		
		P_corr = P_corr + r['gen'][reg_gen,PG]
		Q_corr = Q_corr + r['gen'][reg_gen,QG]
	 
	
	#print(P_corr)
	
	

	
	if(k==0):
		ppc["bus"]=bus     
		ppc["branch"]=branch
		ppc["gen"]=gen
		ppc["success"]=success
		prices=prc
		bus_ids=ids
		
	else:
		ppc["bus"]= np.unique(np.concatenate((ppc["bus"],bus),0),  axis=0)
		ppc["branch"]= np.unique(np.concatenate((ppc["branch"],branch),0),  axis=0)
		ppc["gen"]= np.unique(np.concatenate((ppc["gen"],gen),0),  axis=0)
		prices = np.concatenate([prices,prc])
		bus_ids = np.concatenate([bus_ids,ids])
		if(ppc["success"]==0):pass
		else:ppc["success"]=success    
	
	
	
	curfig = fig.add_subplot(1, 3, k + 1)
	if(k==2): curfig.plot(pgap, color = 'red', linewidth = 2.5, label = 'primal residual')
	else: curfig.plot(pgap, color = 'red', linewidth = 2.5)
	
	if(k==2): curfig.plot(dgap, color = 'blue', linewidth = 2.5, label = 'dual residual')
	else: curfig.plot(dgap, color = 'blue', linewidth = 2.5)
	
	curfig.set_yscale('log')
	curfig.legend(loc='upper right')
	curfig.xaxis.set_tick_params(labelsize=12)
	curfig.yaxis.set_tick_params(labelsize=12)
	
	
	
	
	
	fig2 = plt.figure()
	
	curfig2 = fig2.add_subplot(1, 1, 1)
	curfig2.plot(P_history, color = 'red', linewidth = 2.5, label = 'P')
	curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')
	curfig2.set_title('Region %d' %(k+1), fontsize=25 )
	curfig2.set_xlabel('iterations', fontsize=25)
	curfig2.set_ylabel('Total P (Generation in MW)  ', fontsize=25)
	
	curfig2.xaxis.set_tick_params(labelsize=25)
	curfig2.yaxis.set_tick_params(labelsize=25)

	curfig2.legend()
	
	
	
	"""
	curfig2 = fig2.add_subplot(1, 1, 1)
	curfig2.plot(Q_history, color = 'blue', linewidth = 2.5, label = 'Q')
	curfig2.hlines(Q_corr, xmin=0, xmax=len(Q_history), color='green',linewidth = 2.5,label='Optimal Q')
	curfig2.set_title('Region %d' %(k+1), fontsize=18 )
	curfig2.set_xlabel('iterations', fontsize=17)
	curfig2.set_ylabel('Total Q (Generation in MVAR)  ', fontsize=17)
	curfig2.xaxis.set_tick_params(labelsize=17)
	curfig2.yaxis.set_tick_params(labelsize=17)

	curfig2.legend()

	"""




fig2 = plt.figure()	
curfig2 = fig2.add_subplot(1, 1, 1)

curfig2.plot(Gen_ID, diff_P, 'bo' , label = 'Error in Power Generation')

#print(diff_P)

#curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')

curfig2.set_xlabel('Generator Bus_ID')
curfig2.set_ylabel('Error (%) ')
curfig2.legend()


price_by_bus=prices

Va =ppc["bus"][:,VA]  * pi / 180
Vm=ppc["bus"][:,VM]
V = Vm * exp(1j * Va)

if shape(branch2)[1] < MU_ANGMAX + 1:
		branch2 = c_[branch2, zeros((nl, MU_ANGMAX + 1 - shape(branch2)[1]))]
	

Ybus, Yf, Yt = makeYbus(baseMVA, ppc["bus"], branch2)
#print(Yf)
## compute branch flows
Sf = V[ branch2[:, F_BUS].astype(int) ] * conj(Yf * V)  ## cplx pwr at "from" bus, p["u"].
St = V[ branch2[:, T_BUS].astype(int) ] * conj(Yt * V)  ## cplx pwr at "to" bus, p["u"].
	
	
branch2[:, PF] = Sf.real * baseMVA
branch2[:, QF] = Sf.imag * baseMVA
branch2[:, PT] = St.real * baseMVA
branch2[:, QT] = St.imag * baseMVA


branch2[:,F_BUS]=branch2[:,F_BUS]+1
branch2[:,T_BUS]=branch2[:,T_BUS]+1
ppc["branch"]=branch2
ppc["baseMVA"]=baseMVA
ppc["et"]=clocktime



ppc=ext2int(ppc)
solution = int2ext(ppc)
#printpf(solution)

print('The convergence time is %f' % (totaltime,))
print('The convergence clock time is %f' % (clocktime,))
print('The objective function value is %f' % (objTotal,))
print('The gap in objective function is %f %%' % (gap,))

#print(price_by_bus)
#print(bus_ids)


fig3=plt.figure()
curfig3 = fig3.add_subplot(1, 1, 1)
curfig3.plot( bus_ids, price_by_bus, 'ro', label='ADMM Prices')
curfig3.plot( bus_true_prices, true_prices, 'g+', label='True Prices')
curfig3.legend(fontsize=12)
curfig3.set_xlabel('Bus_ID', fontsize=14)
curfig3.set_ylabel('Price ($/MWh) ', fontsize=14)
curfig3.xaxis.set_tick_params(labelsize=14)
curfig3.yaxis.set_tick_params(labelsize=14)
    
#print(bus_true_prices, true_prices)
#diff = [(100*abs(true_prices[i]-price_by_bus[i])/true_prices[i]) for i in range(min(len(true_prices), len(price_by_bus)))]

#print(diff)


plt.legend()
#plt.ylim(ymin=0)
plt.show()  

##--------some quick commands to check the gradient and hessian of subproblem--------

# f, df = s1.admmopf_costfcn(result['x'])
# h, g, dh, dg = s1.admmopf_consfcn(result['x'])
# H = s1.admmopf_hessfcn(result['x'])

