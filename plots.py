import matplotlib.pyplot as plt
from numpy import array, exp, pi,logical_or



Gen_ID = [ 30 ,31, 32,33,34,35,36,37,38,39]

#flat
diff_P1 = [array([5.19510892e-06]), array([0.83152992]), array([0.60835272]), array([2.27387991]), array([1.99638759e-06]), array([5.16829286e-05]), array([1.32079795e-06]), array([0.66482115]), array([2.223209e-06]), array([0.03514308])]


#warm1
diff_P2 = [array([7.24556462]), array([2.43666109]), array([3.2053196]), array([1.28875631]), array([1.99675401e-06]), array([5.16975803e-05]), array([1.32038623e-06]), array([2.48155608]), array([2.22274955e-06]), array([1.44719645])]

#perf1
diff_P22 = [array([5.19580942e-06]), array([0.49627107]), array([0.25236618]), array([0.86520207]), array([1.99949741e-06]), array([5.1675828e-05]), array([1.32098063e-06]), array([0.31899904]), array([2.22346766e-06]), array([0.03928142])]


#perf and 5% up
diff_P23=[array([4.79779646e-07]), array([7.13204236e-07]), array([1.72355344]), array([12.37780483]), array([8.1666866e-07]), array([5.49450704e-07]), array([3.86845013e-07]), array([6.2388767e-07]), array([4.3727153e-07]), array([4.78454441])]


#warm2
diff_P3 = [array([5.76142116]), array([2.2875509]), array([2.88541334]), array([0.52491458]), array([1.99787943e-06]), array([5.16920342e-05]), array([1.31907956e-06]), array([2.96336725]), array([2.22081852e-06]), array([0.93475173])]


# 5% up
diff_P4=[array([4.75832636e-07]), array([0.19855336]), array([3.52713024]), array([5.12176774]), array([8.17139835e-07]), array([5.5266567e-07]), array([3.88889462e-07]), array([6.27647345e-07]), array([4.39661521e-07]), array([3.07397978])]

#5% down from up
diff_P5=[array([15.19239799]), array([1.92663531]), array([3.41631633]), array([5.41653615]), array([3.17583397e-05]), array([3.13425492]), array([1.3960928e-05]), array([3.1806891]), array([4.28203703e-05]), array([3.1092555])]


#5% down 
diff_P6=[array([9.90850883]), array([1.91997686]), array([3.30843927]), array([0.91604823]), array([3.1752871e-05]), array([4.48131786]), array([1.39617433e-05]), array([4.56153228]), array([4.28231283e-05]), array([1.02708165])]


#rho = 0.8
diff_P7 = [array([5.84115172]), array([0.94473987]), array([1.82003269]), array([0.37800902]), array([3.17492193e-05]), array([1.47506583]), array([1.39615091e-05]), array([1.50136551]), array([4.2819766e-05]), array([1.8195046])]



#rho 0.8 from warm 
diff_P8 = [array([19.65939445]), array([1.0968132]), array([7.90727308]), array([7.52746343]), array([2.00000151e-06]), array([5.17102827e-05]), array([1.3217206e-06]), array([3.86201854]), array([2.22510279e-06]), array([2.71066351])]

#1.01
diff_P9 =[array([10.93636861]), array([2.11737317]), array([8.92122901]), array([2.75705073]), array([9.94628516e-07]), array([5.21258138e-06]), array([6.90165516e-07]), array([1.92171093]), array([1.08565484e-06]), array([2.38478681])]

#1.03
diff_P10 = [array([18.29385479]), array([1.12604551]), array([9.1398954]), array([7.68959596]), array([9.94889008e-07]), array([5.21222597e-06]), array([6.8865855e-07]), array([1.92171096]), array([1.08400869e-06]), array([0.42664159])]

#0,95   0.8 rho
diff_P11 = [array([13.98534656]), array([1.35394714]), array([2.21123971]), array([5.08696273]), array([3.17583406e-05]), array([2.88810846]), array([1.39607496e-05]), array([2.93047164]), array([4.28197229e-05]), array([4.28319269])]




"""
fig2 = plt.figure()
curfig2 = fig2.add_subplot(1, 1, 1)

curfig2.plot(Gen_ID, diff_P1, 'bo' , label = '(a) Flat Start')
curfig2.plot(Gen_ID, diff_P4, 'r*' , label = '(b) Warm Start and 5% Demand Increase')
curfig2.plot(Gen_ID, diff_P6, 'g+' , label = '(c) Warm Start and 5% Demand Reduction')

#curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')

curfig2.set_xlabel('Generator Bus_ID')
curfig2.set_ylabel('Error (%) ')
curfig2.set_title('39-Bus Case')
curfig2.legend()
plt.show() 

"""




fig2 = plt.figure()
curfig2 = fig2.add_subplot(1, 1, 1)


curfig2.plot(Gen_ID, diff_P1, 'b*' , label = 'Flat Start ' )
curfig2.plot(Gen_ID, diff_P8, 'r+' , label = 'Warm Start and ρ = 0.8*ρ ')


#curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')

curfig2.set_xlabel('Generator Bus_ID' , fontsize=14)
curfig2.set_ylabel('Error (%) ', fontsize=14)
curfig2.set_title('39-Bus Case', fontsize=15)
curfig2.xaxis.set_tick_params(labelsize=14)
curfig2.yaxis.set_tick_params(labelsize=14)


curfig2.legend(fontsize=12)
plt.show() 


Gen_ID2 = [ 1,2,22,27,23,13]



res1 = [array([2.50029068]), array([3.58238753]), array([11.34914709]), array([4.15889821]), array([7.82740928]), array([21.1513774])]
res2=[array([3.90446644]), array([3.53687278]), array([10.16714495]), array([7.03487471]), array([1.41885114]), array([17.09005106])]

res3 = [array([4.65764592]), array([4.09171752]), array([7.58481535]), array([3.85825967]), array([1.54360587]), array([16.9808505])]
"""
fig2 = plt.figure()
curfig2 = fig2.add_subplot(1, 1, 1)

curfig2.plot(Gen_ID2, res1, 'bo' , label = 'Flat Start and Primal Residual = 10^-3')
curfig2.plot(Gen_ID2, res2, 'r*' , label = 'Warm Start and Primal Residual = 10^-4')
curfig2.plot(Gen_ID2, res3, 'g+' , label = 'Warm Start and Primal Residual = 10^-5')

#curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')

curfig2.set_xlabel('Generator Bus_ID')
curfig2.set_ylabel('Error (%) ')
curfig2.set_title('30-Bus Case')
curfig2.legend()
plt.show() 
"""

"""
diffs = [sum(res1) , sum(res2), sum(res3) ]

print(diffs)

diffs2 = [50.56951019, 43.15226108, 38.71689483]
 
y_pos = ['Residual = 10^-3', ' Residual = 10^-4', ' Residual = 10^-5']

pm, pc, pn = plt.bar( range(len(y_pos)), diffs2 , tick_label=y_pos)
pm.set_facecolor('b')
pc.set_facecolor('r')
pn.set_facecolor('g')

plt.ylabel('Total Error(%)')
plt.title('30-Bus Case')


plt.rcdefaults()


plt.show()
"""


#flat 10^-3
diff_P1 = [array([5.19510892e-06]), array([0.83152992]), array([0.60835272]), array([2.27387991]), array([1.99638759e-06]), array([5.16829286e-05]), array([1.32079795e-06]), array([0.66482115]), array([2.223209e-06]), array([0.03514308])]


#flat 10^-4
r2 = [array([5.19997034e-06]), array([0.36960968]), array([0.27139403]), array([2.67968229]), array([1.99065082e-06]), array([5.16599914e-05]), array([1.31778971e-06]), array([1.17391097]), array([2.21842737e-06]), array([0.04672429])]

r3 = [array([40.42334875]), array([0.82066504]), array([23.81847635]), array([22.16460622]), array([2.00041225e-06]), array([5.17120582e-05]), array([1.32239287e-06]), array([3.86201854]), array([2.22602244e-06]), array([7.45946607])]



fig2 = plt.figure()
curfig2 = fig2.add_subplot(1, 1, 1)

curfig2.plot(Gen_ID, diff_P1, 'bo' , label = 'ρ = 10^4 ')
curfig2.plot(Gen_ID, r2, 'r*' , label = 'ρ = 10^5 ')
curfig2.plot(Gen_ID, r3, 'g+' , label = 'ρ = 10^6 ')

#curfig2.hlines(P_corr, xmin=0, xmax=len(P_history), color='green',linewidth = 2.5, label='Optimal P')

curfig2.set_xlabel('Generator Bus_ID' , fontsize=11)
curfig2.set_ylabel('Error (%) ', fontsize=11)
curfig2.set_title('39-Bus Case',fontsize=16)
curfig2.xaxis.set_tick_params(labelsize=11)
curfig2.yaxis.set_tick_params(labelsize=11)
curfig2.legend()
plt.show() 
