# -*- coding: utf-8 -*-
"""
Created on 2024.06.24
@author: Tianhui Meng
"""
# This is the code for the TCL agent to calculate the available power


import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

Ts = 129 # 20：00-23：45 0：00-23：45 0：00-04：00, total 129 point
N = 3000 
Tout_ori = [33,32,32,31.2,31.2,31.2,32,33,35,36,35,35,33,32,31.2,31.2,31.2] # Outside temperature
Tset0 = 24
R=1.75
delta_Tup = 3 # Tmax-Tset0, Maximum temperature gap that the users can accept
delta_Tdn = 3 # Tset0-Tmin
k1 = 0.03 # paramenters from reference[26]
k2 = 0.06
l1 = -0.4
l2 = -0.3
s = k2/k1
l = l2-k2*l1/k1

# Temperature interpolation
x1 = np.linspace(0,128,17) # 2+13+2=17
fx1 = interp1d(x1, Tout_ori, kind='linear')
x2 = np.linspace(0,128,129) 
Tout = fx1(x2) # Call the interpolation function to get the temperature every 15 min
plt.plot(x2,Tout)
plt.show()
Tmax = Tset0+delta_Tup # Allowable upper temperature
Tmin = Tset0-delta_Tdn # Allowable lower temperature

zg_ability = 1/10 # 
ptcl_set0 = np.zeros(Ts) # total initial power
ptcl_max = np.zeros(Ts) # max(total power)
ptcl_min = np.zeros(Ts) # min(total power)
ptcl_b1_max = np.zeros(Ts) # power that can be supplemented by DN, corresponding to func(A.3) in paper
ptcl_b2_max = np.zeros(Ts) # power that can be reduced by DN
ptcl_con_max = np.zeros(Ts) # corresponding to (Ltcl,2) in paper
ptcl_con_min = np.zeros(Ts) # corresponding to (Ltcl,1) in paper
ptcl_acp_max = np.zeros(Ts) # corresponding to (Ltcl,max) in paper
ptcl_acp_min = np.zeros(Ts) # corresponding to (Ptcl,con,min) in paper

for i in range(Ts):
    ptcl_set0[i] = N*((Tout[i]-Tset0)/(s*R)-l/s)/1000
    ptcl_max[i] = N*((Tout[i]-Tmin)/(s*R)-l/s)/1000
    ptcl_min[i] = N*((Tout[i]-Tmax)/(s*R)-l/s)/1000
    ptcl_b1_max[i] = 1/(zg_ability+1)*zg_ability*(ptcl_min[i]-2)
    ptcl_b2_max[i] = zg_ability*(ptcl_min[i]-2-ptcl_b1_max[i])
    ptcl_con_max[i] = ptcl_set0[i]-ptcl_min[i]+ptcl_b1_max[i]+2 
    ptcl_con_min[i] = 2+ptcl_b1_max[i] 
    ptcl_acp_max[i] = ptcl_max[i]-ptcl_min[i]+ptcl_b2_max[i]+2+ptcl_b1_max[i] 
    ptcl_acp_min[i] = 2 
   

