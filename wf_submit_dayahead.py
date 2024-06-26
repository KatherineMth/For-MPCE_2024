# -*- coding: utf-8 -*-
"""
Created on 2024.06.24
@author: Tianhui Meng
"""
# This is the code of how to get dayahead schedule curve
# Use Gurobi solver

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import openpyxl  

# TCL paramenter
from TCLs_paramenter import ptcl_b1_max
from TCLs_paramenter import ptcl_b2_max
from TCLs_paramenter import ptcl_con_max 
from TCLs_paramenter import ptcl_con_min
from TCLs_paramenter import ptcl_acp_max 
from TCLs_paramenter import ptcl_acp_min
# MG paramenter
pmg_con_pos_max = 5*np.ones(127)
pmg_con_neg_max = 5*np.ones(127)

#%% main func
file1 = openpyxl.load_workbook('File path/wind_data.xlsx') # wind forecast
winddata_all = file1['Sheet1']
windp = winddata_all['B3:DX3']
windp_value = [[cell.value for cell in rows] for rows in windp]
windp_ori = np.array(windp_value)
windp = windp_ori/12.5
windp1 = np.transpose(windp) 
pwf = windp1.flatten() 

windscenario_all = file1['Sheet2'] # scenario
winds = windscenario_all['B1:DX6']
winds_value = [[cell.value for cell in rows] for rows in winds]
winds_value1= np.transpose(winds_value)
pwa = np.array(winds_value1)

pr = windscenario_all['B7:G7'] # Probability of each scenario
pr_value = [[cell.value for cell in rows] for rows in pr]
pr_value = np.transpose(pr_value)
pr = np.array(pr_value)

ptcl_b1_max = ptcl_b1_max[1:128]
ptcl_b2_max = ptcl_b2_max[1:128]
ptcl_con_max = ptcl_con_max[1:128]
ptcl_con_min = ptcl_con_min[1:128]
ptcl_acp_max = ptcl_acp_max[1:128]
ptcl_acp_min = ptcl_acp_min[1:128]

#%%
model=gp.Model("wf_submit_dayahead")

Ts = 127
ws = 6
pwmax = 100
adr = 0.1 # allowable deviation rate

#%%
#  valley-0：00-07：45; flat-08：00-09：45 12：00-13：45 19：00-23：45; peak-10：00-11：45 14：00-18：45
# feed-in tariff
cgt_ori = 0.4*np.ones(96)*1000 # 
cgt_ori[0:32] = 0.2*1000
cgt_ori[40:48] = cgt_ori[56:76] = 0.8*1000

cact_ori = 1.25*cgt_ori # WF to TCL original

cw2mgt_ori = 0.8*cgt_ori # WF to MG original
cmg2wt_ori = 2*cgt_ori # MG to WF original

cgt = np.concatenate((cgt_ori[81:96],cgt_ori[:],cgt_ori[0:16]),axis = 0)
cact = np.concatenate((cact_ori[81:96],cact_ori[:],cact_ori[0:16]),axis = 0)
cact1 = cact # WF to TCL seg_1
cact2 = 8.5/12.5*cact # WF to TCL seg_2

cw2mgt = np.concatenate((cw2mgt_ori[81:96],cw2mgt_ori[:],cw2mgt_ori[0:16]),axis = 0)
cmg2wt = np.concatenate((cmg2wt_ori[81:96],cmg2wt_ori[:],cmg2wt_ori[0:16]),axis = 0)
cw2mgt1 = cw2mgt # WF to MG seg_1
cw2mgt2 = cgt*0.7 # WF to MG seg_2
cw2mgt3 = cgt*0.5 # WF to MG seg_3

cmg2wt1 = cmg2wt*0.7 # MG to WF seg_1
cmg2wt2 = cmg2wt*0.8 # MG to WF seg_2
cmg2wt3 = cmg2wt*0.9 # MG to WF seg_3

#  Penalty coefficient for wind curtailment
cwc= 0.05*1000   

#  Penalty factors for positive/negative deviations
cdev_ins = 1.8*1000 
cdev_sup = 1.5*1000

#%%
pws = model.addVars(Ts, vtype=GRB.CONTINUOUS, lb=0, ub=pwmax, name="pws") # submit curve
pwg = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, ub=pwmax, name="pwg") # W-FJOS to main grid
pwc = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, ub=pwmax, name="pwc") # wind curtailment
pwac_inc = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac_inc") # TCLs power increase
pwac_dec = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac_dec") # TCLs power decrease
pwac = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac") # WF to TCLs
pwac1 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac1") # wind power to 1th segment
pwac2 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac2") # wind power to 2th segment
pwac3 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pwac3") # wind power to 3th segment

cwac = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="cwac") # Wind farm electricity sales revenue to TCLs
cmg = model.addVar(vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="cmg")
cwg = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="cwg") # electricity sales revenue of the W-FJOS
model.addConstr(cwg == 0.25*gp.quicksum(cgt[i]*pr[j,0]*pwg[i,j] for i in range(Ts) for j in range(ws)))   

# limits on the rate of change of the schedule curve (can be ignored)
for i in range(1,Ts):        
    model.addConstr(pws[i]-pws[i-1]<=5)
    model.addConstr(pws[i]-pws[i-1]>=-5)

# TCLs constraints
for i in range(Ts):
    for j in range(ws):
        model.addConstr(pwac[i,j] <= ptcl_acp_max[i])
        model.addConstr(pwac_inc[i,j] == pwac2[i,j]+pwac3[i,j])
        model.addConstr(pwac_dec[i,j] <= ptcl_b1_max[i]) 
        
for i in range(Ts):
    for j in range(ws):
        model.addConstr(pwac1[i,j] == ptcl_con_min[i]-pwac_dec[i,j])
        model.addConstr(pwac1[i,j] <= ptcl_con_min[i])
        model.addConstr(pwac2[i,j] <= ptcl_con_max[i]-ptcl_con_min[i])
        model.addConstr(pwac3[i,j] <= ptcl_acp_max[i]-ptcl_con_max[i])
        model.addConstr(pwac[i,j] == pwac1[i,j]+pwac2[i,j]+pwac3[i,j])        
        model.addConstr(pwac[i,j] <= 7.31) # Maximum power flow constraint        
model.addConstr(cwac == 0.25*gp.quicksum(pr[j,0]*(cact1[i]*pwac1[i,j]+cact2[i]*pwac2[i,j]) for i in range(Ts) for j in range(ws)))             

# MG constraints
pmg = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY,  name="pmg") # power between WF and MG
pw2mg = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pw2mg")
pw2mg1 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pw2mg1") # seg_1
pw2mg2 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pw2mg2") # seg_2
pw2mg3 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pw2mg3") # seg_3
pmg2w = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pmg2w")
pmg2w1 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pmg2w1") # seg_1
pmg2w2 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pmg2w2") # seg_2
pmg2w3 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0, name="pmg2w3") # seg_3

for i in range(Ts):
    for j in range(ws): 
        model.addConstr(pw2mg[i,j] == pw2mg1[i,j]+pw2mg2[i,j]+pw2mg3[i,j])
        model.addConstr(pw2mg1[i,j] <= pmg_con_pos_max[i]/3)
        model.addConstr(pw2mg2[i,j] <= pmg_con_pos_max[i]/3)
        model.addConstr(pw2mg3[i,j] <= pmg_con_pos_max[i]/3)        
        model.addConstr(pmg2w[i,j] == pmg2w1[i,j]+pmg2w2[i,j]+pmg2w3[i,j])        
        model.addConstr(pmg2w1[i,j] <= pmg_con_neg_max[i]/3)
        model.addConstr(pmg2w2[i,j] <= pmg_con_neg_max[i]/3)
        model.addConstr(pmg2w3[i,j] <= pmg_con_neg_max[i]/3)        
        model.addConstr(pmg2w[i,j]*pw2mg[i,j] == 0) 
        model.addConstr(pw2mg[i,j]-pmg2w[i,j] == pmg[i,j])     
        
model.addConstr(cmg == 0.25*gp.quicksum((pr[j,0]*(cw2mgt1[i]*pw2mg1[i,j]+cw2mgt2[i]*pw2mg2[i,j]+cw2mgt3[i]*pw2mg3[i,j])-pr[j,0]*(cmg2wt1[i]*pmg2w1[i,j]+cmg2wt2[i]*pmg2w2[i,j]+cmg2wt3[i]*pmg2w3[i,j]) for i in range(Ts) for j in range(ws))))

# power balance
for i in range(Ts):
    for j in range(ws):
        model.addConstr(pwa[i,j]-pwg[i,j] == pwc[i,j]+pwac[i,j]+pmg[i,j])

# wind curtailment penalty
cwcur = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="cwcur")
model.addConstr(cwcur == 0.25*gp.quicksum(cwc*pr[j,0]*pwc[i,j] for i in range(Ts) for j in range(ws)))     

# deviation penalty
pdev_sup =  model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pbias_sup") 
pdev_ins =  model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=0,  name="pbias_ins")
A1 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=-100, ub=100, name="A1") # Auxiliary variables
A2 = model.addVars(Ts, ws, vtype=GRB.CONTINUOUS, lb=-100, ub=100, name="A2") # Auxiliary variables
cdev =  model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="cdev")
for i in range(Ts):
    for j in range(ws):
        model.addConstr(A1[i,j] == (1-adr)*pws[i]-pwg[i,j])
        model.addConstr(A2[i,j] == pwg[i,j]-(1+adr)*pws[i])
        model.addConstr(pdev_ins[i,j] == gp.max_(A1[i,j], 0))
        model.addConstr(pdev_sup[i,j] == gp.max_(A2[i,j], 0)) 

model.addConstr(cdev == 0.25*gp.quicksum(pr[j,0]*(cdev_ins*pdev_ins[i,j]+cdev_sup*pdev_sup[i,j]) for i in range(Ts) for j in range(ws)))


#%% Solving the Model
ctotal = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="ctotal")
model.addConstr(ctotal == cwg+cwac+cmg-cwcur-cdev)

model.setObjective(ctotal, GRB.MAXIMIZE)
model.setParam('nonconvex', 2) 
model.Params.TimeLimit = 50

model.optimize()

#%% optimize result
print('Obj: %g' % model.ObjVal)
print(ctotal.x)
print('%s = %g' % (cwg.VarName, cwg.X))
print('%s = %g' % (cwac.VarName, cwac.X))
print('%s = %g' % (cmg.VarName, cmg.X))
print('%s = %g' % (cwcur.VarName, cwcur.X))

