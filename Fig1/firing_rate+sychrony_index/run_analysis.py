from NeuroTools.parameters import ParameterSet
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterSpace
import numpy as np
import straitum as st
import params_pd as pd
import params_pd_GluInh as pdr
import params_control_Glu as glu
import nest
import matplotlib.pyplot as plt
import pandas as pds

p = ParameterSpace({})
p.outpath = '.'

# Parameters for running
p.timestep = 0.1
p.min_delay = 0.1
p.max_delay = 5.1
p.runtime = 1200

# Parameters for neuronal features
#D1 and D2
p.gl1 = 12.5
p.th1 = -45
p.el1 = -80
p.E_tau = 0.3
p.I_tau = 2
p.cm1 = 200
p.E_ex1 = 0
p.E_in1 = -64
p.reset1 = -80
p.d1_size = 2000
p.d2_size = 2000
p.fsi_size = 80
#FSI
p.gl_f = 25
p.th_f = -54
p.elf = -80
p.cmf = 500
p.E_exf = 0
p.E_inf = -76
p.reset1 = -80

#Parameters for connection
#d12d1
p.d12d1_p = 0.26
p.d12d1_w = -0.5
p.d12d1_d = 2
#d12d2
p.d12d2_p = 0.07
p.d12d2_w = -1
p.d12d2_d = 2
#d22d2
p.d22d2_p = 0.36
p.d22d2_w = -1
p.d22d2_d = 2
#d22d1
p.d22d1_p = 0.27
p.d22d1_w = -1.2
p.d22d1_d = 2
#fsi2d1
p.fsi2d1_p = 0.54
p.fsi2d1_w = -2.5
p.fsi2d1_d = 1
#fsi2d2 
p.fsi2d2_p = 0.36
p.fsi2d2_w = -2.5
p.fsi2d2_d = 1
# background input
p.d1_bk = 2500
p.d2_bk = 2500
p.fsi_bk = 2500
p.cor2d1_w = 3.55
p.cor2d2_w = 3.50
p.cor2fsi_w = 3.8

#Parameters for MIP
p.N = 1000
p.r = 400     # MIP rate
p.c = 0.02 # with-pool correlation in MIP
p.q = 0
p.edge = 0 


msn_sham = []
msn_pd=[]
msn_pdr=[]
msn_glu=[]
sham_d1_firing=[]
sham_d2_firing=[]
pd_d1_firing=[]
pd_d2_firing=[]
pdr_d1_firing=[]
pdr_d2_firing=[]
glu_d1_firing=[]
glu_d2_firing=[]
for i in range(1,11):
    r_msn_syn=[]
    for experiment in p.iter_inner():
        nest.ResetKernel()
        print(experiment)
        model = st.Striatum()
        r_msn,r_d1,r_d2 = model.run(p)
        r_msn_syn.append([r_msn,r_d1,r_d2])
       

    r_msn_syn_pd = pd.get_pd_results()
    r_msn_syn_pd_recover = pdr.get_pdr_results()
    r_msn_syn_sham_glu = glu.get_sham_glu_results()
    msn_sham.append(r_msn_syn[0][0])
    msn_pd.append(r_msn_syn_pd[0][0])
    msn_pdr.append(r_msn_syn_pd_recover[0][0])
    msn_glu.append(r_msn_syn_sham_glu[0][0])

    sham_d1_firing.append(r_msn_syn[0][1])
    sham_d2_firing.append(r_msn_syn[0][2])
    pd_d1_firing.append(r_msn_syn_pd[0][1])
    pd_d2_firing.append(r_msn_syn_pd[0][2])
    pdr_d1_firing.append(r_msn_syn_pd_recover[0][1])
    pdr_d2_firing.append(r_msn_syn_pd_recover[0][2])
    glu_d1_firing.append(r_msn_syn_sham_glu[0][1])
    glu_d2_firing.append(r_msn_syn_sham_glu[0][2])


msn_sham_ = pds.DataFrame(msn_sham)
msn_pd_ = pds.DataFrame(msn_pd)
msn_pdr_ = pds.DataFrame(msn_pdr)
msn_glu_=pds.DataFrame(msn_glu)
msn_sham_f = list(msn_sham_.mean(axis=0))
msn_pd_f = list(msn_pd_.mean(axis=0))
msn_pdr_f = list(msn_pdr_.mean(axis=0))
msn_glu_f = list(msn_glu_.mean(axis=0))

sham_d1 = np.mean(sham_d1_firing)
sham_d2 = np.mean(sham_d2_firing)
pd_d1 = np.mean(pd_d1_firing)
pd_d2 = np.mean(pd_d2_firing)
pdr_d1 = np.mean(pdr_d1_firing)
pdr_d2 = np.mean(pdr_d2_firing)
glu_d1 = np.mean(glu_d1_firing)
glu_d2 = np.mean(glu_d2_firing)

# saving output to excel sync_dync
msn_result = pds.DataFrame(msn_sham_f)
msn_result['PD'] = pds.DataFrame(msn_pd_f)
msn_result['recover'] = pds.DataFrame(msn_pdr_f)
msn_result['Glu'] = pds.DataFrame(msn_glu_f)
msn_result.columns=['Sham','PD','PD Glu block','Sham Glu activate']
#msn_result.to_excel('result_syn.xlsx')

# saving output to excel firing rate
firing_result = pds.DataFrame(sham_d1_firing)
firing_result['d2'] = pds.DataFrame(sham_d2_firing)
firing_result['d1_pd'] = pds.DataFrame(pd_d1_firing)
firing_result['d2_pd'] = pds.DataFrame(pd_d2_firing)
firing_result['d1_pdr'] = pds.DataFrame(pdr_d1_firing)
firing_result['d2_pdr'] = pds.DataFrame(pdr_d2_firing)
firing_result['d1_glu'] = pds.DataFrame(glu_d1_firing)
firing_result['d2_glu'] = pds.DataFrame(glu_d2_firing)
firing_result.columns=['d1','d2','d1_pd','d2_pd','d1_pdr','d2_pdr','d1_glu','d2_glu']
#firing_result.to_excel('firing_result_q01.xlsx')


'''def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')
window_size=2
smoothed_sham = moving_average(msn_sham_f, window_size)
smoothed_pd = moving_average(msn_pd_f, window_size)
smoothed_pdr = moving_average(msn_pdr_f, window_size)
plt.plot(np.arange(0,875,25),smoothed_sham[:-4],linewidth=2,color='black',label='Sham')
plt.plot(np.arange(0,875,25),smoothed_pd[:-4],linewidth=2,color='red',label='PD')
plt.plot(np.arange(0,875,25),smoothed_pdr[:-4],linewidth=2,color='blue',label='Glu block')
plt.title('MSNs',fontweight='bold')
plt.xticks(fontweight='bold')
plt.ylabel('synchrony index',fontweight='bold')
plt.yticks(fontweight='bold')
plt.xlabel('Time(ms)',fontweight='bold')
legend = plt.legend()
for text in legend.get_texts():
    text.set_fontweight('bold')'''
#plt.show()
'''
r_msn_syn_pd = pd.get_pd_results()
plt.plot([0.3,0.6],r_d1_l,color='green',alpha=0.5,linestyle='dashed',marker='o',label='D1_control')
plt.plot([0.3,0.6],r_d2_l,color='green',alpha=0.5,linestyle='solid',marker='o',label='D2_control')
plt.plot([0.3,0.6],r_d1_pd,color='grey',alpha=0.5,linestyle='dashed',marker='o',label='D1_PD')
plt.plot([0.3,0.6],r_d2_pd,color='grey',alpha=0.5,linestyle='solid',marker='o',label='D2_PD')
plt.title('Mean firing rate under different cortical input correlation')
plt.ylabel('mean firing rate')
plt.xlabel('input correlation')
plt.legend()
plt.show()
#plt.savefig('./mean_firing_rate_normal.jpg')
'''





