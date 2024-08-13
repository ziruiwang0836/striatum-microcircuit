from NeuroTools.parameters import ParameterSet
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterSpace
import numpy as np
import straitum as st
import nest

def get_sham_glu_results():

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
    p.cor2d1_w = 3.55*1.1
    p.cor2d2_w = 3.5*1.1
    p.cor2fsi_w = 3.8*1.1

    #Parameters for MIP
    p.N = 1000
    p.r = 400     # MIP rate
    p.c = 0.02#ParameterRange(np.arange(0,0.06,0.01)) # with-pool correlation in MIP
    p.q = 0
    p.edge = 0 

    # Record mean firing rate
    r_msn_syn_sham_glu=[]
    nest.ResetKernel()
    model = st.Striatum()
    r_msn_sham_glu,r_d1_glu,r_d2_glu = model.run(p)
    r_msn_syn_sham_glu.append([r_msn_sham_glu,r_d1_glu,r_d2_glu])
    return r_msn_syn_sham_glu

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






