from NeuroTools.parameters import ParameterSet
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterSpace

def get_parameters():
    p = ParameterSpace({})
    p.outpath = '.'

    # Parameters for running
    p.timestep = 0.1
    p.min_delay = 0.1
    p.max_delay = 5.1
    p.runtime = 4000

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

    p.prefix = 'Control'

    return p
