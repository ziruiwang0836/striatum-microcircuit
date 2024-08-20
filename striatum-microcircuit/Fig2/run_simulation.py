from NeuroTools.parameters import ParameterSet
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterSpace
import numpy as np
import striatum as st
import params_control as sham
import params_pd as pd
import params_pd_glu as pdr
import params_control_Glu as glu
import nest
import matplotlib.pyplot as plt
import pandas as pds
import os
import pickle

output_dir = '/Users/peirui/code/striatum-microcircuit/striatum-microcircuit/Fig2/output/'
os.makedirs(output_dir, exist_ok=True)
#run simualation
def run_simulation(sim_num,parms):
    for i in range(1,sim_num+1):
        id = i
        for experiment in parms.get_parameters().iter_inner():
            print(experiment)
            nest.ResetKernel()
            st.run(experiment,id,output_dir)

run_simulation(2,sham)
run_simulation(2,pd)
run_simulation(2,pdr)
run_simulation(2,glu)






