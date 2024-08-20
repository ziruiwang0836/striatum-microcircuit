import numbers
from operator import le
import random
import numpy as np
from numpy.lib.npyio import _save_dispatcher
import pd5 as model_s# import parameters file
from netpyne import sim,analysis  # import netpyne init module
import matplotlib as mpl
import matplotlib.pyplot as plt

for i in range(11,20):
    model_s.model(i)



'''Ach_list =[0.9,0.7,0.5,0.3,0.1]
for Ach in Ach_list:
    SI_list=[]   
    for i in range(6,16):

        SI = model_s.model(i,Ach)
        SI_list.append(SI)

    filename = 'SI_pd_Ach21_'+str(Ach*100)+'.txt'
    np.savetxt(filename, SI_list, delimiter='\t')'''
'''cell_data_list = []
for i in range(220):  
    key = f'cell_{i}'
    if key in v_trace_data:
        cell_data_list.append(v_trace_data[key])

cell_data_array = np.array(cell_data_list)

np.savetxt(filename_v, cell_data_array, delimiter='\t')
'''