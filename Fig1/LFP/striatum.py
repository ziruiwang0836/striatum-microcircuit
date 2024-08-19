#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 17:43:00 2022

@author: peirui
"""

import numpy as np
import pylab
import matplotlib.pyplot as plt
import nest  
import time                   
import pickle
import os

def run(p,id,output_dir):  
    #Building networks: neurons
    #Parameters for neuronal features
    d1_n_dict = {"g_L": p.gl1,
                    "V_th":p.th1,
                    "E_L":p.el1,
                    "tau_syn_ex":p.E_tau,
                    "tau_syn_in":p.I_tau,
                    "C_m":p.cm1,
                    "E_in":p.E_in1,
                    "E_ex":p.E_ex1,
                    "V_reset":p.reset1}
    d2_n_dict = {"g_L": p.gl1,
                    "V_th":p.th1,
                    "E_L":p.el1,
                    "tau_syn_ex":p.E_tau,
                    "tau_syn_in":p.I_tau,
                    "C_m":p.cm1,
                    "E_in":p.E_in1,
                    "E_ex":p.E_ex1,
                    "V_reset":p.reset1}
    fsi_n_dict = {"g_L": p.gl_f,
                    "V_th":p.th_f,
                    "E_L":p.elf,
                    "tau_syn_ex":p.E_tau,
                    "tau_syn_in":p.I_tau,
                    "C_m":p.cmf,
                "E_in":p.E_inf,
                    "E_ex":p.E_exf,
                    "V_reset":p.reset1}
    #%%
    #membrane potential Initializations 
    d1_pop = nest.Create("iaf_cond_alpha", p.d1_size, params=d1_n_dict)
    d2_pop = nest.Create("iaf_cond_alpha", p.d2_size, params=d2_n_dict)
    fsi_pop = nest.Create("iaf_cond_alpha", p.fsi_size, params=fsi_n_dict)
    d1_pop.set({"V_m":nest.random.uniform(-80.0, -55.0)})
    d2_pop.set({"V_m":nest.random.uniform(-80.0, -55.0)})
    fsi_pop.set({"V_m":nest.random.uniform(-80.0, -55.0)})
    #%%
    #Building networks: connections
    #IF_cond_alpha
    ##fb
    #d12d1
    fb_d12d1_connector = {'rule': 'pairwise_bernoulli', 'p': p.d12d1_p}
    fb_d12d1_syn_dict = {
                        'weight': p.d12d1_w, 'delay': p.d12d1_d}
    nest.Connect(d1_pop, d1_pop, conn_spec=fb_d12d1_connector, syn_spec=fb_d12d1_syn_dict)
    #d12d2
    fb_d12d2_connector = {'rule': 'pairwise_bernoulli', 'p': p.d12d2_p}
    fb_d12d2_syn_dict = {
                        'weight': p.d12d2_w, 'delay': p.d12d2_d}
    nest.Connect(d1_pop, d2_pop, conn_spec=fb_d12d2_connector, syn_spec=fb_d12d2_syn_dict)
    #d22d2
    fb_d22d2_connector = {'rule': 'pairwise_bernoulli', 'p': p.d22d2_p}
    fb_d22d2_syn_dict = {
                        'weight': p.d22d2_w, 'delay': p.d22d2_d}
    nest.Connect(d2_pop, d2_pop, conn_spec=fb_d22d2_connector, syn_spec=fb_d22d2_syn_dict)
    #d22d1
    fb_d22d1_connector = {'rule': 'pairwise_bernoulli', 'p': p.d22d1_p}
    fb_d22d1_syn_dict = {
                        'weight': p.d22d1_w, 'delay': p.d22d1_d}
    nest.Connect(d2_pop, d1_pop, conn_spec=fb_d22d1_connector, syn_spec=fb_d22d1_syn_dict)
    ##ff
    #fsi2d1
    ff_fsi2d1_connector = {'rule': 'pairwise_bernoulli', 'p': p.fsi2d1_p}
    fb_fsi2d1_syn_dict = {
                        'weight': p.fsi2d1_w, 'delay': p.fsi2d1_d}
    nest.Connect(fsi_pop, d1_pop, conn_spec=ff_fsi2d1_connector, syn_spec=fb_fsi2d1_syn_dict)
    #fsi2d2
    ff_fsi2d2_connector = {'rule': 'pairwise_bernoulli', 'p': p.fsi2d2_p}
    fb_fsi2d2_syn_dict = {
                        'weight': p.fsi2d2_w, 'delay': p.fsi2d2_d}
    nest.Connect(fsi_pop, d2_pop, conn_spec=ff_fsi2d2_connector, syn_spec=fb_fsi2d2_syn_dict)
    #%%
    #Define background input
    cort_d1 = nest.Create('poisson_generator',1,{'rate':p.d1_bk})
    cort_d2 = nest.Create('poisson_generator',1,{'rate':p.d2_bk})
    cort_fsi = nest.Create('poisson_generator',1,{'rate':p.fsi_bk})
    cort_d1_connector = 'all_to_all'
    cort_d1_syn_dict= {'weight':p.cor2d1_w}                   
    cort_d2_connector = 'all_to_all'
    cort_d2_syn_dict= {'weight':p.cor2d2_w}
    cort_fsi_connector = 'all_to_all'
    cort_fsi_syn_dict =  {'synapse_model': 'static_synapse',
                        'weight': p.cor2fsi_w}
    nest.Connect(cort_d1,d1_pop, conn_spec=cort_d1_connector, syn_spec=cort_d1_syn_dict)
    nest.Connect(cort_d2,d2_pop, conn_spec=cort_d2_connector, syn_spec=cort_d2_syn_dict)
    nest.Connect(cort_fsi,fsi_pop, conn_spec=cort_fsi_connector, syn_spec=cort_fsi_syn_dict)
    #%%
    # Define MIP
    #c with-pool correlation in MIP
    #q between-cell correlation of MIP input     
    def mip(num):
        spk_times_list=[]
        cort_e = nest.Create('spike_generator',num)
        proc_c = np.around(pylab.rand(int(round(p.q*p.r*p.runtime/(1000.*p.N*p.c))))*p.runtime+3.1)
        for i in range(num):
            proc_ind_n = np.around(pylab.rand(int(round((1-p.q)*p.r*p.runtime/(1000.*p.N*p.c))))*p.runtime+3.1) # this value is added to avoid negative spike times
            spk_times = np.empty(0)
            for j in range(int(round((1-p.q)*p.r*p.runtime/(1000.*p.N*p.c)))):  
                pulse_size = pylab.binomial(p.N,p.c)
                pulse = proc_ind_n[j] + p.edge*(pylab.rand(pulse_size)-0.5)*2
                spk_times = np.append(spk_times,pulse)
            for j in range(int(round(p.q*p.r*p.runtime/(1000.*p.N*p.c)))):
                pulse_size = pylab.binomial(p.N,p.c)
                pulse = proc_c[j] + p.edge*(pylab.rand(pulse_size)-0.5)*2
                spk_times = np.append(spk_times,pulse)
            spk_times.sort()
            cort_e[i].spike_times = spk_times
        return cort_e

    if p.c ==0:
        #msn
        cort_e = nest.Create('poisson_generator',p.d1_size,{'rate':(1-p.q)*p.r})
        cort_c = nest.Create('poisson_generator',p.d1_size,{'rate':p.q*p.r})
        par = nest.Create('parrot_neuron',p.d1_size)
        nest.Connect(cort_e,d1_pop,'one_to_one')
        nest.Connect(cort_e,d2_pop,'one_to_one')           
        nest.Connect(cort_c,par,'one_to_one')
        nest.Connect(par,d1_pop,'one_to_one')
        nest.Connect(par,d2_pop,'one_to_one')
        #fsi
        cort_e_f = nest.Create('poisson_generator',p.fsi_size,{'rate':(1-p.q)*p.r})
        cort_c_f = nest.Create('poisson_generator',p.fsi_size,{'rate':p.q*p.r})
        par_f = nest.Create('parrot_neuron',p.fsi_size)
        nest.Connect(cort_e_f,fsi_pop,'one_to_one')
        nest.Connect(cort_c_f,par_f,'one_to_one')
        nest.Connect(par_f,fsi_pop,'one_to_one')
    else:
        cort_m_mip_msn = mip(2000)
        cort_f_mip = mip(80)
        nest.Connect(cort_m_mip_msn,d1_pop,'one_to_one')
        nest.Connect(cort_m_mip_msn,d2_pop,'one_to_one')
        nest.Connect(cort_f_mip,fsi_pop,'one_to_one')

    #multimeter
    mul_d1 = nest.Create('multimeter', {'record_from': ['V_m', 'g_ex','g_in']})
    mul_d2 = nest.Create('multimeter', {'record_from': ['V_m', 'g_ex','g_in']})   
    #nest.SetStatus(mul_d1, {'interval': 0.1})
    #nest.SetStatus(mul_d2, {'interval': 0.1})
    nest.Connect(mul_d1, d1_pop)
    nest.Connect(mul_d2, d2_pop)
    
    #spike decetor
    spike_det_d1 = nest.Create("spike_recorder")
    spike_det_d2 = nest.Create("spike_recorder")
    spike_det_fsi = nest.Create("spike_recorder")
    nest.Connect(d1_pop, spike_det_d1)
    nest.Connect(d2_pop, spike_det_d2)
    nest.Connect(fsi_pop, spike_det_fsi)

    #simulate
    nest.print_time = True
    nest.Simulate(p.runtime)

    #spk data
    dSD_d1 = nest.GetStatus(spike_det_d1)[0]
    dSD_d2 = nest.GetStatus(spike_det_d2)[0]
    #synaptic current data used to calcualte LFP
    data_d1 = nest.GetStatus(mul_d1,'events')[0]
    data_d2 = nest.GetStatus(mul_d2,'events')[0]

    
    # Define file names
    name_spk_d1 = f'spk_{p.prefix}_detectors_d1_{id}.pickle'
    name_spk_d2 = f'spk_{p.prefix}_detectors_d2_{id}.pickle'
    name_lfp_d1 = f'lfp_{p.prefix}_detectors_d1_{id}.pickle'
    name_lfp_d2 = f'lfp_{p.prefix}_detectors_d2_{id}.pickle'

    # Save data to files
    try:
        with open(os.path.join(output_dir, name_spk_d1), "wb") as file:
            pickle.dump(dSD_d1, file)

        with open(os.path.join(output_dir, name_spk_d2), "wb") as file:
            pickle.dump(dSD_d2, file)

        with open(os.path.join(output_dir, name_lfp_d1), "wb") as file:
            pickle.dump(data_d1, file)

        with open(os.path.join(output_dir, name_lfp_d2), "wb") as file:
            pickle.dump(data_d2, file)

    except Exception as e:
        print(f"An error occurred while saving files: {e}")




    '''# Record mean firing rate
    spike_senders_d1=[]
    spike_senders_d2=[]
    spike_senders_fsi=[]

    spike_times_d1=[]
    spike_times_d2=[]
    spike_times_fsi=[]

    # Read from spike detectors
    #d1
    dSD_d1 = nest.GetStatus(spike_det_d1)[0]
    evs_d1 = dSD_d1['events']['senders']
    ts_d1 = dSD_d1['events'] ['times']
    #d2
    dSD_d2 = nest.GetStatus(spike_det_d2)[0]
    evs_d2 = dSD_d2['events']['senders']
    ts_d2 = dSD_d2['events'] ['times']
    #fsi
    dSD_fsi = nest.GetStatus(spike_det_fsi)[0]
    evs_fsi = dSD_fsi['events']['senders']
    ts_fsi = dSD_fsi['events'] ['times']

    # Store spike times and IDs in list
    spike_senders_d1.append(evs_d1)
    spike_senders_d2.append(evs_d2)
    spike_senders_fsi.append(evs_fsi)

    spike_times_d1.append(ts_d1)
    spike_times_d2.append(ts_d2)
    spike_times_fsi.append(ts_fsi)


    # Calculate the mean firing rate of D1,D2and FSI populations
    secs = float(p.runtime)/1000.
    rateD1 = (len(ts_d1)/secs)/(float(p.d1_size))
    rateD2 = (len(ts_d2)/secs)/(float(p.d2_size))
    rateFSI = (len(ts_fsi)/secs)/(float(p.fsi_size))


    # calculate synchrony index
    # 5 ms as one time bin
    # refer to Yim et al. 2011
    def calculate_synchrony_index(ts_pop):
        count_lis=[]
        for i in range(0,int(p.runtime+1),5):
            count=0
            for t in ts_pop:
                if i<t<=t+5:
                    count=count+1
            count_lis.append(count)
        mean_count = np.mean(count_lis)
        var_count = np.var(count_lis)
        synchrony_index = var_count/mean_count
        return synchrony_index
    syn_inx_d1 = calculate_synchrony_index(ts_d1)
    syn_inx_d2 = calculate_synchrony_index(ts_d2)
    syn_inx_msn = (syn_inx_d2+syn_inx_d1)/2

    return syn_inx_msn,rateD1,rateD2'''















