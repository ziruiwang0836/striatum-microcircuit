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

spike_senders = dict()
spike_times = dict()
rates = dict()

rates["d1"] =[]
rates["d2"] =[]
rates["fsi"]=[]

spike_senders["d1"]=[]
spike_senders["d2"]=[]
spike_senders["fsi"]=[]

spike_times["d1"]=[]
spike_times["d2"]=[]
spike_times["fsi"]=[]

def run(prefix,p):   
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

    #%%spike decetor
    spike_det_d1 = nest.Create("spike_recorder")
    spike_det_d2 = nest.Create("spike_recorder")
    spike_det_fsi = nest.Create("spike_recorder")
    nest.Connect(d1_pop, spike_det_d1)
    nest.Connect(d2_pop, spike_det_d2)
    nest.Connect(fsi_pop, spike_det_fsi)

    #simulate
    nest.print_time = True
    nest.Simulate(p.runtime)

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
    spike_senders['d1'].append(evs_d1)
    spike_senders['d2'].append(evs_d2)
    spike_senders['fsi'].append(evs_fsi)

    spike_times['d1'].append(ts_d1)
    spike_times['d2'].append(ts_d2)
    spike_times['fsi'].append(ts_fsi)

    # Calculate the mean firing rate of D1,D2and FSI populations
    secs = float(p.runtime)/1000.
    rateD1 = (len(ts_d1)/secs)/(float(p.d1_size))
    rateD2 = (len(ts_d2)/secs)/(float(p.d2_size))
    rateFSI = (len(ts_fsi)/secs)/(float(p.fsi_size))

    rates['d1'].append(rateD1)
    rates['d2'].append(rateD2)
    rates['fsi'].append(rateFSI)
        


    # save data
    name = 'Senders_'+prefix+'.pickle'
    pickle.dump(spike_senders,open(name,"wb"))

    name = 'Times_'+prefix+'.pickle' 
    pickle.dump(spike_times,open(name,"wb"))

    name = 'Rates_'+prefix+'.pickle' 
    pickle.dump(rates,open(name,"wb"))














