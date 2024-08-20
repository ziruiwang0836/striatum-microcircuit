from netpyne import specs,analysis
import matplotlib.pyplot as plt
import math
import netpyne
from netpyne.analysis import lfp
import random
import numpy as np
from neuron import test_rxd
from netpyne import sim,analysis  # import netpyne init module


netParams = specs.NetParams()   # object of class NetParams to store the network parameters
simConfig = specs.SimConfig()   # object of class SimConfig to store the simulation configuration

###############################################################################
#
# M1 6-LAYER ynorm-BASED MODEL
#
###############################################################################

###############################################################################
# SIMULATION CONFIGURATION
###############################################################################

# Simulation parameters
simConfig.duration =1500 # Duration of the simulation, in ms
simConfig.dt = 0.05 # Internal integration timestep to use
simConfig.seeds = {'conn': 1, 'stim': 1, 'loc': 1} # Seeds for randomizers (connectivity, input stimulation and cell locations)
simConfig.createNEURONObj = 1  # create HOC objects when instantiating network
simConfig.createPyStruct = 1  # create Python structure (simulator-independent) when instantiating network
simConfig.verbose = 0 # Whether to write diagnostic information on events 
pd=0
simConfig.printRunTime = 0.1
#simConfig.parallel = True
#simConfig.nhosts = 4  # Set the number of MPI processes (equal to the number of cores)

# Recording 
simConfig.recordCells = ['all']  # list of cells to record from 
simConfig.recordLFP= ['all']  # list of cells to record from 
simConfig.recordTraces = {'V':{'sec':'soma','loc':0.5,'var':'v'}} 
simConfig.recordStim = True  # record spikes of cell stims
simConfig.recordStep = 0.1 # Step size in ms to save data (eg. V traces, LFP, etc)

# Saving
simConfig.filename = 'Health'  # Set file output name
simConfig.saveFileStep = 1000 # step size in ms to save data to disk
simConfig.savePickle = False # save to pickle file
simConfig.saveJson =False # save to json file
simConfig.saveMat =  False  # save to mat file
simConfig.saveTxt =   False # save to txt file
simConfig.saveDpk = False # save to .dpk pickled file
simConfig.saveHDF5 = False # save to HDF5 file 


# Analysis and plotting 
simConfig.analysis['plotRaster'] = {'include': [ 'I23','I23L','E23','I5','I5L','E5A','E5B',
'E5P','StrD1','StrD2','StrFSI','STN','GPi','GPe'],'orderInverse': True,'timeRange': [500,1000], 'figSize': (12,10), 'legend':False,
'lw': 0.3,'markerSize':12, 'marker': '.', 'dpi': 500, 'xlabel':[],'ylabel':[],'saveFig': 'm1bgt.png'} # Whether or not to plot a raster

#simConfig.analysis['plotRaster'] = {'include': [ 'StrD1','StrD2'],'orderInverse': True,'timeRange': [200,600], 'figSize': (12,10), 'legend':False,
#'lw': 0.3,'markerSize':12, 'marker': '.', 'dpi': 500, 'xlabel':[],'ylabel':[],'saveFig': 'm1bgt.png'} # Whether or not to plot a raster

#simConfig.analysis['plotRaster'] = {'include':  ['I23','I23L','E23','I5','I5L','E5A','E5B','E5P','StrD1','StrD2','StrFSI'],'orderInverse': True,
#'timeRange': [100,1000], 'figSize': (12,10), 'lw': 0.3,'markerSize':12, 'marker': '.', 'dpi': 500, 'saveFig': 'raster.png','syncLines': False} # Whether or not to plot a raster

#simConfig.analysis['plot2Dnet'] = {'include':['GPi','Th'],'showConns': True}
#simConfig.analysis['plotConn'] = {'saveFig': True}                                                 # plot connectivity matrix
#simConfig.recordLFP = [[60,215, 60],[60,590,60],[60,1020,60],[60,1220,60]]
simConfig.recordLFP = [[60,215,60],[60,590,60],[60,805,60],[60,905,60],[60,1050,60],[60,1150,60],[60,1220,60],[60,1350,60]]
simConfig.analysis['plotLFP'] = {'includeAxon': True, 'figSize': (6,10), 'timeRange': [300,1300], 'saveFig': False,'showFig':True} 
#simConfig.analysis['plotTraces'] = {'include': [('StrD1',20),('StrD2',20),('StrFSI',10)],'timeRange': [100,500],'rerun':False,'saveFig':True}
#simConfig.analysis['plot2Dfiring'] = {'saveFig': True, 'showFig': True}
# NETWORK PARAMETERS 
###############################################################################


# General network parameters
netParams.scale = 1 # Scale factor for number of cells
netParams.sizeX = 120*math.sqrt(netParams.scale) # x-dimension (horizontal length) size in um
netParams.sizeY = 1370 # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = 120*math.sqrt(netParams.scale) # z-dimension (horizontal depth) size in um

# 
netParams.scaleConnWeight = 0.2 # Connection weight scale factor
netParams.scaleConnWeightNetStims = 0.7# Connection weight scale factor for NetStims
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 100.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)

netParams.popParams['I23'] = { 'cellType': 'PV',  'yRange': [100, 310], 'density': 10e3} #  L2/3 PV (FS)
netParams.popParams['I23L'] ={'cellType': 'SOM', 'yRange':[100, 310], 'density': 10e3} #  L2/3 SOM (LTS)
netParams.popParams['E23'] = { 'cellType': 'IT',  'yRange': [120,310], 'density': 80e3} #  L2/3 IT
netParams.popParams['I5'] =  {'cellType': 'PV',  'yRange': [310,770], 'density': 10e3} #  L5 PV (FS)
netParams.popParams['I5L'] = { 'cellType': 'SOM', 'yRange': [310,770], 'density': 10e3} #  L5 SOM (LTS)
netParams.popParams['E5A'] = {'cellType': 'IT',  'yRange': [410,520], 'density': 80e3} #  L5A IT
netParams.popParams['E5B'] = {'cellType': 'IT',  'yRange': [520, 770], 'density': 40e3} #  L5B IT
netParams.popParams['E5P'] = { 'cellType': 'PT',  'yRange': [520, 770], 'density': 40e3} #  L5B PT

# for key in netParams.popParams.items:
netParams.popParams['StrD1'] = {'cellModel': 'StrD1', 'cellType': 'StrD1', 'numCells':200,'yRange': [770, 870]}
netParams.popParams['StrD2'] = {'cellModel': 'StrD2', 'cellType': 'StrD2', 'numCells': 200,'yRange': [870, 970]}
netParams.popParams['StrFSI'] = {'cellModel': 'StrFSI', 'cellType': 'StrFSI', 'numCells': 20,'yRange': [970, 1000]}
netParams.popParams['GPi'] = {'cellModel': 'GPi', 'cellType': 'GPi', 'numCells': 200,'yRange': [1000, 1100]}
netParams.popParams['GPe'] = {'cellModel': 'GPe', 'cellType': 'GPe',  'numCells': 200,'yRange': [1100, 1200]}
netParams.popParams['STN'] = {'cellModel': 'STN', 'cellType': 'STN', 'numCells': 200,'yRange': [1200, 1300]}
netParams.popParams['Th'] = {'cellModel': 'Th', 'cellType': 'Th', 'numCells': 200,'yRange': [1300, 1400]}


## IT cell params
RScell = {'conds': {'cellType': 'IT'}, 'secs': {}}
RScell['secs']['soma'] = {'geom': {}, 'mechs': {}}  
RScell['secs']['soma']['geom'] = {'diam': 61.4, 'L': 61.4, 'cm': 1,'Ra':150}
RScell['secs']['soma']['mechs']['cor_RS_IL'] = {}
RScell['secs']['soma']['mechs']['cor_RS_IK'] = {}
RScell['secs']['soma']['mechs']['cor_RS_INA'] = {}
RScell['secs']['soma']['mechs']['cor_RS_IM'] = {}
RScell['secs']['soma']['vinit']= -71.9
netParams.cellParams['IT'] = RScell # add dict to list of cell properties

## PT cell params
IBcell = {'conds': {'cellType': 'PT'}, 'secs': {}}
IBcell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
IBcell['secs']['soma']['geom'] = {'diam':96, 'L': 96,'Ra':150}
IBcell['secs']['soma']['mechs']['cor_IB_IK'] = {}
IBcell['secs']['soma']['mechs']['cor_IB_IL'] = {}
IBcell['secs']['soma']['mechs']['cor_IB_INA'] = {}
IBcell['secs']['soma']['mechs']['cor_IB_IM'] = {}
IBcell['secs']['soma']['mechs']['cor_IB_ICA'] = {}
IBcell['secs']['soma']['vinit']= -71.4
netParams.cellParams['PT'] = IBcell  # add dict to list of cell properties


## SOM cell params
LTScell = {'conds': {'cellType': 'SOM'}, 'secs': {}}
LTScell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
LTScell['secs']['soma']['geom'] = {'diam':89.2, 'L': 89.2,'Ra':150}
LTScell['secs']['soma']['mechs']['cor_LTS_IK'] = {}
LTScell['secs']['soma']['mechs']['cor_LTS_IL'] = {}
LTScell['secs']['soma']['mechs']['cor_LTS_INA'] = {}
LTScell['secs']['soma']['mechs']['cor_LTS_IM'] = {}
LTScell['secs']['soma']['mechs']['cor_LTS_ICA'] = {}
LTScell['secs']['soma']['vinit']= -54
netParams.cellParams['SOM'] = LTScell

## PV cell params
FScell = {'conds': {'cellType': 'PV'}, 'secs': {}}
FScell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
FScell['secs']['soma']['geom'] = {'diam': 56.9, 'L': 56.9, 'cm': 1,'Ra':150}
FScell['secs']['soma']['mechs']['cor_FS_Ik'] = {}
FScell['secs']['soma']['mechs']['cor_FS_IL'] = {}
FScell['secs']['soma']['mechs']['cor_FS_INA'] = {}
FScell['secs']['soma']['mechs']['cor_FS_IM'] = {}
FScell['secs']['soma']['vinit']= -71.4
netParams.cellParams['PV'] = FScell

#D1 and D2 neurons
cellRule = {'conds': {'cellModel': 'StrD1', 'cellType': 'StrD1'}, 'secs': {}}
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}
cellRule['secs']['soma']['geom'] = {'diam': 12.6157,
                                    'L': 12.6157,
                                    'Ra': 1,
                                    'nseg': 1}
cellRule['secs']['soma']['mechs']['Str'] = {'gmbar': (2.6e-3)}
cellRule['secs']['soma']['vinit'] = random.gauss(-63.8, 5)
netParams.cellParams['StrD1'] = cellRule

cellRule = {'conds': {'cellModel': 'StrD2', 'cellType': 'StrD2'}, 'secs': {}}
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}
cellRule['secs']['soma']['geom'] = {'diam': 12.6157,
                                    'L': 12.6157,
                                    'Ra': 1,
                                    'nseg': 1}
cellRule['secs']['soma']['mechs']['Str'] = {'gmbar': (2.6e-3)}
cellRule['secs']['soma']['vinit'] = random.gauss(-63.8, 5)
netParams.cellParams['StrD2'] = cellRule

#FSI neurons
#refer to Mccarthy et al. 2011
#FSI=MSN without m current + gap junction
#refer to cortical fast spiking interneurons, adding corresponding ion channels
#fsi threshold is higher than msn, tuned as -11
FSIcell = {'conds': {'cellModel': 'StrFSI', 'cellType': 'StrFSI'}, 'secs': {}}
FSIcell['secs']['soma'] = {'geom': {}, 'mechs': {}}
FSIcell['secs']['soma']['geom'] = {'diam': 12.6157,
                                    'L': 12.6157,
                                    'Ra': 1,
                                    'nseg': 1}
FSIcell['secs']['soma']['mechs']['cor_FS_Ik'] = {}
FSIcell['secs']['soma']['mechs']['cor_FS_IL'] = {}
FSIcell['secs']['soma']['mechs']['cor_FS_INA'] = {}
FSIcell['secs']['soma']['mechs']['cor_FS_IM'] = {}
FSIcell['secs']['soma']['vinit'] = random.gauss(-63.8, 5)
netParams.cellParams['StrFSI'] = FSIcell 

# GPe neruons
GPcell = {'secs': {}}
GPcell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
GPcell['secs']['soma']['geom'] = {'diam': 12.6157, 'L':12.6157, 'cm': 1,'Ra':100}
GPcell['secs']['soma']['mechs']['I_L_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_K_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_Na_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_T_gpe'] = {}
GPcell['secs']['soma']['mechs']['I_Ca_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_AHP_GPe']= {}
GPcell['secs']['soma']['vinit']= random.gauss(-62, 5)
netParams.cellParams['GPi'] =GPcell

#GPi neurons
GPcell = {'secs': {}}
GPcell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
GPcell['secs']['soma']['geom'] = {'diam': 12.6157, 'L':12.6157, 'cm': 1,'Ra':100}
GPcell['secs']['soma']['mechs']['I_L_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_K_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_Na_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_T_gpe'] = {}
GPcell['secs']['soma']['mechs']['I_Ca_GPe'] = {}
GPcell['secs']['soma']['mechs']['I_AHP_GPe']= {}
GPcell['secs']['soma']['vinit']=random.gauss(-62, 5)
netParams.cellParams['GPe'] =GPcell


# STN neurons
PYRcell = {'secs': {}}
PYRcell['secs']['soma'] = {'geom': {}, 'mechs': {}}  # soma params dict
PYRcell['secs']['soma']['geom'] = {'diam': 12.6157, 'L':12.6157, 'cm': 0.1,'Ra':100}
PYRcell['secs']['soma']['mechs']['I_L_STN'] = {}
PYRcell['secs']['soma']['mechs']['I_Na_STN'] = {}
PYRcell['secs']['soma']['mechs']['I_T_stn'] = {}
PYRcell['secs']['soma']['mechs']['I_Ca_STN'] = {}
PYRcell['secs']['soma']['mechs']['I_K_STN'] = {}
PYRcell['secs']['soma']['mechs']['I_AHP_STN'] = {}
PYRcell['secs']['soma']['vinit']= random.gauss(-62, 5)
netParams.cellParams['STN'] = PYRcell

# TH neurons
cellRule = {'conds': {'cellModel': 'Th', 'cellType': 'Th'}, 'secs': {}}
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}
cellRule['secs']['soma']['geom'] ={'diam': 5.642,
                                    'L': 5.642,
                                    'Ra': 1,
                                    'nseg': 1}
cellRule['secs']['soma']['mechs']['thalamus'] = {}
cellRule['secs']['soma']['vinit'] = random.gauss(-62, 5)
netParams.cellParams['Th'] = cellRule

# Synaptic mechanism parameters
netParams.synMechParams['AMPA'] = {'mod': 'AMPA_S'}  # excitatory synaptic mechanism
netParams.synMechParams['GABAB'] = {'mod': 'GABAa_S'}  # inhibitory synaptic mechanism

# refer to Ding et al. 2008, cortical-striartal projection NMDA/AMPA=2.75/1,
# thamalo-striatal projectioni NNDA/AMPA=2.04/1
# refer to Moyer et al. 2007, DA acting on MSN D1R is excitory(*130%),
# DA acting on MSN D2R is inhibitory(*80%)refer to Lindahl et al. 2016 convert DA modulation 
#into synaptic current I(.mod file)-DA level [0,1];0 stands for high DA, 1 stands for low DA
# refer to Evans et al. 2012 specify Glu2B as the main NMDA channel
# refer to Wall et al. 2013 D1 receives 20% fewer input than D2
# refer to Damodaran et al. 2015 gapjunction between FSI, prob=0.3
glu_f = 1
netParams.synMechParams['Glu_AMPA_d1'] = {'mod':'AMPAk','gAmax':3420*glu_f,'alpha_DA':0.8}
netParams.synMechParams['Glu_NMDA_d1'] = {'mod':'NMDAk','gNmax':9400*glu_f, 'tauon':2.25,'tauoff':150,
'mg':3.57,'alpha_DA':1,'beta_nmda':0.7}
netParams.synMechParams['Glu_AMPA_d2'] = {'mod':'AMPAk','gAmax':3420*glu_f,'alpha_DA':1,'beta_ampa':-0.8}
netParams.synMechParams['Glu_NMDA_d2'] = {'mod':'NMDAk','gNmax':9400*glu_f, 'tauon':2.25,'tauoff':150,
'mg':3.57,'alpha_DA':0.8}
netParams.synMechParams['esyn'] = {'mod': 'ElectSyn'}

# GPe
netParams.synMechParams['Insge,ampa'] = {'mod': 'Exp2Syn',
                                                'tau1': 0.4,
                                                'tau2': 2.5,
                                                'e': 0}  # stn -> gpe
netParams.synMechParams['Insge,nmda'] = {'mod': 'Exp2Syn',
                                                'tau1': 2,
                                                'tau2': 67,
                                                'e': 0}  # stn -> gpe
netParams.synMechParams['Igege'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': -85}  # gpe -< gpe
netParams.synMechParams['Istrgpe'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': -85}  # D2 -> gpe
netParams.synMechParams['Ithco'] = {'mod': 'Exp2Syn',
                                                 'tau1': 5,
                                                 'tau2': 5,
                                                 'e': 0}  # th->rs

netParams.synMechParams['Igith'] = {'mod': 'Exp2Syn',
                                                    'tau1': 5,
                                                    'tau2': 5,
                                                    'e': -85}  # gpi -<th
netParams.synMechParams['Igestr'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': -85}  # gpe -< str
# GPi
netParams.synMechParams['Igegi'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': -85}  # gpe -< gp
netParams.synMechParams['Isngi'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': 0}  # stn -> gpi
netParams.synMechParams['Istrgpi'] = {'mod': 'Exp2Syn',
                                            'tau1': 5,
                                            'tau2': 5,
                                            'e': -85}  # D1 -> gpi


# STN
netParams.synMechParams['Igesn'] = {'mod': 'Exp2Syn',
                                            'tau1': 0.4,
                                            'tau2': 7.7,
                                            'e': -85}  # gpe -< stn
netParams.synMechParams['Icosn,ampa'] = {'mod': 'Exp2Syn',
                                                'tau1': 0.5,
                                                'tau2': 2.49,
                                                'e': 0}  # ctx -> gpe
netParams.synMechParams['Icosn,nmda'] = {'mod': 'Exp2Syn',
                                                'tau1': 2,
                                                'tau2': 90,
                                                'e': 0}  # ctx -> gpe

# Str
# refer to Damodaran et al. 2015 MSN-MSN gaba synapse,FSI-MSN gaba synapse
netParams.synMechParams['Igaba_Str'] = {'mod': 'Exp2Syn',
                                            'tau1': 0.75,
                                            'tau2': 6.7,
                                            'e': -60}  # str -< str
netParams.synMechParams['Igaba_fsi'] = {'mod': 'Exp2Syn',
                                            'tau1': 1.33,
                                            'tau2': 4.0,
                                            'e': -60}  # fsi -> str
netParams.synMechParams['Glu_AMPA_th'] = {'mod':'AMPAk','gAmax':3420}
netParams.synMechParams['Glu_NMDA_th'] = {'mod':'NMDAk','gNmax':6980, 'tauon':2.25,'tauoff':150,
'mg':3.57} # th -> str

# Stimulation parameters
# fsi receieve external input(cortical input),meicated by ampa receptors
# Stimulation parameters
netParams.stimSourceParams['background_E']  = {'type': 'NetStim', 'rate': 100, 'noise': 1.0} # background inputs to Exc
netParams.stimSourceParams['background_I']  = {'type': 'NetStim', 'rate': 100, 'noise': 1.0} # background inputs to Inh
netParams.stimSourceParams['background_S']  = {'type': 'NetStim', 'rate': 100, 'noise': 1.0} # background inputs to Inh

netParams.stimTargetParams['bgE->IT'] = {'source': 'background_E', 'conds': {'cellType': ['IT']}, 
                                            'synMech': 'AMPA', 'weight': 0.35, 'delay': '2+normal(5,3)'}  
netParams.stimTargetParams['bgE->PT'] = {'source': 'background_E', 'conds': {'cellType': ['PT']}, 
                                            'synMech': 'AMPA', 'weight': 0.4, 'delay': '2+normal(5,3)'}  
netParams.stimTargetParams['bgI->PV'] = {'source': 'background_E', 'conds': {'cellType': ['PV']}, 
                                            'synMech': 'AMPA', 'weight': 0.2, 'delay': '2+normal(5,3)'}  
netParams.stimTargetParams['bgI->SOM'] = {'source': 'background_E', 'conds': {'cellType': ['SOM']}, 
                                            'synMech': 'AMPA', 'weight': 0.27, 'delay': '2+normal(5,3)'}   
netParams.stimTargetParams['bgI->StrFSI'] = {'source': 'background_E', 'conds': {'cellType': ['StrFSI']}, 
                                            'synMech': 'AMPA', 'weight': 0.12, 'delay': '2+normal(5,3)'}  
netParams.stimTargetParams['bgS->StrD1'] = {'source': 'background_S', 'conds': {'cellType': ['StrD1']}, 
                                             'synMech': 'AMPA', 'weight': 0.2, 'probability':0.4,'delay': '2+normal(5,3)'}   
netParams.stimTargetParams['bgS->StrD2'] = {'source': 'background_S', 'conds': {'cellType': ['StrD2']}, 
                                             'synMech': 'AMPA', 'weight': 0.2, 'probability':0.4,'delay': '2+normal(5,3)'} 
netParams.stimTargetParams['bgI->GPe'] = {'source': 'background_I', 'conds': {'pop': 'GPe'}, 
                                            'synMech': 'AMPA', 'weight':0.8, 'delay': '2+normal(5,3)'}  
netParams.stimTargetParams['bgI->stn'] = {'source': 'background_I', 'conds': {'pop': 'STN'}, 
                                            'synMech': 'AMPA', 'weight':1, 'delay': '2+normal(5,3)'}  


amp_gpe=0#7
amp_gpi=0.001#6
amp_stn=0#3
amp_dstr=0.001
amp_istr=0.001
amp_e5a=0
amp_e5p=0
bin_gpe = 0;
bin_gpi = 0;
bin_stn = 0;
bin_dstr = 0;
bin_istr = 0
bin_e5a= 0
bin_e5p=0
t_sim=1000

netParams.stimSourceParams['Input_th'] = {'type': 'IClamp',
                                                'delay': 0, 'dur': t_sim,
                                                'amp': bin_gpi * -1 + amp_gpi}
netParams.stimTargetParams['Input_th->th'] = {'source': 'Input_th',
                                                        'conds': {'pop': 'Th'},
                                                        'sec': 'Th',
                                                        'loc': 0}

# GPe no stimulation
netParams.stimSourceParams['Input_GPe'] = {'type': 'IClamp',
                                                'delay': 4,
                                                'dur': t_sim,
                                                'amp': bin_gpe * -1 + amp_gpe}
netParams.stimTargetParams['Input_GPe->GPe'] = {'source': 'Input_GPe',
                                                        'conds': {'pop': 'GPe'},
                                                        'sec': 'GPe',
                                                        'loc': 0}

# GPi receve stimulation
netParams.stimSourceParams['Input_GPi'] = {'type': 'IClamp',
                                                'delay': 0, 'dur': t_sim,
                                                'amp': bin_gpi * -1 + amp_gpi}
netParams.stimTargetParams['Input_GPi->GPi'] = {'source': 'Input_GPi',
                                                        'conds': {'pop': 'GPi'},
                                                        'sec': 'GPi',
                                                        'loc': 0}

# STN no stimulation
netParams.stimSourceParams['Input_STN'] = {'type': 'IClamp',
                                                'delay': 3,
                                                'dur': t_sim,
                                                'amp': bin_stn * -1 + amp_stn}
netParams.stimTargetParams['Input_STN->STN'] = {'source': 'Input_STN',
                                                        'conds': {'pop': 'STN'},
                                                        'sec': 'STN',
                                                        'loc': 0}

# dStr no stimulation
netParams.stimSourceParams['Input_StrD1'] = {'type': 'IClamp',
                                                    'delay': 0,
                                                    'dur': t_sim,
                                                    'amp': bin_dstr * -1 + amp_dstr}
netParams.stimTargetParams['Input_StrD1->StrD1'] = {'source': 'Input_StrD1',
                                                            'conds': {'pop': 'StrD1'},
                                                            'sec': 'StrD1',
                                                            'loc': 0}

# iStr  no stimulation
netParams.stimSourceParams['Input_StrD2'] = {'type': 'IClamp',
                                                    'delay': 0, 'dur': t_sim,
                                                    'amp': bin_istr * -1 + amp_istr}
netParams.stimTargetParams['Input_StrD2->StrD2'] = {'source': 'Input_StrD2',
                                                            'conds': {'pop': 'StrD2'},
                                                            'sec': 'StrD2',
                                                            'loc': 0}
# E5A  no stimulation
netParams.stimSourceParams['Input_E5A'] = {'type': 'IClamp',
                                                    'delay': 0, 'dur': t_sim,
                                                    'amp': bin_e5a * -1 + amp_e5a}
netParams.stimTargetParams['Input_E5A->E5A'] = {'source': 'Input_E5A',
                                                            'conds': {'pop': 'E5A'},
                                                            'sec': 'E5A',
                                                            'loc': 0}

# E5P  no stimulation
netParams.stimSourceParams['Input_E5P'] = {'type': 'IClamp',
                                                    'delay': 0, 'dur': t_sim,
                                                    'amp': bin_e5p * -1 + amp_e5p}
netParams.stimTargetParams['Input_E5P->E5P'] = {'source': 'Input_E5P',
                                                            'conds': {'pop': 'E5P'},
                                                            'sec': 'E5P',
                                                            'loc': 0}


scale_f = 0.75
##I2L-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'I23L'},
'postConds': {'pop': 'I23L'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})

##I2L-E2
netParams.addConnParams(None, {'preConds': {'pop': 'I23L'},
'postConds': {'pop': 'E23'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.45,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I2L-I2
netParams.addConnParams(None, {'preConds': {'pop': 'I23L'},
'postConds': {'pop': 'I23'},
'synMech': 'GABAB',
'probability': 1.0,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#I2-E2
netParams.addConnParams(None, {'preConds': {'pop': 'I23'},
'postConds': {'pop': 'E23'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.225,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#I2-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'I23'},
'postConds': {'pop': 'I23L'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I2-I2
netParams.addConnParams(None, {'preConds': {'pop': 'I23'},
'postConds': {'pop': 'I23'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E2-E2
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'E23'},
'synMech': 'AMPA',
'probability': 0.15,
'weight': 0.16,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E2-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'I23L'},
'synMech': 'AMPA',
'probability': 0.19,
'weight': 0.117,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E2-I2
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'I23'},
'synMech': 'AMPA',
'probability': 0.19,
'weight': 0.117,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E2-I5
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'I5'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.017,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E2-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'I5L'},
'synMech': 'AMPA',
'probability': 0.22,
'weight': 0.151,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E2-E5a
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'E5A'},
'synMech': 'AMPA',
'probability': 0.1,
'weight': 0.226,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E2-E5B
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'E5B'},
'synMech': 'AMPA',
'probability': 0.05,
'weight': 0.211,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E2-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'E23'},
'postConds': {'pop': 'E5P'},
'synMech': 'AMPA',
'probability': 0.11,
'weight': 0.211,
'delay': 'defaultDelay+dist_3D/propVelocity'})




#E5A-E2
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'E23'},
'synMech': 'AMPA',
'probability': 0.04,
'weight': 0.131,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5A-I2
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'I23'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5a-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'I23L'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5a-E5a
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'E5A'},
'synMech': 'AMPA',
'probability': 0.18,
'weight': 0.143,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5A-E5B
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'E5B'},
'synMech': 'AMPA',
'probability': 0.01,
'weight': 0.208,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5A-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'E5P'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.208,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5A-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'I5L'},
'synMech': 'AMPA',
'probability': 0.03,
'weight': 0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5A-I5
netParams.addConnParams(None, {'preConds': {'pop': 'E5A'},
'postConds': {'pop': 'I5'},
'synMech': 'AMPA',
'probability': 0.19,
'weight':0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5b-E2
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'E23'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.059,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5B-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'I23L'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5B-I2
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'I23'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5B-E5A
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'E5A'},
'synMech': 'AMPA',
'probability': 0.05,
'weight': 0.08 ,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5b-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'E5P'},
'synMech': 'AMPA',
'probability': 0.04,
'weight': 0.171,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5B-E5B
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'E5B'},
'synMech': 'AMPA',
'probability': 0.18,
'weight': 0.171,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5B-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'I5L'},
'synMech': 'AMPA',
'probability': 0.03,
'weight': 0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5b-I5
netParams.addConnParams(None, {'preConds': {'pop': 'E5B'},
'postConds': {'pop': 'I5'},
'synMech': 'AMPA',
'probability': 0.19,
'weight':0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})



#E5p-I2
netParams.addConnParams(None, {'preConds': {'pop': 'E5P'},
'postConds': {'pop': 'I23'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#E5P-I2L
netParams.addConnParams(None, {'preConds': {'pop': 'E5P'},
'postConds': {'pop': 'I23L'},
'synMech': 'AMPA',
'probability': 0.02,
'weight': 0.054,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5p-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'E5P'},
'postConds': {'pop': 'E5P'},
'synMech': 'AMPA',
'probability': 0.18,
'weight': 0.171,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5P-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'E5P'},
'postConds': {'pop': 'I5L'},
'synMech': 'AMPA',
'probability': 0.03,
'weight': 0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#E5P-I5
netParams.addConnParams(None, {'preConds': {'pop': 'E5P'},
'postConds': {'pop': 'I5'},
'synMech': 'AMPA',
'probability': 0.19,
'weight': 0.162,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I5L-E5A
netParams.addConnParams(None, {'preConds': {'pop': 'I5L'},
'postConds': {'pop': 'E5A'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.35,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I5L-E5b
netParams.addConnParams(None, {'preConds': {'pop': 'I5L'},
'postConds': {'pop': 'E5B'},
'synMech': 'GABAB',
'probability': 1*scale_f,
'weight': 0.35,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#I5L-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'I5L'},
'postConds': {'pop': 'E5P'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.35,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I5L-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'I5L'},
'postConds': {'pop': 'I5L'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})



#I5L-I5
netParams.addConnParams(None, {'preConds': {'pop': 'I5L'},
'postConds': {'pop': 'I5'},
'synMech': 'GABAB',
'probability': 1.0*0.3,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#I5-E5a
netParams.addConnParams(None, {'preConds': {'pop': 'I5'},
'postConds': {'pop': 'E5A'},
'synMech': 'GABAB',
'probability': 1*scale_f,
'weight': 0.2,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I5-E5B
netParams.addConnParams(None, {'preConds': {'pop': 'I5'},
'postConds': {'pop': 'E5B'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.2,
'delay': 'defaultDelay+dist_3D/propVelocity'})

#I5-E5P
netParams.addConnParams(None, {'preConds': {'pop': 'I5'},
'postConds': {'pop': 'E5P'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.2,
'delay': 'defaultDelay+dist_3D/propVelocity'})


#I5-I5L
netParams.addConnParams(None, {'preConds': {'pop': 'I5'},
'postConds': {'pop': 'I5L'},
'synMech': 'GABAB',
'probability': 1.0*0.3,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})



#I5-I5
netParams.addConnParams(None, {'preConds': {'pop': 'I5'},
'postConds': {'pop': 'I5'},
'synMech': 'GABAB',
'probability': 1.0*scale_f,
'weight': 0.15,
'delay': 'defaultDelay+dist_3D/propVelocity'})


# STN->GPe connections

netParams.connParams['STN->GPe'] = {
    'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'GPe'},  # STN-> GPe
    'probability': 0.3,
    'weight':0.05,  # synaptic weight (conductance)0.007
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Insge,ampa'}  # target synaptic mechanism

# GPe-< GPe connections
netParams.connParams['GPe->GPe'] = {
    'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'GPe'},  # GPe-< GPe
    'probability': 0.5,
    'weight': 0.05,  # synaptic weight (conductance) 0.05+0.75*pd
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igege'}  # target synaptic mechanism

# StrD2>GPe connections
netParams.connParams['StrD2->GPe'] = {
    'preConds': {'pop': 'StrD2'}, 'postConds': {'pop': 'GPe'},  # StrD2-> GPe
    'probability': 1,
    'weight': 0.02,  # synaptic weight (conductance) 0.02
    'delay':7,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Istrgpe'}  # target synaptic mechanism


netParams.connParams['STN->GPi'] = {
    'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'GPi'},
    'probability': 0.6,
    'weight': 0.03,  # synaptic weight (conductance)
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Isngi'}  # target synaptic mechanism

# GPe-< GPi connections 

netParams.connParams['GPe->GPi'] = {
    'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'GPi'},
    'probability': 0.4,
    # [ [ idx, (idx+2) % n_neurons ] for idx in range( n_neurons ) ] + \
    # [ [ (idx+1) % n_neurons, idx ] for idx in range( n_neurons ) ],
    'weight': 0.036,  # synaptic weight (conductance) 
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igegi'}  # target synaptic mechanism

# GPe -> D1 connections
netParams.connParams['GPe->D1'] = {
    'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'StrD1'},
    'probability': 0.2,
    # [ [ idx, (idx+2) % n_neurons ] for idx in range( n_neurons ) ] + \
    # [ [ (idx+1) % n_neurons, idx ] for idx in range( n_neurons ) ],
    'weight': 0.015,  # synaptic weight (conductance)
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igestr'}  # target synaptic mechanism

# GPe -> D2 connections
netParams.connParams['GPe->D2'] = {
    'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'StrD2'},
    'probability': 0.2,
    # [ [ idx, (idx+2) % n_neurons ] for idx in range( n_neurons ) ] + \
    # [ [ (idx+1) % n_neurons, idx ] for idx in range( n_neurons ) ],
    'weight': 0.018,  # synaptic weight (conductance)
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igestr'}  # target synaptic mechanism


# StrD1>GPi connections

netParams.connParams['StrD1->GPi'] = {
    'preConds': {'pop': 'StrD1'}, 'postConds': {'pop': 'GPi'},  # StrD1-> GPi
    'probability': 0.7,
    'weight': 0.025,  # synaptic weight (conductance)0.015
    'delay': 7,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Istrgpi'}  # target synaptic mechanism

# GPe-> STN connections 

netParams.connParams['GPe->STN'] = {
    'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'STN'},  # GPe-< STN
    'probability': 0.6,
    'weight': 0.02,  # synaptic weight (conductance)0.01
    'delay':6,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igesn'}  # target synaptic mechanism


# CTX-> STN connections

netParams.connParams['CTX->STN'] = {
    'preConds': {'pop': 'E5B'}, 'postConds': {'pop': 'STN'},  # CTX-> STN
    'probability': 0.2,
    'weight': 0.02,  # synaptic weight (conductance)
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Icosn,ampa'}  # target synaptic mechanism

netParams.connParams['CTX->STN2'] = {
    'preConds': {'pop': 'E5B'}, 'postConds': {'pop': 'STN'},  # CTX-> STN
    'probability': 0.2,
    'weight': 0.02,  # synaptic weight (conductance)
    'delay': 2,  # transmission delay (ms) 0.01
    'loc': 1,  # location of synapse
    'synMech': 'Icosn,nmda'}  # target synaptic mechanism


# refer to Bahuguan et al. 2015
# StrD2-< StrD2 connections
netParams.connParams['StrD2-> StrD2'] = {
    'preConds': {'pop': 'StrD2'}, 'postConds': {'pop': 'StrD2'},  # StrD2-< StrD2
    'probability': 0.36,
    'weight': 0.1,  # synaptic weight (conductance) -> mudar essa maluquisse 0.05e-2
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_Str'}  # target synaptic mechanism

# StrD1-< StrD1 connections
netParams.connParams['StrD1-> StrD1'] = {
    'preConds': {'pop': 'StrD1'}, 'postConds': {'pop': 'StrD1'},  # StrD1-< StrD1
    'probability': 0.26,
    'weight':0.05 ,  # synaptic weight (conductance) -> mudar aqui tb
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_Str'}  # target synaptic mechanism

# StrD1-< StrD2 connections
netParams.connParams['StrD1-> StrD2'] = {
    'preConds': {'pop': 'StrD1'}, 'postConds': {'pop': 'StrD2'},  # StrD1-< StrD2
    'probability': 0.07,
    'weight':0.1 ,  # synaptic weight (conductance) -> mudar aqui tb
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_Str'}  # target synaptic mechanism

# StrD2-< StrD1 connections
netParams.connParams['StrD2-> StrD1'] = {
    'preConds': {'pop': 'StrD2'}, 'postConds': {'pop': 'StrD1'},  # StrD1-< StrD1
    'probability': 0.27,
    'weight':0.12 ,  # synaptic weight (conductance) -> mudar aqui tb
    'delay': 2,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_Str'}  # target synaptic mechanism

# StrFSI-< StrD1 connections
netParams.connParams['StrFSI-> StrD1'] = {
    'preConds': {'pop': 'StrFSI'}, 'postConds': {'pop': 'StrD1'},  # StrD1-< StrD1
    'probability': 0.54,
    'weight':0.25 ,  # synaptic weight (conductance) -> mudar aqui tb
    'delay': 1,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_fsi'}  # target synaptic mechanism

# StrFSI-< StrD2 connections
netParams.connParams['StrFSI-> StrD2'] = {
    'preConds': {'pop': 'StrFSI'}, 'postConds': {'pop': 'StrD2'},  # StrD1-< StrD1
    'probability': 0.36,
    'weight':0.25 ,  # synaptic weight (conductance) -> mudar aqui tb
    'delay': 1,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igaba_fsi'}  # target synaptic mechanism

# StrFSI -> StrFSI connections(gap junction)
# refer to Damodaran et al. 2015
netParams.connParams['StrFSI->StrFSI'] = {
    'preConds': {'pop': 'StrFSI'}, 'postConds': {'pop': 'StrFSI'}, 
    'probability': 0.3,
    'weight': 0.3,  #gap junction weight does not contribute to FSI firing pattern
    'delay': 1,            
    'synMech': 'esyn',                   
    'gapJunction': True,
    'loc': 1}  

# RS-> StrD1 connections 
# refer to Deng et al. 2015
# D1: 44.7%IT type +18%PT type
# D2: 24.2%IT type + 50%PT type
netParams.connParams['RS->StrD1'] = {
    'preConds': {'pop': 'E5A'}, 'postConds': {'pop': 'StrD1'},  # RS-> StrD1
   'probability': 0.5,
    'weight': 0.68,  # synaptic weight (conductance) (0.07 - 0.044 *pd)* 0.43 e-2* gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_AMPA_d1'}  # target synaptic mechanism


netParams.connParams['RS->StrD1'] = {
    'preConds': {'pop': 'E5A'}, 'postConds': {'pop': 'StrD1'},  # RS-> StrD1
   'probability': 0.5,
    'weight': 0.68,  # synaptic weight (conductance) (0.07 - 0.044 *pd)* 0.43 e-2* gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_NMDA_d1'}  # target synaptic mechanism

netParams.connParams['RS->StrD1'] = {
    'preConds': {'pop': 'E5P'}, 'postConds': {'pop': 'StrD1'},  # RS-> StrD1
   'probability': 0.2,
    'weight': 0.68,  # synaptic weight (conductance) (0.07 - 0.044 *pd)* 0.43 e-2* gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_AMPA_d1'}  # target synaptic mechanism


netParams.connParams['RS->StrD1'] = {
    'preConds': {'pop': 'E5P'}, 'postConds': {'pop': 'StrD1'},  # RS-> StrD1
   'probability': 0.2,
    'weight': 0.68 ,  # synaptic weight (conductance) (0.07 - 0.044 *pd)* 0.43 e-2* gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_NMDA_d1'}  # target synaptic mechanism

# RS-> StrD2 connections 
# n_neurons = min(netParams.popParams['E5B']['numCells'],
#                 netParams.popParams['StrD2']['numCells'])
netParams.connParams['RS->StrD2'] = {
    'preConds': {'pop': 'E5P'}, 'postConds': {'pop': 'StrD2'},  # RS-> StrD2 
    'probability': 0.6,
    'weight':0.12,  # synaptic weight (conductance) 0.43e-2e-2 0.07  * 0.43e-2*gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_AMPA_d2'}  # target synaptic mechanism

netParams.connParams['RS->StrD2'] = {
    'preConds': {'pop': 'E5P'}, 'postConds': {'pop': 'StrD2'},  # RS-> StrD2 
    'probability': 0.6,
    'weight':0.12,  # synaptic weight (conductance) 0.43e-2e-2 0.07  * 0.43e-2*gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_NMDA_d2'}  # target synaptic mechanism

netParams.connParams['RS->StrD2'] = {
    'preConds': {'pop': 'E5A'}, 'postConds': {'pop': 'StrD2'},  # RS-> StrD2 
    'probability': 0.3,
    'weight':0.12,  # synaptic weight (conductance) 0.43e-2e-2 0.07  * 0.43e-2*gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_AMPA_d2'}  # target synaptic mechanism

netParams.connParams['RS->StrD2'] = {
    'preConds': {'pop': 'E5A'}, 'postConds': {'pop': 'StrD2'},  # RS-> StrD2 
    'probability': 0.3,
    'weight':0.12,  # synaptic weight (conductance) 0.43e-2e-2 0.07  * 0.43e-2*gsynmod
    'delay': 5,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Glu_NMDA_d2'}  # target synaptic mechanism

#the connections from STN to cortex
netParams.connParams['STN->E23'] = {
    'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'E23'},
    'probability': 0.5,
    'weight': 0.002,  # synaptic weight (conductance)
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Isngi'}  # target synaptic mechanism

netParams.connParams['STN->I23'] = {
    'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'I23'},
    'probability': 0.5,
    'weight': 0.001,  # synaptic weight (conductance)
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Isngi'}  # target synaptic mechanism

netParams.connParams['STN->I23L'] = {
    'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'I23L'},
    'probability': 0.5,
    'weight': 0.002,  # synaptic weight (conductance)
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Isngi'}  # target synaptic mechanism

n_neurons = 100
netParams.connParams['GPi->th'] = {
    'preConds': {'pop': 'GPi'}, 'postConds': {'pop': 'Th'},  # GPi-> th
    'connList': [[i, i] for i in range(n_neurons)],
    'weight': 0.0003,  # synaptic weight (conductance)
    'delay': 1,  # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Igith'}  # target synaptic mechanism

# Or probability connection
# netParams.connParams['GPi->th'] = {
#     'preConds': {'pop': 'GPi'}, 'postConds': {'pop': 'Th'},
#     'probability': 0.2,
#     'weight': 0.0001,  # synaptic weight (conductance)
#     'delay': 3 , # transmission delay (ms)
#     'loc': 1,  # location of synapse
#     'synMech':'Igith'}  # target synaptic mechanism

netParams.connParams['TH->E23'] = {
    'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'E23'},
    'probability': 0.2,
    'weight': 0.002,  # synaptic weight (conductance)
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Ithco'}  # target synaptic mechanism


netParams.connParams['TH->E5A'] = {
    'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'E5A'},
    'probability': 0.2,
    'weight': 0.002,  # synaptic weight (conductance)#0.0002
    'delay': 3 , # transmission delay (ms)
    'loc': 1,  # location of synapse
    'synMech': 'Ithco'}  # target synaptic mechanism

netParams.connParams['TH->E5B'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'E5B'},
'probability': 0.2,
'weight': 0.002,  # synaptic weight (conductance)
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Ithco'}  # target synaptic mechanism

netParams.connParams['TH->E5P'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'E5P'},
'probability': 0.2,
'weight': 0.002,  # synaptic weight (conductance)#0.005
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Ithco'}  # target synaptic mechanism

netParams.connParams['TH->StrD1'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'StrD1'},
'probability': 0.2,
'weight': 0.03,  # synaptic weight (conductance)#0.005
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Glu_AMPA_th'}  # target synaptic mechanism

netParams.connParams['TH->StrD1'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'StrD1'},
'probability': 0.2,
'weight': 0.03,  # synaptic weight (conductance)#0.005
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Glu_NMDA_th'}  # target synaptic mechanism

netParams.connParams['TH->StrD2'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'StrD2'},
'probability': 0.2,
'weight': 0.02,  # synaptic weight (conductance)#0.005
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Glu_AMPA_th'}  # target synaptic mechanism

netParams.connParams['TH->StrD2'] = {
'preConds': {'pop': 'Th'}, 'postConds': {'pop': 'StrD2'},
'probability': 0.2,
'weight': 0.02,  # synaptic weight (conductance)#0.005
'delay': 3 , # transmission delay (ms)
'loc': 1,  # location of synapse
'synMech': 'Glu_NMDA_th'}  # target synaptic mechanism

# Dictionary of annotations
netParams.annots = {}

(pops, cells, conns,  stims, simData)=sim.createSimulateAnalyze(netParams = netParams, 
simConfig = simConfig, output=True)  # create and simulate netw
matrix = analysis.spikes_legacy.calculateRate(include=['allCells', 'eachPop'], peakBin=5, timeRange=None)
print(matrix)
sim.gather.gatherData(gatherLFP=True)
# retrieve the LFP data
lfp_data = sim.allSimData['LFP']
filename_lfp = 'lfp_sham_.txt'
np.savetxt(filename_lfp, lfp_data, delimiter='\t')

netpyne.plotting.plotRaster()

'''
netpyne.analysis.plotTraces([('I23', [17,26])],oneFigPer='cell')

v_trace_data = sim.allSimData['V']
v_list=[17,45,100,294,374,430,636,730,834,954,1104,1194,1304,1394]
filename_v = 'v_data_pd.txt'
cell_data_list = []
for v in v_list:  
    key = f'cell_{v}'
    if key in v_trace_data:
        cell_data_list.append(v_trace_data[key])

cell_data_array = np.array(cell_data_list)

np.savetxt(filename_v, cell_data_array, delimiter='\t')
'''
