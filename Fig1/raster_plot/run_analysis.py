import numpy as np
import matplotlib.cm as cm
import pickle
import sys
import itertools
import matplotlib
import pylab as pl
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import params_sham as par
#import params_pd as par
import straitum as st

fontP = FontProperties()

p = par.get_parameters()
postfix = "Sham"
#postfix = "PD"
st.run(postfix,p)


# Plot the Rasters 
evsAll = pickle.load(open("Senders_Sham.pickle","rb"))
tsAll = pickle.load(open("Times_Sham.pickle","rb"))
Rates = pickle.load(open("Rates_Sham.pickle","rb"))

evs = dict()
ts = dict()

evs['d1'] = np.array(evsAll['d1'])
evs['d2'] = np.array(evsAll['d2'])
evs['fsi'] = np.array(evsAll['fsi'])

ts['d1'] =  np.array(tsAll['d1'])
ts['d2'] =  np.array(tsAll['d2'])
ts['fsi'] = np.array(tsAll['fsi'])


# Print a smaller snapshot of simulation, 100ms
stoptime=600
binsize=2.
starttime = 500
binning = np.arange(starttime,stoptime,binsize)
fig, ax1 = plt.subplots(1, 1)
ax2 = ax1.twinx()
# Raster plot 
#d1
ind = np.where(ts['d1']<=stoptime) and np.where(ts['d1']>starttime)
semi_ts = ts['d1'][ind]
semi_evs = evs['d1'][ind]
ax1.plot(semi_ts,semi_evs,'.',color='#B9D2BD',label='D1',markersize=4)
#d2
ind = np.where(ts['d2']<=stoptime) and np.where(ts['d2']>starttime)
semi_ts = ts['d2'][ind]
semi_evs = evs['d2'][ind]
ax1.plot(semi_ts,semi_evs,'.',color='#658C8D',label='D2',markersize=4)
#plotting parameters
ax1.set_yticks([1,2000,4000])
ax1.set_yticklabels(['0','2k','4k'],fontsize=10,fontweight='bold')#,fontsize=8,stretch='ultra-condensed')
for i in ax1.get_xticklabels():
    i.set_fontsize(10)
    i.set_visible(False)
    i.set_fontweight('bold')
ax1.set_ylabel(' Neuron ID',fontsize=12,fontweight='bold',fontname='Computer Modern')
ax1.set_xlim(starttime,stoptime)
# Mean firing rate plot   
divider = make_axes_locatable(ax1)  
ax2 = divider.append_axes("bottom", size="100%", pad=0.0)  # To makes subplots stick to each other but not others
ind = np.where(ts['d1']<=stoptime) and np.where(ts['d1'] > starttime)
a1,b1 = np.histogram(ts['d1'][ind],binning)#d1
ind = np.where(ts['d2']<=stoptime) and np.where(ts['d2'] > starttime)
a2,b2 = np.histogram(ts['d2'][ind],binning)#d2
sizeInsecs = binsize/1000.
ax2.step(b1[:-1],(a1/2000.)/sizeInsecs,'-',color='#B9D2BD',label='D1',linewidth=1.5)
ax2.step(b2[:-1],(a2/2000.)/sizeInsecs,'-',color='#658C8D',label='D2',linewidth=1.5)
ax2.set_xlim(starttime,stoptime)
for i in ax2.get_xticklabels():
    i.set_fontsize(10)
    i.set_visible(True)
    i.set_fontweight('bold')
for i in ax2.get_yticklabels():
    i.set_visible(True)
    i.set_fontweight('bold')
#ax2.set_ylim(0,5)
ax2.set_xlabel('Time(ms)',fontsize=12,fontweight='bold',fontname='Computer Modern')
ax2.set_ylabel('Rate(Hz)',fontsize=12,fontweight='bold',fontname='Computer Modern')
#ax2.legend(('D1','D2'),loc='upper left',prop=fontP)
fig.show()
fig.savefig('raster_sham.pdf',dpi=500)