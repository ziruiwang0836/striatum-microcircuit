# striatum-microcircuit
## striatum-microcircuit using I-F neurons

The code of M1-BG-Th circuit is based on the work of Yim et al. 2011

first install nest-simulator following https://nest-simulator.readthedocs.io/en/stable/installation/index.html

### Reference:

1.The connectivity of the striatal-microcircuit is based on Yim et al. 2011 and Bahuguna et al. 2015

2.Simulation of PD conncectivity alternations is  based on Damodaran et al. 2015

3.Cortical MIP input is based on Yim et al. 2011 and Yim et al. 2014

### Instructions to setup the model
1.params.py define the params used in different settings(Control, PD, PD+Glu Inh, Control + Glu Exc)

2.Striatum.py define the striatal micro-circuit structure

3.run_simulation.py runs the simulation and save the output into output_dir

4.test_run_analysis.ipynb analysis the  output .pickle file and plot the results 

## M1-BG-Th circuit using H-H neuron

The code of M1-BG-Th circuit is based on the work of Yu et al. 2022

first install the netpyne and neuron following https://www.netpyne.org/documentation/installation


### Reference:

1. The connectivity of the M1 layer2/3 and layer 5 is based on the work of Neymotin SA 2016 and Dura-Bernal S et al 2019. 

2. Cortico-striatal DA-dependednt Glu modulation is based on Humphries et al. 2009, Evans et al. 2012 and Lindahl et al. 2016

3. The connectivity of striatum is based on Yim et al. 2011, Bahuguna et al. 2015 and Damodaran et al. 2015

4. Simualtion of PD state include the alternations within striatal connectivity (Damodaran et al. 2015) and GPe-STN loop (Yu et al. 2022), and the strengthening of cortico-STN connection (West et al. 2018)

### Instructions to setup the model

1. The nrnivmodl command will produce an architecture-dependent folder with a script called special. On 64 bit systems the folder is x86_64.

2. params.py define the params used in different settings(Control, PD, PD+Glu Inh, Control + Glu Exc)

3. run_simulation.py runs the simulation and save the output into output_dir

4. test_run_analysis.ipynb analysis the  output .pickle file and plot the results 

PAC analysis related code refer to Tort et al.:
tortlab/phase-amplitude-coupling: Matlab routines for computing the Modulation Index and Comodulogram, as described in Tort et al., J Neurophysiol 2010 (github.com)



