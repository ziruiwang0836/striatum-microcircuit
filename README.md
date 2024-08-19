# striatum-microcircuit
striatum-microcircuit using I-F neurons

1. first install nest-simulator
The connectivity of the striatal-microcircuit is based on Yim et al. 2011
Simulation of PD conncectivity alternations is  based on Damodaran et al. 2015
Cortical MIP input is based on Yim et al. 2011 and Yim et al. 2014

params.py define the params used in different settings(Control, PD, PD+Glu Inh, Control + Glu Exc)
Striatum.py define the striatal micro-circuit structure
run_simulation.py runs the simulation and save the output into output_dir
test_run_analysis.ipynb analysis the  output .pickle file and plot the results 
