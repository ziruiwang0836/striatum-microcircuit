# striatum-microcircuit
striatum-microcircuit using I-F neurons

first install nest-simulator

1.The connectivity of the striatal-microcircuit is based on Yim et al. 2011

2.Simulation of PD conncectivity alternations is  based on Damodaran et al. 2015

3.Cortical MIP input is based on Yim et al. 2011 and Yim et al. 2014


1.params.py define the params used in different settings(Control, PD, PD+Glu Inh, Control + Glu Exc)

2.Striatum.py define the striatal micro-circuit structure

3.run_simulation.py runs the simulation and save the output into output_dir

4.test_run_analysis.ipynb analysis the  output .pickle file and plot the results 
