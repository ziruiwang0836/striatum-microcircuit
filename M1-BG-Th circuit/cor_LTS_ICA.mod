TITLE Delayed-rectifier Potassium Current for Cortical Neuron Soma

COMMENT
  
  Model Reference: 
  
  Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
  Bal, T., Frégnac, Y., Markram, H. and Destexhe, A., 2008. 
  "Minimal Hodgkin–Huxley type models for different classes of 
  cortical and thalamic neurons." 
  Biological cybernetics, 99(4-5), pp.427-441.
  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cor_LTS_ICA
	USEION ca WRITE ica				: Using k ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_Ca, i_Ca					: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	eca = 120 (mV)
	i_Ca = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_Ca = 4e-4 (S/cm2)
	V_x = -7(mV)
}

ASSIGNED {
	v (mV)
	ica (mA/cm2)
	s_inf
	u_inf
	tau_u (ms)
}

STATE {
	u
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica= g_Ca*s_inf*s_inf*u*(v - eca)
	i_Ca = ica						: Record i_M (just this potassium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)					
	u = u_inf
}

DERIVATIVE states {
	settables(v)
	u' = (u_inf - u)/tau_u
}

PROCEDURE settables(v) {
	TABLE s_inf, u_inf, tau_u DEPEND V_x FROM -100 TO 100 WITH 400
	
	s_inf= 1/(1+exp(-(v+V_x+57)/6.2))
    u_inf= 1/(1+exp((v+V_x+81)/4))
	tau_u = (30.8+(211.4+exp((v+V_x+113.2)/5)))/(3.7*(1+exp((v+V_x+84)/3.2)))
}

UNITSON 