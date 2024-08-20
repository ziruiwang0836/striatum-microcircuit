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
	SUFFIX cor_IB_ICA
	USEION ca WRITE ica				: Using k ion, treat the reversal potential as a parameter and write to ik so the total k current can be tracked
	RANGE g_Ca, i_Ca					: Potassium current, specific conductance and equilibrium potential
}

PARAMETER {
	eca = 120 (mV)
	i_Ca = 0.0 (mA/cm2)				: Parameter to record this current separately to total sodium current
	g_Ca = 1e-4 (S/cm2)
	tau_max = 608 (ms)
}

ASSIGNED {
	v (mV)
	ica (mA/cm2)
	alpha_q
	beta_q
	alpha_r
	beta_r
 }

STATE {
	q r
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = g_Ca*q*q*r*(v - eca)
	i_Ca = ica 						: Record i_Na (just this sodium current) to check it is working
}

UNITSOFF

INITIAL {
	settables(v)					: ** Need to double check these intials are correct
	q= 0
	r = 0
}

DERIVATIVE states {
	settables(v)
	q' = alpha_q*(1-q)-beta_q*q
	r' = alpha_r*(1-r)-beta_r*r
}

PROCEDURE settables(v) {
	TABLE alpha_q, beta_q, alpha_r, beta_r DEPEND tau_max FROM -100 TO 100 WITH 400
	
	alpha_q = 0.055*vtrap((-27-v),3.8)
	beta_q = 0.94*exp((-75-v)/17)
	alpha_r = 0.000457*exp(-(13+v)/50)
	beta_r = 0.0065/(exp(-(15+v)/28)+1)
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON
