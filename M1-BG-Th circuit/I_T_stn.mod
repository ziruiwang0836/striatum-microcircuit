COMMENT

I_T_stn.mod as Tstn

Hodgkin-Huxley like T current for STN cells 0as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_T_stn
	NONSPECIFIC_CURRENT I
	RANGE  i_t, g0, v0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	 
	g0 =0.05		(S/cm2)
	v0 =140		(mV)
    i_t = 0.0 (mA/cm2)	
	phi_r=0.2
	theta_tr=68
	theta_r=-67
	sigma_tr=-2.2
	sigma_r=-2.0
	tau_0r=40	(ms)
	tau_1r =17.5	(ms)

	theta_a=-63
	sigma_a=7.8

	theta_b=0.4
	sigma_b=-0.1
}

ASSIGNED {
    v 		(mV)
	I 		(mA/cm2)
	r_inf
	tau_r	(ms)
	a_inf
	b_inf
}

STATE {
	r
}
UNITSOFF

INITIAL {
	rates(v)
	r = r_inf
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

	b_inf = 1/(1+exp((r-theta_b)/sigma_b))-1/(1+exp(-theta_b/sigma_b))

 	I=g0*(a_inf*a_inf*a_inf)*(b_inf*b_inf)*(v-v0)
}

DERIVATIVE states {
	rates(v)
	r' = phi_r*((r_inf-r)/tau_r)
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
    TABLE tau_r , r_inf, a_inf  DEPEND tau_0r, tau_1r, theta_tr, theta_r, sigma_tr, sigma_r, theta_a,sigma_a FROM -100 TO 100 WITH 400

	tau_r = tau_0r + tau_1r/(1+exp(-(v-theta_tr)/sigma_tr))
	r_inf = 1/(1+exp(-(v-theta_r)/sigma_r))

	a_inf = 1/(1+exp(-(v-theta_a)/sigma_a))
}

UNITSON