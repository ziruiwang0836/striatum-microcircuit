COMMENT

add CA dependence?

I_K.mod as K

Hodgkin-Huxley like potassium current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_K_STN
	NONSPECIFIC_CURRENT I
	RANGE g0, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
		
	g0=4.5 		(S/cm2)
	v0=-80 		(mV)
    ik=0.0   (mA/cm2)	
	phi_n=0.75
	theta_tn=-80
	theta_n=-32
	sigma_tn=-26
	sigma_n=8
	tau_0n=1	(ms)
	tau_1n=100 	(ms)
}

ASSIGNED {	
   v 		(mV)
	I 		(mA/cm2)
	n_inf
	tau_n	(ms)
}

STATE {
	n
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(n*n*n*n)*(v-v0)
	ik=I
}

UNITSOFF

INITIAL {
	rates(v)
	n = n_inf
}

DERIVATIVE states {
	rates(v)
	n' = phi_n*((n_inf-n)/tau_n)
}




PROCEDURE rates(v) {  :Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
   TABLE tau_n, n_inf  DEPEND  tau_0n, tau_1n, theta_tn, sigma_tn, theta_n, sigma_n FROM -100 TO 100 WITH 400
	tau_n = tau_0n + tau_1n/(1+exp(-(v-theta_tn)/sigma_tn))
	n_inf = 1/(1+exp(-(v-theta_n)/sigma_n))
}

UNITSON