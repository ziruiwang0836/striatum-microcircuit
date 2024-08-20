COMMENT

I_Na.mod as Na

Hodgkin-Huxley like sodium current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_Na_STN
	NONSPECIFIC_CURRENT I
	RANGE g0, v0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
    
	g0=3.75 		(S/cm2)
	v0=55 		(mV)

	phi_h=0.75
	theta_th=-57
	theta_h=-39
	sigma_th=-3.0
	sigma_h=-3.1
	tau_0h=1.0	(ms)
	tau_1h=500 	(ms)

	theta_m=-30
	sigma_m=15
}

ASSIGNED {
	v 		(mV)
	I 		(mA/cm2)
	h_inf
	tau_h	(ms)
	m_inf
	
}

STATE {
	h
}

UNITSOFF

INITIAL {
	rates(v)
	h = h_inf
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(m_inf*m_inf*m_inf)*h*(v-v0)
}

DERIVATIVE states {
	rates(v)
	h' = phi_h*((h_inf-h)/tau_h)
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
TABLE tau_h , h_inf, m_inf DEPEND tau_0h, tau_1h,  theta_th, theta_h, sigma_th, sigma_h, theta_m, sigma_m FROM -100 TO 100 WITH 400

	tau_h = tau_0h + tau_1h/(1+exp(-(v-theta_th)/sigma_th))
	h_inf = 1/(1+exp(-(v-theta_h)/sigma_h))

	m_inf = 1/(1+exp(-(v-theta_m)/sigma_m))
}

UNITSON