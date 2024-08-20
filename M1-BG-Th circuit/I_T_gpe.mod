COMMENT

I_T_gpe.mod as Tgpe

Hodgkin-Huxley like T current for STN cells as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_T_gpe
	NONSPECIFIC_CURRENT I
	RANGE g0, v0, tau_r, phi_r, theta_r, sigma_r, theta_a, sigma_a,r,a_inf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v 		(mV)
	g0 	=0.05	(S/cm2)
	v0 	=120	(mV)

	phi_r=1.0
	theta_tr=-70
	theta_r=-70
	sigma_tr=-2.0
	sigma_r=-2.0
    tau_r	=30  (ms)

	theta_a =-57
	sigma_a=2.0
}

ASSIGNED {
	I 		(mA/cm2)
	r_inf
	a_inf
}

STATE {
	r
}

INITIAL {
	rates(v)
	r = r_inf
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(a_inf*a_inf*a_inf)*r*(v-v0)
}

DERIVATIVE states {
	rates(v)
	r' = phi_r*((r_inf-r)/tau_r)
}

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.

UNITSOFF
	r_inf = 1/(1+exp(-(v-theta_r)/sigma_r))

	a_inf = 1/(1+exp(-(v-theta_a)/sigma_a))
}

UNITSON