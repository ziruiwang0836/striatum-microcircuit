COMMENT

I_Ca.mod as Ca

Hodgkin-Huxley like calcium current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_Ca_STN
	NONSPECIFIC_CURRENT I
	RANGE g0, v0
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	
	g0=0.05 		(S/cm2)
	v0= 140		(mV)

	theta_s=-39
	sigma_s=8.0
}

ASSIGNED {
	v 		(mV)
	I 		(mA/cm2)
	s_inf
}

UNITSOFF

INITIAL {
	rates(v)
}

BREAKPOINT {
	rates(v)

 	I=g0*s_inf*s_inf*(v-v0)
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
   TABLE s_inf  DEPEND theta_s, sigma_s FROM -100 TO 100 WITH 400

	s_inf = 1/(1+exp(-(v-theta_s)/sigma_s))
}

UNITSON