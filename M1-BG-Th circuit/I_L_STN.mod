COMMENT

I_L.mod as l

Hodgkin-Huxley like leak current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76

ENDCOMMENT

NEURON {
	SUFFIX I_L_STN
	NONSPECIFIC_CURRENT I
	RANGE g0, v0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {

	g0=0.225 	(S/cm2)
	v0=-60 	(mV)
}

ASSIGNED {
	v 	(mV)
	I 	(mA/cm2)
}

BREAKPOINT {
 	I=g0*(v-v0)
}