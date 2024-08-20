COMMENT

I_AHP.mod as AHP

Hodgkin-Huxley like sodium current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76POINTER I_T, I_Ca

ENDCOMMENT

NEURON {
	SUFFIX I_AHP_STN
	NONSPECIFIC_CURRENT I
	RANGE g0, v0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	
	g0=0.99		(S/cm2)
	v0= -80		(mV)

	Ca0=0.3
	k1=15.0		
	kca =22.5
	e=5e-5
	Iscale = -10
}

ASSIGNED {
	v 		(mV)
	I 		(mA/cm2)
	I_T		(mA/cm2)
	I_Ca 	(mA/cm2)
	
}

STATE {
	Ca		
}

UNITSOFF

INITIAL {
	
	Ca = Ca0
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(v-v0)*(Ca/(Ca+k1))
}


DERIVATIVE states {
	Ca' = e*(Iscale*(I_Ca+I_T)-kca*Ca)
}



UNITSON