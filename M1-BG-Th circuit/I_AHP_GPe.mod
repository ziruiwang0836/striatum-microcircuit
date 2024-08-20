COMMENT

I_AHP.mod as AHP

Hodgkin-Huxley like sodium current as decribed in:
Terman D, Rubin JE, Yew AC, Wilson CJ (2002) Activity patterns in a model
for the subthalamopallidal network of the basal ganglia. J Neurosci 22:2963-76	POINTER I_T, I_Ca

ENDCOMMENT

NEURON {
	SUFFIX  I_AHP_GPe
	NONSPECIFIC_CURRENT I
    RANGE g0, v0, Ca0, k1, kca, epsilon, Ca
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v 		(mV)
	g0=3 		(S/cm2)
	v0 =-80		(mV)

	Ca0=0.02
	k1	=30	
	kca =20	
	epsilon=1e-4
	Iscale = 10
}

ASSIGNED {
	I 		(mA/cm2)
	I_T		(mA/cm2)
	I_Ca 	(mA/cm2)
}

STATE {
	Ca		
}

INITIAL {
	Ca = Ca0
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(v-v0)*(Ca/(Ca+k1))
}

UNITSOFF

DERIVATIVE states {
	Ca' = epsilon*(-Iscale*(I_Ca+I_T)-kca*Ca)
}

UNITSON