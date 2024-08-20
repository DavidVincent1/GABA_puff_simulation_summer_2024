TITLE leak current :passive electrical properties

NEURON {
	SUFFIX leak
	USEION k READ ek WRITE ik VALENCE 1
	USEION na READ ena WRITE ina VALENCE 1
	USEION cl READ ecl WRITE icl VALENCE -1
    NONSPECIFIC_CURRENT ifix
	RANGE gk, ik, gna, ina, gcl, icl, ecl, gfix, ifix, gnaother
    RANGE qk, qna, qcl
}

UNITS { 
	(mV) = (millivolt)
    (mA) = (milliamp)
    PI		= (pi) (1)
	FARADAY		= 96485.309 (coul/mole)
}

PARAMETER {
	gk	 = 5e-5 (mho/cm2)  :potassium leak conductance
	gna	 = 1e-5 (mho/cm2)  :sodium leak conductance : 1
	gnaother = 1e-5 (mho/cm2) :sodium other conductance
	gcl	 = 0.05e-5 (mho/cm2)  :chloride leak conductance : 1
	gfix = 0 (mho/cm2)  :fixed leak conductance
}

ASSIGNED {
	v (mV)
    v_init (mV)
	ik (mA/cm2)
	ek (mV)
	ina (mA/cm2)
	ena (mV)
	icl (mA/cm2)
	ecl (mV)
	ifix (mA/cm2) :fixed leak current
    diam (um)
}

BREAKPOINT {
	ik = gk*(v-ek)
	ina = gna*(v-ena) + gnaother*(v-ena)
	icl = gcl*(v-ecl)
	:ifix = gfix*(v-v_init)
	ifix = gfix*(v+70)
    SOLVE integrate METHOD sparse
}

STATE { qk qna qcl}

INITIAL {
	ik = gk*(v-ek)
	ina = gna*(v-ena) + gnaother*(v-ena)
	icl = gcl*(v-ecl)
	qk = 0
	qna = 0
	qcl = 0
}

KINETIC integrate {
	COMPARTMENT diam*diam*PI/4 { qna qk qcl }
	~ qna << ((-ina*diam)*PI*(1e4)/FARADAY )
	~ qk  << ((-ik*diam)*PI*(1e4)/FARADAY )
	~ qcl << ((-icl*diam)*PI*(1e4)/FARADAY )
}
