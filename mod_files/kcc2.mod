TITLE K-Cl cotransporter KCC2

NEURON {
	SUFFIX kcc2
	USEION k READ ko, ki WRITE ik VALENCE 1
	USEION cl READ clo, cli WRITE icl VALENCE -1
	RANGE ik, icl, S, Vi, U
	:GLOBAL U
}

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	FARADAY	= 96485.309 (faraday)
	PI	= (pi) (1)
}

PARAMETER {
	S = 5654.87 (um2)
  	Vi = 8913.48 (um3)
	U = 0.0003    (mM/ms)
}

ASSIGNED {
	ik		(mA/cm2)
	icl		(mA/cm2)
	ko		(mM)
	ki		(mM)
  	clo   (mM)
  	cli   (mM)
}

BREAKPOINT {
	LOCAL rate
	rate = pumprate(ki,ko,cli,clo)
	ik =  rate
	icl = -rate
}

FUNCTION pumprate(ki,ko,cli,clo) {
	pumprate = U*log((ki*cli)/(ko*clo))*(FARADAY*Vi/(S*1e4)) : unite verifiee
}
