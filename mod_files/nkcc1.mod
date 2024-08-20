TITLE N-K-Cl cotransporter NKCC1

NEURON {
	SUFFIX nkcc1
	USEION k READ ko, ki WRITE ik VALENCE 1
	USEION na READ nao, nai WRITE ina VALENCE 1
	USEION cl READ clo, cli WRITE icl VALENCE -1
	RANGE ina, ik, icl, S, Vi, U
	:GLOBAL U
}

UNITS {
	(mV)	= (millivolt)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(mol)	= (1)
	FARADAY	= 96485.309 (coul/mole)
	PI	= (pi) (1)
}

PARAMETER {
	S = 5654.87 (um2)
	Vi = 8913.48 (um3)
	U = 0.0001    (mM/ms)
}

ASSIGNED {
    ina		(mA/cm2)
	ik		(mA/cm2)
	icl		(mA/cm2)
    nao     (mM)
    nai     (mM)
	ko		(mM)
	ki		(mM)
  	clo   (mM)
  	cli   (mM)
}

BREAKPOINT {
	LOCAL rate
	rate = pumprate(nai,nao,ki,ko,cli,clo)
    ina  = rate
	ik   = rate
	icl  = -2*rate
}

FUNCTION pumprate(nai,nao,ki,ko,cli,clo) {
	pumprate = U*log((nai*ki*cli)/(nao*ko*clo))*(FARADAY*Vi/(S*1e4)) : unite verifiee
}
