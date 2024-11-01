COMMENT

Chloride accumulation and diffusion with decay (time constant tau) to resting level cli0.
The decay approximates a reversible chloride pump with first order kinetics.
To eliminate the chloride pump, just use this hoc statement
To make the time constant effectively "infinite".
tau and the resting level are both RANGE variables

Diffusion model is modified from Ca diffusion model in Hines & Carnevale: 
Expanding NEURON with NMODL, Neural Computation 12: 839-851, 2000 (Example 8)

Nannuli correspond to the number of shells*

ENDCOMMENT

NEURON {
	SUFFIX iondifus
	USEION cl READ icl WRITE cli, clo VALENCE -1
	USEION hco3 READ hco3i, hco3o VALENCE -1
	USEION na READ ina WRITE nai, nao VALENCE 1
	USEION k READ ik WRITE ki VALENCE 1
	USEION gab READ gabo WRITE gabo VALENCE 0
	USEION mess READ messi VALENCE 0
	GLOBAL vrat, DGab, taugaba, fhspace, areaext		:vrat must be GLOBAL
	GLOBAL tau, clipip, kipip, naipip
	RANGE cli0, clo0, nai0, nao0, ki0, ko0, egaba, delta_egaba, init_egaba, ehco3_help, ecl_help
	RANGE gabo0, temp, DCl
	RANGE clamp, vclamp, hco3i0, hco3o0
}

DEFINE Nannuli 4

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	(mV)    = (millivolt)
	FARADAY = (faraday) (10000 coulomb)
	PI = (pi) (1) : le 1 entre parenthese est necessaire pour dire que c'est adimensionnel
	F = (faraday) (coulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	bath 	= 0		(mM)
	DCl 	= 2 	(um2/ms) : Kuner & Augustine, Neuron 27: 447
	DK 		= 1.96	(um2/ms)
	DNa	 	= 1.3   (um2/ms) :0.1 in original
	DGab 	= 0.6	(um2/ms)
	fhspace = 0.03	(um)
	tau 	= 100 	(ms)
	cli0 	= 3.5 	(mM) : 8 mM in original file
	clo0 	= 133.5 (mM)
	nai0 	= 10 	(mM)
	gabo0 	= 0 	(mM)
	nao0 	= 147.5 (mM)
	ki0 	= 135 	(mM)
	ko0 	= 3.5 	(mM)
	hco3i0 	= 16	(mM)
	hco3o0 	= 26	(mM)
	P_help 	= 0.18
	celsius = 37    (degC)
	taugaba = 100	(ms)
	clamp   = 0
	vclamp  = 1000	(mV)
	clipip  = 8		(mM)
	kipip   = 140	(mM)
	naipip  = 12	(mM)
}

ASSIGNED {
	v 		(mV)
	diam 	(um)
	icl 	(mA/cm2)
	cli 	(mM)
	clo		(mM)
	ina 	(mA/cm2)
	nai 	(mM)
	nao		(mM)
	ik 	(mA/cm2)
	ki 	(mM)
	ko		(mM)
	hco3i	(mM)
	hco3o	(mM)
	vrat[Nannuli]	: numeric value of vrat[i] equals the volume
			: of annulus i of a 1um diameter cylinder
			: multiply by diam^2 to get volume per um length
	areaext
	areapip
	egaba 	(mV)
	ehco3_help 	(mV)
	ecl_help	(mV)
	init_egaba  (mV)
	delta_egaba (mV)
	messi 		(mM)
	temp
}

STATE {
	: cl[0] is equivalent to cli
	: cl[] are very small, so specify absolute tolerance
	cl[Nannuli]	(mM) <1e-10>
	na[Nannuli]	(mM) <1e-10>
	k[Nannuli]	(mM) <1e-10>
	gabo		(mM) <1e-10>
}

BREAKPOINT { 
		if (messi > 0) {
			if (temp == 0) {
				gabo = messi
				temp = 1
			}
		}
		
		SOLVE state METHOD sparse
		ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
		egaba = P_help*ehco3_help + (1-P_help)*ecl_help
		delta_egaba = egaba - init_egaba
}

LOCAL factors_done

INITIAL {
	if (factors_done == 0) {  	: flag becomes 1 in the first segment	
		factors_done = 1	: all subsequent segments will have
		factors()		: vrat = 0 unless vrat is GLOBAL
	}

	temp = 0
	messi = 0
	cli = cli0
	clo = clo0
	nai = nai0
	nao = nao0
	ki = ki0
	:ko = ko0
	gabo = 0
	hco3i = hco3i0
	hco3o = hco3o0
	FROM i=0 TO Nannuli-1 {
		cl[i] = cli
		na[i] = nai
		k[i] = ki
	}

	ehco3_help = log(hco3i/hco3o)*(1000)*(celsius + 273.15)*R/F
	ecl_help = log(cli/clo0)*(1000)*(celsius + 273.15)*R/F
	egaba = P_help*ehco3_help + (1-P_help)*ecl_help
	init_egaba = egaba
	delta_egaba = egaba - init_egaba 
}

LOCAL frat[Nannuli]	: scales the rate constants for model geometry

PROCEDURE factors() {
	LOCAL r, dr2, rgab
	r = 1/2			: starts at edge (half diam), diam = 1, length = 1
	dr2 = r/(Nannuli-1)/2	: full thickness of outermost annulus,
				: half thickness of all other annuli
	vrat[0] = 0
	frat[0] = 2*r		: = diam

	FROM i=0 TO Nannuli-2 {
		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2	: interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	: outer radius of annulus Ai+1/delta_r=2PI*r*1/delta_r
						: div by distance between centers 
		r = r - dr2
		vrat[i+1] = PI*(r+dr2/2)*2*dr2	: outer half of annulus
	}

	rgab = (diam+fhspace)/2
	areaext = PI*rgab*rgab - PI*(diam*diam)/4
	areapip = PI*2*2
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

KINETIC state {
	COMPARTMENT areaext {gabo}
	LONGITUDINAL_DIFFUSION DGab*areaext  {gabo}
	~ gabo <-> bath (1/(taugaba), 1/(taugaba))

	COMPARTMENT i, diam*diam*vrat[i] {cl na k}
	LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}
	LONGITUDINAL_DIFFUSION i, DNa*diam*diam*vrat[i] {na}
	LONGITUDINAL_DIFFUSION i, DK*diam*diam*vrat[i] {k}

	~ cl[0] << ((icl*PI*diam/FARADAY))
	~ na[0] << ((-ina*PI*diam/FARADAY))
	~ k[0] << ((-ik*PI*diam/FARADAY))

	FROM i=0 TO Nannuli-2 {
		~ cl[i] <-> cl[i+1]	(DCl*frat[i+1], DCl*frat[i+1])
		~ na[i] <-> na[i+1]	(DNa*frat[i+1], DNa*frat[i+1])
		~ k[i] <-> k[i+1]	(DK*frat[i+1], DK*frat[i+1])
	}

	if (clamp == 1) {
		~ cl[0] <-> clipip (1/(tau), 1/(tau))
		~ na[0] <-> naipip (1/(tau), 1/(tau))
		~ k[0] <-> kipip (1/(tau), 1/(tau))
		cli = cl[0]
		nai = na[0]
		ki  = k[0]
	}
	else {
		cli = cl[0]
		nai = na[0]
		ki  = k[0]
	}

	clo = clo0
	nao = nao0
	:ko = ko0
}
