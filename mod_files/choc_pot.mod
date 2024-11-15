TITLE chocpotassique.mod


NEURON {
	POINT_PROCESS CHOCpot
	USEION k READ ko WRITE ko
	RANGE tauchoc, ko, kchoc, tchoc, ko0, tdur
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	tdur  (ms)
	kchoc (mM)
	tchoc (ms)
	ko0	  (ms)
}


ASSIGNED {
	tauchoc   (ms)
	kchoctemp (mM)
	t2		  (ms)
}


INITIAL {
	ko = ko0
	kchoctemp = kchoc
	t2 = tchoc + tdur
}


BREAKPOINT {
	at_time(tchoc)
	at_time(t2)

	if (t < tchoc) {
		ko = ko0
	}
	else {
		if (t < t2) {
			kchoctemp = kchoc
		}
		else {
			kchoctemp = ko0
		}
		SOLVE state METHOD sparse
	}

	:if (t >= tchoc && t <= t2) {
	:	SOLVE state METHOD sparse
	:}
	:else if (t > t2) {
	:	SOLVE stateFinal METHOD sparse
	:}
	:else {
	:	ko = ko0
	:}
}

STATE { ko (mM) }


KINETIC state {
	~ ko <-> kchoctemp (1/(tauchoc), 1/(tauchoc))
}

:KINETIC stateFinal {
:	~ ko <-> ko0 (1/(tauchoc), 1/(tauchoc))
:}