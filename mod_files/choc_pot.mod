TITLE chocpotassique.mod


NEURON {
	POINT_PROCESS CHOCpot
	USEION k READ ko WRITE ko
	RANGE tauchoc, ko, kchoc, tchoc, ko0
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	kchoc (mM)
	tchoc (ms)
}


ASSIGNED {
	tauchoc  (ms)
}


INITIAL {
	ko = ko0
}


BREAKPOINT {
	at_time(tchoc)
	if (t >= tchoc) {
		SOLVE state METHOD sparse
	}
	else {
		ko = ko0
	}
}

STATE { ko (mM) }


KINETIC state {
	~ ko <-> kchoc (1/(tauchoc), 1/(tauchoc))
}