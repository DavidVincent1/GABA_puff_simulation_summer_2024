TITLE CLC2 channels

COMMENT
 Equation comes from Ratté S, Prescott SA. (2011)
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

? interface
NEURON {
        SUFFIX clc2
        USEION cl READ ecl WRITE icl
        RANGE gclc2, vhalf, vslope, ptau, icl
        : `GLOBAL minf` will be replaced with `RANGE minf` if CoreNEURON enabled
        GLOBAL pinf
        THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        gclc2 = .001 (S/cm2)    <0,1e9> : .040
        ptau = 300   (ms)
        vhalf = 15   (mV)
        vslope = -14 (mV)
}

STATE {
        p
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ecl (mV)

        gclc (S/cm2)
        icl (mA/cm2)
        pinf
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gclc = gclc2*p
        icl = gclc*(v - ecl)
}


INITIAL {
        rates(v)
        p = pinf
}

? states
DERIVATIVE states {
        rates(v)
        p' =  (pinf-p)/ptau
}


? rates
PROCEDURE rates(v(mV)) {
        :TABLE pinf DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        :q10 = 3^((celsius - 23)/10)

                :"p" clc2 activation system
        pinf = 1 / (1 + exp((ecl - vhalf - v) / vslope))
}
UNITSON