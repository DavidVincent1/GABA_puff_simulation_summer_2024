TITLE hh.mod   squid sodium, potassium, and leak channels

COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium,
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set a squid-appropriate temperature
 (e.g. in HOC: "celsius=6.3" or in Python: "h.celsius=6.3").
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

? interface
NEURON {
        SUFFIX hhrat
        REPRESENTS NCIT:C17145   : sodium channel
        REPRESENTS NCIT:C17008   : potassium channel
        USEION na READ ena WRITE ina REPRESENTS CHEBI:29101
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, Vt, Vs, ik, ina, il
        : `GLOBAL minf` will be replaced with `RANGE minf` if CoreNEURON enabled
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
        THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        gnabar = .12 (S/cm2)    <0,1e9> : .040
        gkbar = .036 (S/cm2)    <0,1e9> : .035
        gl = .0003 (S/cm2)      <0,1e9>
        el = -54.3 (mV)
        Vt = -58 (mV) :??
        Vs = -10 (mV) :??
}

STATE {
        m
        h
        n
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
        mtau (ms)
        htau (ms)
        ntau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)
        il = gl*(v - el)
}


INITIAL {
        rates(v)
        m = minf
        h = hinf
        n = ninf
}

? states
DERIVATIVE states {
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}

:LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum, q10
        : `TABLE minf` will be replaced with `:TABLE minf` if CoreNEURON enabled
        TABLE minf, mtau, hinf, htau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        :q10 = 3^((celsius - 6.3)/10)

        :        :"m" sodium activation system
        :alpha = .32 * vtrap(-(v-Vt-13),5)
        :beta =  0.28 * exp(-(v-Vt-40)/5)
        :sum = alpha + beta
        :mtau = 1/(q10*sum)
        :minf = alpha/sum

                :"h" sodium inactivation system
        :alpha = .128 * exp(-(v-Vt-Vs-17)/18)
        :beta = 4 / (exp(-(v-Vt-Vs-40)/5) + 1)
        :sum = alpha + beta
        :htau = 1/(q10*sum)
        :hinf = alpha/sum

                :"n" potassium activation system
        :alpha = .0032*vtrap(-(v-Vt-15),5)
        :beta = .5*exp(-(v-Vt-10)/40)
        :sum = alpha + beta
        :ntau = 1/(q10*sum)
        :ninf = alpha/sum

        q10 = 3^((celsius - 23)/10)

                :"m" sodium activation system
        alpha = -0.182 * vtrap(-(v+35), 9)
        beta =  -0.124 * vtrap((v+35),9)
        sum = alpha + beta
        mtau = 1/(q10*sum)
        minf = alpha/sum

                :"h" sodium inactivation system
        alpha = 0.25 * exp(-(v+90)/12)
        beta = 0.25 * exp((v+62)/6 - (v+90)/12)
        sum = alpha + beta
        htau = 1/(q10*sum)
        hinf = alpha/sum

                :"n" potassium activation system
        alpha = -0.02 * vtrap(-(v-25), 9)
        beta = -0.002 * vtrap((v-25),9)
        sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = -y + x/2
        }else{
                vtrap = x/(1-exp(x/y))
        }
}

UNITSON