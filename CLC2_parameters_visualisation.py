import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate import solve_ivp
rcParams.update({'font.size': 22}) # Graph parameters


# Range of potential and Ecl -------------
v = np.linspace(-100, 100, 1000) # mV
ecl = [-40, -50, -60, -70, -80]  # mV
gclc2 = 1e-6                     # S/cm^2


# function of the opening variable p of ClC-2 -----------
def pfunc(t, y, x=-70, ecl=-60, v1=15, v2=-14, taup=300):
    return (1/(1 + np.exp((ecl - x - v1)/v2)) - y)/taup


# Calculate Iclc2 for different value of v and for one value of Ecl -------------------
for Ecl in ecl:
    p = []
    iclc2 = []
    for i, val in enumerate(v):
        # Integration
        result = solve_ivp(pfunc, (0, 20000), y0=[0.5], args=(val, Ecl, 15, -14, 300))

        # mean value of p for the last 10 evaluation points
        p_ = np.mean(result.y[0][-11:-1])
        p.append(p_)

        # ClC-2 current
        iclc2.append(gclc2*p_*(val - Ecl))

    # Graphs
    plt.figure("Courbe I-V")
    plt.plot(v, iclc2, label=r'$E_{cl} = $' + f'{Ecl}')

    plt.figure("Courbe p-V")
    plt.plot(v, p, label=r'$E_{cl} = $' + f'{Ecl}')

plt.figure("Courbe I-V")
plt.ylabel('I [mA/cm2]')
plt.xlabel('v [mV]')
plt.legend()

plt.figure("Courbe p-V")
plt.ylabel('p [-]')
plt.xlabel('v [mV]')
plt.legend()

plt.figure('Last integration')
plt.plot(result.t, result.y[0])
plt.xlabel('t [ms]')
plt.ylabel('p [-]')

plt.show()