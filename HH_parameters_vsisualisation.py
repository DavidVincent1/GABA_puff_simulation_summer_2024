import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'font.size': 22}) # Graph parameters


# Range of potential ----------------------
# We want m, n, h = 0.5 at -40, -20, -20 mV 
v = np.linspace(-100, 100, 10000)

def expcc(x, y):
    return x/(np.exp(x/y)-1)

def expcc2(x, y):
    return x/(1-np.exp(x/y))


# From Neuronal Dynamics ---------------------
# sodium
alpham2 = -0.182 * expcc2(-(v+35), 9)
betam2 = -0.124 * expcc2((v+35),9)
taum2 = 1/(alpham2 + betam2)
m2 = alpham2*taum2

alphah2 = 0.25 * np.exp(-(v+90)/12)
betah2 = 0.25 * np.exp((v+62)/6 - (v+90)/12)
tauh2 = 1/(alphah2 + betah2)
h2 = alphah2*tauh2

# potassium
alphan2 = -0.02 * expcc2(-(v-25), 9)
betan2 = -0.002 * expcc2((v-25),9)
taun2 = 1/(alphan2 + betan2)
n2 = alphan2*taun2


# Graphs --------------------------------------------------------------
fig, ax = plt.subplots(1, 2)
fig.canvas.manager.set_window_title('HH_channels_parameters')

ax[0].plot(v, m2, color='black', linestyle='-', label='m')
ax[0].plot(v, n2, color='black', linestyle='--', label='n')
ax[0].plot(v, h2, color='black', linestyle=':', label='h')

ax[1].plot(v, taum2, color='black', linestyle='-', label=r'$\tau_m$')
ax[1].plot(v, taun2, color='black', linestyle='--', label=r'$\tau_n$')
ax[1].plot(v, tauh2, color='black', linestyle=':', label=r'$\tau_h$')

ax[0].set_ylabel(r"$x_{0}$ [-]")
ax[1].set_ylabel(r"$\tau_{x}$ [ms]")
ax[0].set_xlabel("Potential [mV]")
ax[1].set_xlabel("Potential [mV]")

ax[0].legend()
ax[1].legend()


plt.show()