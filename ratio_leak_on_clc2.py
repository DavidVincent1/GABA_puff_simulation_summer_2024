import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rcParams


# Graph parameters -----------------------------
rcParams.update({'font.size': 22})
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'


ratio_clc2_on_leak = [0, 1/3, 1, 3]
ratio = np.linspace(0, 3, 10000)

# high NKCC1, unclamp
cl1 = [27.86719739270441, 27.9378788425, 28.007839168976663, 28.076992748780828]

# high NKCC1, clamped at -90 mV
cl2 = [21.317966721450592, 21.34835782913345, 21.39857611255515, 21.448618058458678]

# high KCC2, unclamp
cl3 = [4.906641152164216, 4.884998256339075, 4.860319580646839, 4.833785803467128]

# high KCC2, clamped at -50 mV
cl4 = [6.564003676848365, 6.445300698017499, 6.323576321577707, 6.198709974914942]

# To fit
def droite(x, a, b):
    return a*x + b

popt1, pcov1 = curve_fit(droite, ratio_clc2_on_leak, cl1)
popt2, pcov2 = curve_fit(droite, ratio_clc2_on_leak, cl2)
popt3, pcov3 = curve_fit(droite, ratio_clc2_on_leak, cl3)
popt4, pcov4 = curve_fit(droite, ratio_clc2_on_leak, cl4)

fit1 = droite(ratio, *popt1)
fit2 = droite(ratio, *popt2)
fit3 = droite(ratio, *popt3)
fit4 = droite(ratio, *popt4)


fig, ax = plt.subplots(2,2)

ax[0,0].set_title('NKCC1 Élevé, Non clampé', fontsize=16)
ax[0,0].scatter(ratio_clc2_on_leak, cl1, color='blue', label=r'$[Cl^-]_i$')
ax[0,0].plot(ratio, fit1, color='blue', label='Ajustement de courbe\n' + f'y = {popt1[0]:.3f}x + {popt1[1]:.3f}', linestyle='--')
ax[0,0].set_xlabel(r"$g_{clc2}/g_{cl,leak}$ (-)")
ax[0,0].set_ylabel(r"$[Cl^-]_i$ (mM)")
ax[0,0].legend(fontsize=14, loc='lower right')

ax[0,1].set_title('High NKCC1, Clampé à -90 mV', fontsize=16)
ax[0,1].scatter(ratio_clc2_on_leak, cl2, color='blue', label=r'$[Cl^-]_i$')
ax[0,1].plot(ratio, fit2, color='blue', label='Ajustement de courbe\n' + f'y = {popt2[0]:.3f}x + {popt2[1]:.3f}', linestyle='--')
ax[0,1].set_xlabel(r"$g_{clc2}/g_{cl,leak}$ (-)")
ax[0,1].set_ylabel(r"$[Cl^-]_i$ (mM)")
ax[0,1].legend(fontsize=14, loc='lower right')

ax[1,0].set_title('KCC2 Élevé, Non Clampé', fontsize=16)
ax[1,0].scatter(ratio_clc2_on_leak, cl3, color='blue', label=r'$[Cl^-]_i$')
ax[1,0].plot(ratio, fit3, color='blue', label='Ajustement de courbe\n' + f'y = {popt3[0]:.3f}x + {popt3[1]:.3f}', linestyle='--')
ax[1,0].set_xlabel(r"$g_{clc2}/g_{cl,leak}$ (-)")
ax[1,0].set_ylabel(r"$[Cl^-]_i$ (mM)")
ax[1,0].legend(fontsize=14)

ax[1,1].set_title('KCC2 Élevé, Clampé à -50 mV', fontsize=16)
ax[1,1].scatter(ratio_clc2_on_leak, cl4, color='blue', label=r'$[Cl^-]_i$')
ax[1,1].plot(ratio, fit4, color='blue', label='Ajustement de courbe\n' + f'y = {popt4[0]:.3f}x + {popt4[1]:.3f}', linestyle='--')
ax[1,1].set_xlabel(r"$g_{clc2}/g_{cl,leak}$ (-)")
ax[1,1].set_ylabel(r"$[Cl^-]_i$ (mM)")
ax[1,1].legend(fontsize=14)


plt.show()