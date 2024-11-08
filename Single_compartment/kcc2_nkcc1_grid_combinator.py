import pickle
import matplotlib.pyplot as plt
import numpy as np
from neuron import h
from matplotlib import rcParams
rcParams.update({'font.size': 22}) # Graph parameters


# g_clc2 value for each graph
g_clc2_1 = 0
g_clc2_2 = 1e-5
g_clc2_3 = 3e-5
g_clc2_4 = 1e-4


# Pickle files ------------------------------------------------------------------------------------------------------------SS
with open('dataset\pickle_kcc2=1e-06-0.0005_nkcc1=1e-06-0.0005_gclc2=0.0_clamp=-90mV\Stable_chloride.pickle', 'rb') as f:
    fig1 = pickle.load(f)
with open('dataset\pickle_kcc2=1e-06-0.0005_nkcc1=1e-06-0.0005_gclc2=1e-06_clamp=-90mV\Stable_chloride.pickle', 'rb') as f:
    fig2 = pickle.load(f)
with open('dataset\pickle_kcc2=1e-06-0.0005_nkcc1=1e-06-0.0005_gclc2=5e-06_clamp=-90mV\Stable_chloride.pickle', 'rb') as f:
    fig3 = pickle.load(f)
with open('dataset\pickle_kcc2=1e-06-0.0005_nkcc1=1e-06-0.0005_gclc2=1e-05_clamp=-90mV\Stable_chloride.pickle', 'rb') as f:
    fig4 = pickle.load(f)


# Obtenir l'axe de la grille (image) pour chaque graphique
ax1 = fig1.axes[0]
ax2 = fig2.axes[0]
ax3 = fig3.axes[0]
ax4 = fig4.axes[0]

# Obtenir les objets image (cela peut être un `AxesImage` ou similaire)
im1 = ax1.get_images()[0]
im2 = ax2.get_images()[0]
im3 = ax3.get_images()[0]
im4 = ax4.get_images()[0]

# Récupérer les données utilisées dans chaque image
data1 = im1.get_array()
data2 = im2.get_array()
data3 = im3.get_array()
data4 = im4.get_array()

# Trouver les valeurs minimales et maximales parmi toutes les grilles
vmin = min(data1.min(), data2.min(), data3.min(), data4.min())
vmax = max(data1.max(), data2.max(), data3.max(), data4.max())


Unkcc1 = np.linspace(1e-6, 5e-4, num=len(data1[0]))
Ukcc2 = np.linspace(5e-4, 1e-6, num=len(data1[0]))
Unkcc1_label = [f'{val:.2e}' for val in Unkcc1]
Ukcc2_label = [f'{val:.2e}' for val in Ukcc2]


# Replication of the graphs considering the new min and max values for the colorbar -------------------
# First graph
fig_1 = plt.figure("Stable chloride - Updated")
im_1 = plt.imshow(data1, cmap='cividis', vmin=vmin, vmax=vmax)  # vmin and vmax
cbar_1 = plt.colorbar(im_1, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')

for i in range(len(data1)):
    for j in range(len(data1[i])):
        text = plt.text(j, i, f'{data1[i, j]:.3f}', ha="center", va="center", color="w")

plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
plt.ylabel(r"$U_{kcc2}$ [mM/ms]")


# Second graph
fig_2 = plt.figure("Stable chloride - Updated 2")
im_2 = plt.imshow(data2, cmap='cividis', vmin=vmin, vmax=vmax)  # vmin and vmax
cbar_2 = plt.colorbar(im_2, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')

for i in range(len(data2)):
    for j in range(len(data2[i])):
        text = plt.text(j, i, f'{data2[i, j]:.3f}', ha="center", va="center", color="w")

plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
plt.ylabel(r"$U_{kcc2}$ [mM/ms]")


# Third graph
fig_3 = plt.figure("Stable chloride - Updated 3")
im_3 = plt.imshow(data3, cmap='cividis', vmin=vmin, vmax=vmax)  # vmin and vmax
cbar_3 = plt.colorbar(im_3, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')

for i in range(len(data3)):
    for j in range(len(data3[i])):
        text = plt.text(j, i, f'{data3[i, j]:.3f}', ha="center", va="center", color="w")

plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
plt.ylabel(r"$U_{kcc2}$ [mM/ms]")


# Fourth graph
fig_4 = plt.figure("Stable chloride - Updated 4")
im_4 = plt.imshow(data4, cmap='cividis', vmin=vmin, vmax=vmax)  # vmin and vmax
cbar_4 = plt.colorbar(im_4, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')

for i in range(len(data4)):
    for j in range(len(data4[i])):
        text = plt.text(j, i, f'{data4[i, j]:.3f}', ha="center", va="center", color="w")

plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
plt.ylabel(r"$U_{kcc2}$ [mM/ms]")



# Combination of the four graphs into a four panels figure ---------------------------------------
fig, axes = plt.subplots(2, 2)

# First graph
im1 = axes[0, 0].imshow(data1, cmap='cividis', vmin=vmin, vmax=vmax)
axes[0, 0].set_title(r'$g_{clc2} = $' + f'{g_clc2_1} ' + r'$S/\rm{cm}^{2}$', fontsize=19)
for i in range(len(data1)):
    for j in range(len(data1[i])):
        text = axes[0, 0].text(j, i, f'{data1[i, j]:.1f}', ha="center", va="center", color="w", fontsize=14)

# Second graph
im2 = axes[0, 1].imshow(data2, cmap='cividis', vmin=vmin, vmax=vmax)
axes[0, 1].set_title(r'$g_{clc2} = $' + f'{g_clc2_2} ' + r'$S/\rm{cm}^{2}$', fontsize=19)
for i in range(len(data2)):
    for j in range(len(data2[i])):
        text = axes[0, 1].text(j, i, f'{data2[i, j]:.1f}', ha="center", va="center", color="w", fontsize=14)

# Third graph
im3 = axes[1, 0].imshow(data3, cmap='cividis', vmin=vmin, vmax=vmax)
axes[1, 0].set_title(r'$g_{clc2} = $' + f'{g_clc2_3} ' + r'$S/\rm{cm}^{2}$', fontsize=19)
for i in range(len(data3)):
    for j in range(len(data3[i])):
        text = axes[1, 0].text(j, i, f'{data3[i, j]:.1f}', ha="center", va="center", color="w", fontsize=14)

# Fourth graph
im4 = axes[1, 1].imshow(data4, cmap='cividis', vmin=vmin, vmax=vmax)
axes[1, 1].set_title(r'$g_{clc2} = $' + f'{g_clc2_4} ' + r'$S/\rm{cm}^{2}$', fontsize=19)
for i in range(len(data4)):
    for j in range(len(data4[i])):
        text = axes[1, 1].text(j, i, f'{data4[i, j]:.1f}', ha="center", va="center", color="w", fontsize=14)

for i in range(2):
    axes[i,0].set_ylabel(r"$U_{kcc2}$ [mM/ms]", fontsize=18)
    axes[i,0].set_yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label, fontsize='16')
    axes[1,i].set_xlabel(r"$U_{nkcc1}$ [mM/ms]", fontsize=18)
    axes[1,i].set_xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor", fontsize='16')
    axes[i,1].set_yticks(np.arange(len(Ukcc2)))
    axes[0,i].set_xticks(np.arange(len(Unkcc1)))
    axes[i,1].set_yticklabels([])
    axes[0,i].set_xticklabels([])

# Colorbar for all graphs
#cbar = fig.colorbar(im1, ax=axes, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]', orientation='vertical')
cbar = fig.colorbar(im1, ax=axes, extend='both', spacing='proportional', label=r'$MP$ before puff [mV]', orientation='vertical')
#cbar = fig.colorbar(im1, ax=axes, extend='both', spacing='proportional', label=r'$E_{cl}$ before puff [mV]', orientation='vertical')

plt.show()