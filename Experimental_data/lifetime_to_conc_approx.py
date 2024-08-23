import numpy as np
import pandas as pd
import matplotlib.pyplot as plt, cmap
from matplotlib import rcParams
from scipy.ndimage import gaussian_filter1d
import h5py


# For graphs --------------------------------------------
colormap = cmap.Colormap("viridis")(np.linspace(0, 1, 9))
rcParams.update({'font.size': 22})


# Experimental data --------------------------------------------------------------------------
# INPUT 1
# the data is a excel file
path = "Experimental_data\Autres\lifetimes.xlsx"
indice = pd.read_excel(path, sheet_name=0)['ROI ID']
data = pd.read_excel(path, sheet_name=1)
data = np.asarray(data) # excel -> np.array


# Simulation results to compare with experimental data ----------------------------------------------------------------------------------------------------------------------------------
# INPUT 2
#f = h5py.file(r"Path of the dataset file", 'r')
f = h5py.File(r"Path of the dataset file", 'r')

decal = 470000                # Offset to skip the initial stabilization of the simulation
t_and_soma_v = f["mp_soma"]   # Time and MP in soma
t_ = t_and_soma_v[0] - decal  # Shifted time
t = t_/1000                   # Time array in second
dend_conc = f["conc_dend"]    # Concentrations in dendrite (at recording pos)
dend_cli = dend_conc[0]       # Chloride concentration in dendrite

# Other values used later
record_pos = dend_conc.attrs["recording_locations"]
syn_pos = dend_conc.attrs["Synapses positions"]
center_ind = int(len(syn_pos)/2-1)
lenght = dend_conc.attrs["dendrite_full_length_um"][0]
puff_time = dend_conc.attrs["puff_time"]
puff_pos = dend_conc.attrs["puff_position"]
sim_time = dend_conc.attrs["simulation_length"]
rnum = dend_conc.attrs["number_of_receptors"]


# Separation of the data into time_array and lifetime array
time_array = data[:,0]
lifetime_arrays = [data[:,i+1] for i in range(9)]


# Filtering the data and calculating the approximate chloride concentration -----------------------------------------------
shift = 0 # Shift to make the chloride concentration be positive
lifetime_arrays_filtered = [gaussian_filter1d(lifetime_arrays[i], 5) for i in range(len(lifetime_arrays))]
k_q = 0.0044587465
un_sur_tau0 = 0.1897855695
conc_approx = [(1/k_q)*(1/lifetime_arrays_filtered[i] - un_sur_tau0) + shift for i in range(len(lifetime_arrays_filtered))]


# Positions of the recording vectors that fit with the experimental data
# Away from soma
syn_pos_imp = [(55 + (45*n)/9)/lenght for n in range(9)]

# Toward soma
syn_pos_imp2 = [(55 - (45*n)/9)/lenght for n in range(9)]


# Finding the recording vectors that are the nearest from the positions
index = []
for j in syn_pos_imp:
    val = 0
    ind = 0
    for i,k in enumerate(record_pos):
        if abs(j-k) < abs(j-val):
            val = k
            ind = i
    index.append(ind)
index.reverse

index2 = []
for j in syn_pos_imp2:
    val = 0
    ind = 0
    for i,k in enumerate(record_pos):
        if abs(j-k) < abs(j-val):
            val = k
            ind = i
    index2.append(ind)
index2.reverse

pos = [55+record_pos[i]*lenght for i in index]
pos2 = [55-record_pos[i]*lenght for i in index2]


# Label for the legend ----------------------------------------------------------------------------------------------------
seg = ["segment 1", "segment 2", "segment 3", "segment 4", "segment 5", "segment 6", "segment 7", "segment 8", "segment 9"]
label = ['~ 2.5 µm', '~ 7.2 µm', '~ 12 µm', '~ 16.8 µm', '~ 21.6 µm', '~ 26.4 µm', '~ 31.2 µm', '~ 36 µm', '~ 40.8 µm']
label.reverse


# Normalization of the experimental data by the same min/max ------------------------------------
maxi = np.max([np.nanmax(conc_approx), np.nanmax(dend_cli)])
mini = np.min([np.nanmin(conc_approx), np.nanmin(dend_cli)])
normalized_dend_cli = [(dend_cli[i] - mini)/(maxi-mini) for i in index2]
normalized_conc_approx = [(conc_approx[i] - mini)/(maxi-mini) for i in range(len(conc_approx))]


# Graphs --------------------------------------------------------------------------------------------------------------
# Comparison between experimental and simulation with simulation results away from soma
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title("chloride_comparison_simulation-exp_away_from_soma")
for i, val in enumerate(conc_approx): 
    ax[0].plot(time_array, val, color=colormap[i], label=label[i])
for i, val in enumerate(indice):
    ax[1].plot(t, dend_cli[val], color=colormap[-(i+1)])
plt.xlim(0, (sim_time-decal)/1000)
ax[1].set_xlabel("Time [s]")
ax[0].set_ylabel(r"$[Cl^-]_i$ [mM]")
ax[1].set_ylabel(r"$[Cl^-]_i$ [mM]")


# Comparison between experimental and simulation with simulation results toward soma
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title("chloride_comparison_simulation-exp_towards_soma")
for i, val in enumerate(conc_approx): 
    ax[0].plot(time_array, val, color=colormap[i])
for i, val in enumerate(index2):
    ax[1].plot(t, dend_cli[val], color=colormap[i])
ax[0].set_xlim(0, (sim_time-decal)/1000)
ax[1].set_xlim(0, (sim_time-decal)/1000)
ax[1].set_xlabel("Time [s]")
ax[0].set_ylabel(r"$[Cl^-]_i$ [mM]")
ax[1].set_ylabel(r"$[Cl^-]_i$ [mM]")


# Relative comparison between experimental and simulation toward soma
fig, ax = plt.subplots(2,1)
fig.canvas.manager.set_window_title("Relative_chloride_comparison_simulation-exp_towards_soma")
for i, val in enumerate(normalized_conc_approx): 
    ax[0].plot(time_array, val, color=colormap[i])
for i, val in enumerate(normalized_dend_cli):
    ax[1].plot(t, val, color=colormap[i])
ax[0].set_xlim(0, (sim_time-decal)/1000)
ax[1].set_xlim(0, (sim_time-decal)/1000)
ax[0].set_ylim(0, 1)
ax[1].set_ylim(0, 1)
ax[1].set_xlabel("Time [s]")
ax[0].set_ylabel('Normalized chloride\n' + "concentration [a.u.]")
ax[1].set_ylabel('Normalized chloride\n' + "concentration [a.u.]")

# Figure for the label of the last two graphs
plt.figure('Label_comparison_graph')
for i, val in enumerate(conc_approx): 
    plt.plot(time_array, np.zeros_like(time_array), color=colormap[-(i+1)], label=label[i])
plt.ylim(1,2)
plt.legend(loc='center')


# Figure showing the lifetime arrays
plt.figure("Lifetime chloride")
for i,j in enumerate(lifetime_arrays_filtered): 
    plt.plot(time_array, j, color=colormap[i], label=seg[i])
plt.xlabel("Time [s]")
plt.ylabel("Lifetime [ns]")
plt.legend(fontsize=18)


# Figure showing the chloride concentration approximation
plt.figure("chloride concentration approx")
for i,j in enumerate(conc_approx): 
    plt.plot(time_array, j, color=colormap[i], label=seg[i])
plt.xlabel("Time [s]")
plt.ylabel("Chloride concentration [mM]")
plt.legend(fontsize=18)



plt.show()