import numpy as np
import matplotlib.pyplot as plt, cmap
from matplotlib import rcParams
import h5py
from matplotlib.gridspec import GridSpec


# Graph parameters -----------------------------
rcParams.update({'font.size': 14})
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'


# Files path ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
path1 = r'dataset\adult_leak=0.2_clc2_varie\-50mV\clamped=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffcon=1_Dcl=2_gclc2=1e-05.hdf5'
path2 = r'dataset\fork\clamp=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffco=1.00_Dcl=2.0_gclc2=1e-05_fork=0.80.hdf5'
path3 = r'dataset\fork\clamp=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffco=1.00_Dcl=2.0_gclc2=1e-05_fork=0.60.hdf5'
path4 = r'dataset\fork\clamp=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffco=1.00_Dcl=2.0_gclc2=1e-05_fork=0.47.hdf5'

title_label = ["No fork", "Fork at 120 µm", "Fork at 90 µm", "Fork at 70 µm"]
title_label_fr = ["", "", "", ""]

cl_conc = False
clc2_leak_currents = True

graph_fr = False
decal = 590000


# First file --------------------------------------------------------------------------------------------
f1 = h5py.File(path1, 'r')
t_and_soma_v1 = f1["mp_soma"] # Time and MP in soma
t_1 = t_and_soma_v1[0]         # Time array [ms]
t1 = (t_1 - decal)/1000        # Time array with shifted zero [s]
dend_conc1 = f1["conc_dend"]
soma_conc1 = f1["conc_soma"]
soma_current_cl1 = f1["soma_current_cl"] # Icl in soma
soma_cli1 = soma_conc1[0]
dend_cli1 = dend_conc1[0]
icl_leak1 = soma_current_cl1[3]
icl_clc21 = soma_current_cl1[4]

record_pos = dend_conc1.attrs["recording_locations"]    # Recording positions for the concentrations [-]
lenght = dend_conc1.attrs["dendrite_full_length_um"][0] # Dendrite (first part) lenght [µm]
puff_time = dend_conc1.attrs["puff_time"]               # Time at wich the GABA puff occurs [ms]
puff_pos = dend_conc1.attrs["puff_position"]            # Position of the GABA puff [-]
sim_time = dend_conc1.attrs["simulation_length"]        # Simulation lenght [ms]
center = 100
for j in record_pos:
    if abs(puff_pos-j) < abs(puff_pos-center):
        center = j

dend_cli_begin1, dend_cli_end1 = [], []
record_pos_begin, record_pos_end = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin1.append(dend_cli1[i])
        record_pos_begin.append(record_pos[i])
    if record_pos[i] >= center:
        dend_cli_end1.append(dend_cli1[i])
        record_pos_end.append(record_pos[i])
print(r'Number of recording verctors toward the soma : ', len(dend_cli_begin1))
print(r'Number of recording verctors away from soma  : ', len(dend_cli_end1))

colormap = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(record_pos)/2)))
colormap_begin = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(dend_cli_begin1))))
colormap_end = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(dend_cli_end1))))
xlimit_min = 0
xlimit_max = (sim_time-decal)/1000

del f1, t_and_soma_v1, t_1, dend_conc1, dend_cli1, soma_conc1, soma_current_cl1


# Second file ------------------------------------------------------------------
f2 = h5py.File(path2, 'r')
t_and_soma_v2 = f2["mp_soma"] # Time and MP in soma
t_2 = t_and_soma_v2[0]         # Time array [ms]
t2 = (t_2 - decal)/1000        # Time array with shifted zero [s]
dend_conc2 = f2["conc_dend"]
soma_conc2 = f2["conc_soma"]
soma_current_cl2 = f2["soma_current_cl"] # Icl in soma
icl_leak2 = soma_current_cl2[3]
icl_clc22 = soma_current_cl2[4]
soma_cli2 = soma_conc2[0]
dend_cli2 = dend_conc2[0]


dend_cli_begin2, dend_cli_end2 = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin2.append(dend_cli2[i])
    if record_pos[i] >= center:
        dend_cli_end2.append(dend_cli2[i])
print(r'Number of recording verctors toward the soma : ', len(dend_cli_begin2))
print(r'Number of recording verctors away from soma  : ', len(dend_cli_end2))

del f2, t_and_soma_v2, t_2, dend_conc2, dend_cli2, soma_conc2


# Third file ---------------------------------------------------------------------
f3 = h5py.File(path3, 'r')
t_and_soma_v3 = f3["mp_soma"] # Time and MP in soma
t_3 = t_and_soma_v3[0]         # Time array [ms]
t3 = (t_3 - decal)/1000        # Time array with shifted zero [s]
dend_conc3 = f3["conc_dend"]
soma_conc3 = f3["conc_soma"]
soma_cli3 = soma_conc3[0]
dend_cli3 = dend_conc3[0]
soma_current_cl3 = f3["soma_current_cl"] # Icl in soma
icl_leak3 = soma_current_cl3[3]
icl_clc23 = soma_current_cl3[4]

dend_cli_begin3, dend_cli_end3 = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin3.append(dend_cli3[i])
    if record_pos[i] >= center:
        dend_cli_end3.append(dend_cli3[i])
print(r'Number of recording vectors toward the soma : ', len(dend_cli_begin3))
print(r'Number of recording vectors away from soma  : ', len(dend_cli_end3))

del f3, t_and_soma_v3, t_3, dend_conc3, dend_cli3, soma_conc3


# Fourth file ------------------------------------------------------------------
f4 = h5py.File(path4, 'r')
t_and_soma_v4 = f4["mp_soma"] # Time and MP in soma
t_4 = t_and_soma_v4[0]         # Time array [ms]
t4 = (t_4 - decal)/1000        # Time array with shifted zero [s]
dend_conc4 = f4["conc_dend"]
soma_conc4 = f4["conc_soma"]
soma_cli4 = soma_conc4[0]
dend_cli4 = dend_conc4[0]
soma_current_cl4 = f4["soma_current_cl"] # Icl in soma
icl_leak4 = soma_current_cl4[3]
icl_clc24 = soma_current_cl4[4]

dend_cli_begin4, dend_cli_end4 = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin4.append(dend_cli4[i])
    if record_pos[i] >= center:
        dend_cli_end4.append(dend_cli4[i])
print(r'Number of recording vectors toward the soma : ', len(dend_cli_begin4))
print(r'Number of recording vectors away from soma  : ', len(dend_cli_end4))

del f4, t_and_soma_v4, t_4, dend_conc4, dend_cli4, soma_conc4


# Maximum and minimum values for y-axis ----------------------
max = np.max([np.max(dend_cli_begin1), np.max(dend_cli_end1),
            np.max(dend_cli_begin2), np.max(dend_cli_end2),
            np.max(dend_cli_begin3), np.max(dend_cli_end3),
            np.max(dend_cli_begin4), np.max(dend_cli_end4)])
min = np.min([np.min(dend_cli_begin1), np.min(dend_cli_end1),
            np.min(dend_cli_begin2), np.min(dend_cli_end2),
            np.min(dend_cli_begin3), np.min(dend_cli_end3),
            np.min(dend_cli_begin4), np.min(dend_cli_end4)])


# Figure time ----------------------------------------------------------------------
# Figure creation
fig = plt.figure(figsize=(19, 11))

# To add an empty space between the central rows
gs = GridSpec(5, 2, figure=fig, height_ratios=[1, 1, 0.1, 1, 1])
ax = [[0,0],[0,0],[0,0],[0,0]]
for i in range(4):
    for j in range(2):
        if i < 2:
            ax_ = fig.add_subplot(gs[i, j])
            ax[i][j] = ax_
        elif i >= 2:
            ax_ = fig.add_subplot(gs[i+1, j])
            ax[i][j] = ax_

# y-axis max and min
max_cl = max
min_cl = min

# First graph
if graph_fr:
    ax[0][0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[0][0].set_title(title_label[0], fontsize=15)

for i, val in enumerate(dend_cli_begin1):
    if i == 0:
        if graph_fr:
            ax[0][0].plot(t1, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[0][0].plot(t1, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin1)-1:
        if graph_fr:
            ax[0][0].plot(t1, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[0][0].plot(t1, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[0][0].plot(t1, val, color=colormap_begin[i])

if graph_fr:
    ax[0][0].plot(t1, soma_cli1, color='red', label=f"dans le soma", zorder=3)
else:
    ax[0][0].plot(t1, soma_cli1, color='red', label=f"in soma", zorder=3)

for i, val in enumerate(dend_cli_end1):
    if i == 0:
        if graph_fr:
            ax[1][0].plot(t1, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[1][0].plot(t1, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end1)-1:
        if graph_fr:
            ax[1][0].plot(t1, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[1][0].plot(t1, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[1][0].plot(t1, val, color=colormap_end[-(i+1)])

ax[0][0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[0][0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[0][0].set_xlim(xlimit_min, xlimit_max)

ax[1][0].set_xlim(xlimit_min, xlimit_max)
ax[1][0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[1][0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[1][0].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[1][0].set_xlabel("Time (s)", fontsize=13)

ax[0][0].legend(loc='upper right', fontsize=14)
ax[1][0].legend(loc='upper right', fontsize=14)

del dend_cli_begin1
del dend_cli_end1
del soma_cli1



# Second graph
if graph_fr:
    ax[0][1].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[0][1].set_title(title_label[1], fontsize=15)

for i, val in enumerate(dend_cli_begin2):
    if i == 0:
        if graph_fr:
            ax[0][1].plot(t2, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[0][1].plot(t2, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin2)-1:
        if graph_fr:
            ax[0][1].plot(t2, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[0][1].plot(t2, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[0][1].plot(t2, val, color=colormap_begin[i])

if graph_fr:
    ax[0][1].plot(t2, soma_cli2, color='red', label=f"dans le soma", zorder=3)
else:
    ax[0][1].plot(t2, soma_cli2, color='red', label=f"in soma", zorder=3)

for i, val in enumerate(dend_cli_end2):
    if i == 0:
        if graph_fr:
            ax[1][1].plot(t2, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[1][1].plot(t2, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end2)-1:
        if graph_fr:
            ax[1][1].plot(t2, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[1][1].plot(t2, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[1][1].plot(t2, val, color=colormap_end[-(i+1)])

ax[0][1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[0][1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[0][1].set_xlim(xlimit_min, xlimit_max)

ax[1][1].set_xlim(xlimit_min, xlimit_max)
ax[1][1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[1][1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[1][1].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[1][1].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin2
del dend_cli_end2
del soma_cli2



# Third graph
if graph_fr:
    ax[2][0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[2][0].set_title(title_label[2], fontsize=15)

for i, val in enumerate(dend_cli_begin3):
    if i == 0:
        if graph_fr:
            ax[2][0].plot(t3, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[2][0].plot(t3, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin3)-1:
        if graph_fr:
            ax[2][0].plot(t3, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[2][0].plot(t3, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[2][0].plot(t3, val, color=colormap_begin[i])

if graph_fr:
    ax[2][0].plot(t3, soma_cli3, color='red', label=f"dans le soma", zorder=3)
else:
    ax[2][0].plot(t3, soma_cli3, color='red', label=f"in soma", zorder=3)

for i, val in enumerate(dend_cli_end3):
    if i == 0:
        if graph_fr:
            ax[3][0].plot(t3, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[3][0].plot(t3, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end3)-1:
        if graph_fr:
            ax[3][0].plot(t3, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[3][0].plot(t3, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[3][0].plot(t3, val, color=colormap_end[-(i+1)])

ax[2][0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[2][0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[2][0].set_xlim(xlimit_min, xlimit_max)

ax[3][0].set_xlim(xlimit_min, xlimit_max)
ax[3][0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[3][0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[3][0].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[3][0].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin3
del dend_cli_end3
del soma_cli3



# Fourth graph
if graph_fr:
    ax[2][1].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[2][1].set_title(title_label[3], fontsize=15)

for i, val in enumerate(dend_cli_begin4):
    if i == 0:
        if graph_fr:
            ax[2][1].plot(t4, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[2][1].plot(t4, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin4)-1:
        if graph_fr:
            ax[2][1].plot(t4, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[2][1].plot(t4, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[2][1].plot(t4, val, color=colormap_begin[i])

if graph_fr:
    ax[2][1].plot(t4, soma_cli4, color='red', label=f"dans le soma", zorder=3)
else:
    ax[2][1].plot(t4, soma_cli4, color='red', label=f"in soma", zorder=3)

for i, val in enumerate(dend_cli_end4):
    if i == 0:
        if graph_fr:
            ax[3][1].plot(t4, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[3][1].plot(t4, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end4)-1:
        if graph_fr:
            ax[3][1].plot(t4, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[3][1].plot(t4, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[3][1].plot(t4, val, color=colormap_end[-(i+1)])

ax[2][1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[2][1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[2][1].set_xlim(xlimit_min, xlimit_max)

ax[3][1].set_xlim(xlimit_min, xlimit_max)
ax[3][1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[3][1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[3][1].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[3][1].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin4
del dend_cli_end4
del soma_cli4

plt.tight_layout()
plt.subplots_adjust(hspace=0.45, wspace=0.175)


plt.show()