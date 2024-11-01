import numpy as np
import matplotlib.pyplot as plt, cmap
from matplotlib import rcParams
import h5py


# Graph parameters -----------------------------
rcParams.update({'font.size': 18})
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'

path1 = r'dataset\adult_leak=0.2_clc2_varie\-50mV\clamped=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffcon=1_Dcl=2_gclc2=5e-06.hdf5'
path2 = r'dataset\adult_leak=0.2_clc2_varie\-50mV\clamped=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffcon=1_Dcl=2_gclc2=1e-05.hdf5'
path3 = r'dataset\adult_leak=0.2_clc2_varie\-50mV\clamped=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffcon=1_Dcl=2_gclc2=2e-05.hdf5'
path4 = r'dataset\adult_leak=0.2_clc2_varie\-50mV\clamped=-50mV_synnb=20_simlen=650000_dt=(5,0.144,0.144)_L=150_kcc2=0.0001_nkcc1=1e-05_rnum=200_puffcon=1_Dcl=2_gclc2=5e-05.hdf5'

graph_fr = False
decal = 590000
max = 9.1
min = 6

f1 = h5py.File(path1, 'r')
t_and_soma_v1 = f1["mp_soma"] # Time and MP in soma
t_1 = t_and_soma_v1[0]         # Time array [ms]
t1 = (t_1 - decal)/1000        # Time array with shifted zero [s]
dend_conc1 = f1["conc_dend"]
soma_conc1 = f1["conc_soma"]
soma_cli1 = soma_conc1[0]
dend_cli1 = dend_conc1[0]

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

del f1, t_and_soma_v1, t_1, dend_conc1, dend_cli1, soma_conc1


f2 = h5py.File(path2, 'r')
t_and_soma_v2 = f2["mp_soma"] # Time and MP in soma
t_2 = t_and_soma_v2[0]         # Time array [ms]
t2 = (t_2 - decal)/1000        # Time array with shifted zero [s]
dend_conc2 = f2["conc_dend"]
soma_conc2 = f2["conc_soma"]
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


f3 = h5py.File(path3, 'r')
t_and_soma_v3 = f3["mp_soma"] # Time and MP in soma
t_3 = t_and_soma_v3[0]         # Time array [ms]
t3 = (t_3 - decal)/1000        # Time array with shifted zero [s]
dend_conc3 = f3["conc_dend"]
soma_conc3 = f3["conc_soma"]
soma_cli3 = soma_conc3[0]
dend_cli3 = dend_conc3[0]

dend_cli_begin3, dend_cli_end3 = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin3.append(dend_cli3[i])
    if record_pos[i] >= center:
        dend_cli_end3.append(dend_cli3[i])
print(r'Number of recording verctors toward the soma : ', len(dend_cli_begin3))
print(r'Number of recording verctors away from soma  : ', len(dend_cli_end3))

del f3, t_and_soma_v3, t_3, dend_conc3, dend_cli3, soma_conc3


f4 = h5py.File(path4, 'r')
t_and_soma_v4 = f4["mp_soma"] # Time and MP in soma
t_4 = t_and_soma_v4[0]         # Time array [ms]
t4 = (t_4 - decal)/1000        # Time array with shifted zero [s]
dend_conc4 = f4["conc_dend"]
soma_conc4 = f4["conc_soma"]
soma_cli4 = soma_conc4[0]
dend_cli4 = dend_conc4[0]

dend_cli_begin4, dend_cli_end4 = [], []
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin4.append(dend_cli4[i])
    if record_pos[i] >= center:
        dend_cli_end4.append(dend_cli4[i])
print(r'Number of recording verctors toward the soma : ', len(dend_cli_begin4))
print(r'Number of recording verctors away from soma  : ', len(dend_cli_end4))

del f4, t_and_soma_v4, t_4, dend_conc4, dend_cli4, soma_conc4


# Figure time ----------------------------------------------------------------------
fig, ax = plt.subplots(4,2)

max_cl = max
min_cl = min

if graph_fr:
    ax[0,0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[0,0].set_title(f"Chloride concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_begin1):
    if i == 0:
        if graph_fr:
            ax[0,0].plot(t1, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[0,0].plot(t1, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin1)-1:
        if graph_fr:
            ax[0,0].plot(t1, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[0,0].plot(t1, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[0,0].plot(t1, val, color=colormap_begin[i])

if graph_fr:
    ax[0,0].plot(t1, soma_cli1, color='red', label=f"dans le soma", zorder=3)
else:
    ax[0,0].plot(t1, soma_cli1, color='red', label=f"in soma", zorder=3)


if graph_fr:
    ax[1,0].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)", fontsize=15)
else:
    ax[1,0].set_title(f"Chloride concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_end1):
    if i == 0:
        if graph_fr:
            ax[1,0].plot(t1, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[1,0].plot(t1, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end1)-1:
        if graph_fr:
            ax[1,0].plot(t1, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[1,0].plot(t1, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[1,0].plot(t1, val, color=colormap_end[-(i+1)])

ax[0,0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[0,0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[0,0].set_xlim(xlimit_min, xlimit_max)

ax[1,0].set_xlim(xlimit_min, xlimit_max)
ax[1,0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[1,0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[1,0].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[1,0].set_xlabel("Time (s)", fontsize=13)

ax[0,0].legend(loc='upper right', fontsize=14)
ax[1,0].legend(loc='upper right', fontsize=14)

del dend_cli_begin1
del dend_cli_end1
del soma_cli1






if graph_fr:
    ax[0,1].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[0,1].set_title(f"Chloride concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_begin2):
    if i == 0:
        if graph_fr:
            ax[0,1].plot(t2, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[0,1].plot(t2, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin2)-1:
        if graph_fr:
            ax[0,1].plot(t2, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[0,1].plot(t2, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[0,1].plot(t2, val, color=colormap_begin[i])

if graph_fr:
    ax[0,1].plot(t2, soma_cli2, color='red', label=f"dans le soma", zorder=3)
else:
    ax[0,1].plot(t2, soma_cli2, color='red', label=f"in soma", zorder=3)


if graph_fr:
    ax[1,1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)", fontsize=15)
else:
    ax[1,1].set_title(f"Chloride concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_end2):
    if i == 0:
        if graph_fr:
            ax[1,1].plot(t2, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[1,1].plot(t2, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end2)-1:
        if graph_fr:
            ax[1,1].plot(t2, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[1,1].plot(t2, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[1,1].plot(t2, val, color=colormap_end[-(i+1)])

ax[0,1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[0,1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[0,1].set_xlim(xlimit_min, xlimit_max)

ax[1,1].set_xlim(xlimit_min, xlimit_max)
ax[1,1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[1,1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[1,1].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[1,1].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin2
del dend_cli_end2
del soma_cli2









if graph_fr:
    ax[2,0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[2,0].set_title(f"Chloride concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_begin3):
    if i == 0:
        if graph_fr:
            ax[2,0].plot(t3, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[2,0].plot(t3, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin3)-1:
        if graph_fr:
            ax[2,0].plot(t3, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[2,0].plot(t3, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[2,0].plot(t3, val, color=colormap_begin[i])

if graph_fr:
    ax[2,0].plot(t3, soma_cli3, color='red', label=f"dans le soma", zorder=3)
else:
    ax[2,0].plot(t3, soma_cli3, color='red', label=f"in soma", zorder=3)


if graph_fr:
    ax[3,0].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)", fontsize=15)
else:
    ax[3,0].set_title(f"Chloride concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_end3):
    if i == 0:
        if graph_fr:
            ax[3,0].plot(t3, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[3,0].plot(t3, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end3)-1:
        if graph_fr:
            ax[3,0].plot(t3, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[3,0].plot(t3, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[3,0].plot(t3, val, color=colormap_end[-(i+1)])

ax[2,0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[2,0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[2,0].set_xlim(xlimit_min, xlimit_max)

ax[3,0].set_xlim(xlimit_min, xlimit_max)
ax[3,0].set_ylim(min_cl-0.025, max_cl+0.025)
ax[3,0].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[3,0].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[3,0].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin3
del dend_cli_end3
del soma_cli3






if graph_fr:
    ax[2,1].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)", fontsize=15)
else:
    ax[2,1].set_title(f"Chloride concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_begin4):
    if i == 0:
        if graph_fr:
            ax[2,1].plot(t4, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[2,1].plot(t4, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    elif i == len(dend_cli_begin4)-1:
        if graph_fr:
            ax[2,1].plot(t4, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[2,1].plot(t4, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    else:
        ax[2,1].plot(t4, val, color=colormap_begin[i])

if graph_fr:
    ax[2,1].plot(t4, soma_cli4, color='red', label=f"dans le soma", zorder=3)
else:
    ax[2,1].plot(t4, soma_cli4, color='red', label=f"in soma", zorder=3)


if graph_fr:
    ax[3,1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)", fontsize=15)
else:
    ax[3,1].set_title(f"Chloride concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm", fontsize=15)

for i, val in enumerate(dend_cli_end4):
    if i == 0:
        if graph_fr:
            ax[3,1].plot(t4, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
        else:
            ax[3,1].plot(t4, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
    elif i == len(dend_cli_end4)-1:
        if graph_fr:
            ax[3,1].plot(t4, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
        else:
            ax[3,1].plot(t4, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
    else:
        ax[3,1].plot(t4, val, color=colormap_end[-(i+1)])

ax[2,1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[2,1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
ax[2,1].set_xlim(xlimit_min, xlimit_max)

ax[3,1].set_xlim(xlimit_min, xlimit_max)
ax[3,1].set_ylim(min_cl-0.025, max_cl+0.025)
ax[3,1].set_ylabel(r"$[Cl^-]_i$ (mM)", fontsize=13)
if graph_fr:
    ax[3,1].set_xlabel("Temps (s)", fontsize=13)
else:
    ax[3,1].set_xlabel("Time (s)", fontsize=13)

del dend_cli_begin4
del dend_cli_end4
del soma_cli4


plt.show()