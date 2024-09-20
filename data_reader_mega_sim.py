import h5py
from matplotlib import rcParams
import matplotlib.pyplot as plt, cmap
import numpy as np
from function import show_info_sim


# Graph parameters -----------------------------
rcParams.update({'font.size': 22})
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'


# Loading the h5py dataset ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# INPUT
#f = h5py.file(r"Path of the dataset file", 'r')
path = r"dataset\voltclamped_-90mV_syn_nb_20_sim_lenght_520000_dt_(5,0.144,0.144)_L_150_kcc2_5.5e-05_nkcc1_5e-07_rnum=200_puffconc=1_2.0.hdf5"
f = h5py.File(path, 'r')


# To see the keys of the datatset
#print(list(f.keys()))


decal = 490000   # Offset to skip the initial stabilization of the simulation
graph_fr = False # If True, the graphs axis, titles and legends will be in french


# Graphs choices. Put 1 if you want the graph and 0 if not. -------------------------------------------------------------------
show_info = 0 # Print information on the simulation

chloride = 1  # Chloride intracellular concentration with multiple curves corresponding to different recording positions
potassium = 0 # Potassium intracellular concentration with multiple curves corresponding to different recording positions
sodium = 0    # Sodium intracellular concentration with multiple curves corresponding to different recording positions
gab = 0       # GABA extracellular concentration with multiple curves corresponding to different recording positions

mp_soma = 0          # 2 graphs : Membrane potential in soma and at one point in dendrite. Membrane potentials at each synapse.
current_soma = 0     # All currents in soma and at one point in dendrite
rev_pot_soma = 0     # Reversal potentials in soma and at one point in dendrite
rev_pot_dend_all = 0 # Reversal potentials at each of the 16 first synapses

icl_dend_all = 0            # Chloride currents at each of the 16 first synapses in dendrite
icl_dend_soma_add_check = 0 # Chloride currents in the soma and at one point in the dendrite
ik_dend_all = 0             # Potassium currents at each of the 16 first synapses in dendrite
ik_dend_soma_add_check = 0  # Potassium currents in the soma and at one point in the dendrite
ina_dend_all = 0            # Sodium currents at each of the 16 first synapses in dendrite
ina_dend_soma_add_check = 0 # Sodium currents in the soma and at one point in the dendrite
isynapse = 0                # Total current caused by each of the 16 first synapses

gcl_ghco3_grel = 0   # Chloride and HCO3- conductance at each of the 16 first synapses
open_states = 0      # Number of open channels at each synapses in the 3 states and combined
open_states_all = 0  # Number of open channels at each of the 16 first synapses


# Separation of the dataset in arrays ----------------------------------------------------
t_and_soma_v = f["mp_soma"] # Time and MP in soma
#print(len(t_and_soma_v[0]))

dend_v = f["mp_dend"]                 # MP on dendrite (at synapses)
dend_conc = f["conc_dend"]            # Concentrations in dendrite (at recording pos)
dend_e = f["e_dend"]                  # Reversal potentials in dendrite (at synapses)
dend_current_cl = f["syn_current_cl"] # Icl in dendrite (at synapses)
dend_current_k = f["syn_current_k"]   # Ik in dendrite (at synapses)
dend_current_na = f["syn_current_na"] # Ina in dendrite (at synapses)

soma_conc = f["conc_soma"]             # Concentrations in soma
soma_e = f["e_soma"]                   # Reversal potentials in soma
soma_current_cl = f["soma_current_cl"] # Icl in soma
soma_current_k = f["soma_current_k"]   # Ik in soma
soma_current_na = f["soma_current_na"] # Ina in soma

synapse_current = f["syn_current_other"] # Ihco3 and Igaba (Itot) in dendrite (at synapses)
gcl_and_other = f["syn_g_and_o"]         # gcl and ghco3 in dendrite (at synapses)


# Show info on the simulation and on the arrays -----------------------
if show_info == 1:
    show_info_sim(time_mp_soma=t_and_soma_v,
                                    mp_dend=dend_v,
                                    dend_concentration=dend_conc,
                                    soma_concentration=soma_conc,
                                    reversal_pot_soma=soma_e,
                                    reversal_pot_dend=dend_e,
                                    dend_currents_cl=dend_current_cl,
                                    dend_currents_k=dend_current_k,
                                    dend_currents_na=dend_current_na,
                                    soma_currents_cl=soma_current_cl,
                                    soma_currents_k=soma_current_k,
                                    soma_currents_na=soma_current_na,
                                    dend_currents_other=synapse_current,
                                    dend_g=gcl_and_other)


# Other values used later -----------------------------------------------------------------------------
record_pos = dend_conc.attrs["recording_locations"]    # Recording positions for the concentrations [-]
syn_pos = dend_conc.attrs["Synapses positions"]        # Synapses positions [-]
center_ind = int(len(syn_pos)/2-1)                     # Central synapse index [-]
lenght = dend_conc.attrs["dendrite_full_length_um"][0] # Dendrite (first part) lenght [µm]
puff_time = dend_conc.attrs["puff_time"]               # Time at wich the GABA puff occurs [ms]
puff_pos = dend_conc.attrs["puff_position"]            # Position of the GABA puff [-]
sim_time = dend_conc.attrs["simulation_length"]        # Simulation lenght [ms]
rnum = dend_conc.attrs["number_of_receptors"]          # Number of GABA receptors per synapse [-]


# Everything separated -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
t_ = t_and_soma_v[0]  # Time array [ms]
t = (t_ - decal)/1000 # Time array with shifted zero [s]

soma_v = t_and_soma_v[1] # Membrane potential in soma

# Concentration arrays in dendrite and in soma
dend_cli, dend_ki, dend_nai, dend_gabo = dend_conc[0], dend_conc[1], dend_conc[2], dend_conc[3]
soma_cli, soma_ki, soma_nai, soma_gabo = soma_conc[0], soma_conc[1], soma_conc[2], soma_conc[3]

# Chloride currents arrays in dendrite and in soma
dend_icl, dend_icl_kcc2, dend_icl_nkcc1, dend_icl_leak, dend_icl_synapses = dend_current_cl[0], dend_current_cl[1], dend_current_cl[2], dend_current_cl[3], dend_current_cl[4]
soma_icl, soma_icl_kcc2, soma_icl_nkcc1, soma_icl_leak = soma_current_cl[0], soma_current_cl[1], soma_current_cl[2], soma_current_cl[3]

# Potassium currents arrays in dendrite and in soma
dend_ik, dend_ik_kcc2, dend_ik_nkcc1, dend_ik_leak, dend_ik_nak, dend_ik_hh = dend_current_k[0], dend_current_k[1], dend_current_k[2], dend_current_k[3], dend_current_k[4], dend_current_k[5]
soma_ik, soma_ik_kcc2, soma_ik_nkcc1, soma_ik_leak, soma_ik_nak, soma_ik_hh = soma_current_k[0], soma_current_k[1], soma_current_k[2], soma_current_k[3], soma_current_k[4], soma_current_k[5]

# Sodium currents arrays in dendrite and in soma
dend_ina, dend_ina_nkcc1, dend_ina_leak, dend_ina_nak, dend_ina_hh = dend_current_na[0], dend_current_na[1], dend_current_na[2], dend_current_na[3], dend_current_na[4]
soma_ina, soma_ina_nkcc1, soma_ina_leak, soma_ina_nak, soma_ina_hh = soma_current_na[0], soma_current_na[1], soma_current_na[2], soma_current_na[3], soma_current_na[4]

# Reversal potentials arrays in dendrite and in soma
dend_ecl, dend_ek, dend_ena = dend_e[0], dend_e[1], dend_e[2]
soma_ecl, soma_ek, soma_ena = soma_e[0], soma_e[1], soma_e[2]

# HCO3 and total currents arrays in dendrite (at synapses)
dend_ihco3, dend_igaba = synapse_current[0], synapse_current[1]

# Chloride and HCO3 conductance at synapses. Open states of the GABA channels.
gcl, ghco3, grel, o1, o2, o3 = gcl_and_other[0], gcl_and_other[1], gcl_and_other[2], gcl_and_other[3], gcl_and_other[4], gcl_and_other[5]


# Find the central recording point (at puff) --
center = 100
for j in record_pos:
    if abs(puff_pos-j) < abs(puff_pos-center):
        center = j
#print(center)


# Separation of recording vectors between the ones that are from the puff toward the soma and
# the ones that are away from the soma ------------------------------------------------------
dend_cli_begin, dend_cli_end = [], []
dend_ki_begin, dend_ki_end = [], []
dend_nai_begin, dend_nai_end = [], []
dend_gabo_begin, dend_gabo_end = [], []
record_pos_begin, record_pos_end = [],[]
for i in range(len(record_pos)):
    if record_pos[i] <= center:
        dend_cli_begin.append(dend_cli[i])
        dend_ki_begin.append(dend_ki[i])
        dend_nai_begin.append(dend_nai[i])
        dend_gabo_begin.append(dend_gabo[i])
        record_pos_begin.append(record_pos[i])
    if record_pos[i] >= center:
        dend_cli_end.append(dend_cli[i])
        dend_ki_end.append(dend_ki[i])
        dend_nai_end.append(dend_nai[i])
        dend_gabo_end.append(dend_gabo[i])
        record_pos_end.append(record_pos[i])
print(r'Number of recording verctors toward the soma : ', len(dend_cli_begin))
print(r'Number of recording verctors away from soma  : ', len(dend_cli_end))


# Creation of the colormap and x-axis limits for the graphs---------------------------
colormap = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(record_pos)/2)))
colormap_syn_pos = cmap.Colormap("viridis")(np.linspace(0, 1, len(syn_pos)))
colormap_begin = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(dend_cli_begin))))
colormap_end = cmap.Colormap("viridis")(np.linspace(0, 1, int(len(dend_cli_end))))
xlimit_min = 0
xlimit_max = (sim_time-decal)/1000


# Deletion of the arrays that won't be used again
del t_and_soma_v
del dend_conc
del dend_e
del dend_current_cl
del dend_current_k
del dend_current_na
del soma_conc
del soma_e
del soma_current_cl
del soma_current_k
del soma_current_na
del synapse_current
del gcl_and_other



if chloride == 1:
    arg_max_cl = np.argmax(dend_cli[:, int(decal/5):])
    arg_min_cl = np.argmin(dend_cli[:, int(decal/5):])
    max_cl = max([np.max(dend_cli[:, int(decal/5):]), np.max(soma_cli[int(decal/5):])])
    min_cl = min([np.min(dend_cli[:, int(decal/5):]), np.min(soma_cli[int(decal/5):])])

    fig1, ax1 = plt.subplots(2,1)
    fig1.canvas.manager.set_window_title('[cl]')
    
    if graph_fr:
        ax1[0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax1[0].set_title(f"Chloride concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_cli_begin):
        if i == 0:
            if graph_fr:
                ax1[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax1[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(dend_cli_begin)-1:
            if graph_fr:
                ax1[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax1[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            ax1[0].plot(t, val, color=colormap_begin[i])

    if graph_fr:
        ax1[0].plot(t, soma_cli, color='red', label=f"dans le soma", zorder=3)
    else:
        ax1[0].plot(t, soma_cli, color='red', label=f"in soma", zorder=3)


    if graph_fr:
        ax1[1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)")
    else:
        ax1[1].set_title(f"Chloride concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_cli_end):
        if i == 0:
            if graph_fr:
                ax1[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax1[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        elif i == len(dend_cli_end)-1:
            if graph_fr:
                ax1[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax1[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        else:
            ax1[1].plot(t, val, color=colormap_end[-(i+1)])

    ax1[0].set_xlim(xlimit_min, xlimit_max)
    ax1[0].set_ylim(min_cl-0.025, max_cl+0.025)
    ax1[0].set_ylabel(r"$[Cl^-]_i$ (mM)")

    ax1[1].set_xlim(xlimit_min, xlimit_max)
    ax1[1].set_ylim(min_cl-0.025, max_cl+0.025)
    ax1[1].set_ylabel(r"$[Cl^-]_i$ (mM)")
    if graph_fr:
        ax1[1].set_xlabel("Temps (s)")
    else:
        ax1[1].set_xlabel("Time (s)")

    if arg_max_cl > arg_min_cl:
        ax1[0].legend(loc='lower right')
        ax1[1].legend(loc='lower right')
    else:
        ax1[0].legend(loc='upper right')
        ax1[1].legend(loc='upper right')
    
    # If you want to manually choose the location of the legend
    manual = False
    if manual:
        ax1[0].legend(loc='lower right')
        ax1[1].legend(loc='upper right')

del dend_cli_begin
del dend_cli_end
del soma_cli


if potassium == 1:
    arg_max_k = np.argmax(dend_ki[:, int(decal/5):])
    arg_min_k = np.argmin(dend_ki[:, int(decal/5):])
    max_k = max([np.max(dend_ki[:, int(decal/5):]), max(soma_ki[int(decal/5):])])
    min_k = min([np.min(dend_ki[:, int(decal/5):]), min(soma_ki[int(decal/5):])])

    fig2, ax2 = plt.subplots(2,1)
    fig2.canvas.manager.set_window_title('[k]')

    if graph_fr:
        ax2[0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax2[0].set_title(f"Potassium concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_ki_begin):
        if i == 0:
            if graph_fr:
                ax2[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax2[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(dend_ki_begin)-1:
            if graph_fr:
                ax2[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax2[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            ax2[0].plot(t, val, color=colormap_begin[i])

    if graph_fr:
        ax2[0].plot(t, soma_ki, color='red', label=f"dans le soma", zorder=3)
    else:
        ax2[0].plot(t, soma_ki, color='red', label=f"in soma", zorder=3)


    if graph_fr:
        ax2[1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans le dendrite)")
    else:
        ax2[1].set_title(f"Potassium concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_ki_end):
        if i == 0:
            if graph_fr:
                ax2[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax2[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        elif i == len(dend_ki_end)-1:
            if graph_fr:
                ax2[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax2[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        else:
            ax2[1].plot(t, val, color=colormap_end[-(i+1)])
    
    ax2[0].set_xlim(xlimit_min, xlimit_max)
    ax2[0].set_ylim(min_k-0.025, max_k+0.025)
    ax2[0].set_ylabel(r"$[K^+]_i$ (mM)")

    ax2[1].set_xlim(xlimit_min, xlimit_max)
    ax2[1].set_ylim(min_k-0.025, max_k+0.025)
    ax2[1].set_ylabel(r"$[K^+]_i$ (mM)")
    if graph_fr:
        ax2[1].set_xlabel("Temps (s)")
    else:
        ax2[1].set_xlabel("Time (s)")

    if arg_max_k > arg_min_k:
        ax2[0].legend(loc='lower right')
        ax2[1].legend(loc='lower right')
    else:
        ax2[0].legend(loc='upper right')
        ax2[1].legend(loc='upper right')
    
    # If you want to manually choose the location of the legend
    manual = False
    if manual:
        ax2[0].legend(loc='lower right')
        ax2[1].legend(loc='upper right')

del dend_ki_begin
del dend_ki_end
del soma_ki


if sodium == 1:
    arg_max_na = np.argmax(dend_nai[:, int(decal/5):])
    arg_min_na = np.argmin(dend_nai[:, int(decal/5):])
    max_na = max([np.max(dend_nai[:, int(decal/5):]), max(soma_nai[int(decal/5):])])
    min_na = min([np.min(dend_nai[:, int(decal/5):]), min(soma_nai[int(decal/5):])])

    fig4, ax4 = plt.subplots(2,1)
    fig4.canvas.manager.set_window_title('[Na]')

    if graph_fr:
        ax4[0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax4[0].set_title(f"Sodium concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_nai_begin):
        if i == 0:
            if graph_fr:
                ax4[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax4[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(dend_nai_begin)-1:
            if graph_fr:
                ax4[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax4[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            ax4[0].plot(t, val, color=colormap_begin[i])

    if graph_fr:
        ax4[0].plot(t, soma_nai, color='red', label=f"dans le soma", zorder=3)
    else:
        ax4[0].plot(t, soma_nai, color='red', label=f"in soma", zorder=3)


    if graph_fr:
        ax4[1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax4[1].set_title(f"Sodium concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_nai_end):
        if i == 0:
            if graph_fr:
                ax4[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                ax4[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        elif i == len(dend_nai_end)-1:
            if graph_fr:
                ax4[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                ax4[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        else:
            ax4[1].plot(t, val, color=colormap_end[-(i+1)])

    ax4[0].set_xlim(xlimit_min, xlimit_max)
    ax4[0].set_ylim(min_na-0.025, max_na+0.025)
    ax4[0].set_ylabel(r"$[Na^+]_i$ (mM)")

    ax4[1].set_xlim(xlimit_min, xlimit_max)
    ax4[1].set_ylim(min_na-0.025, max_na+0.025)
    ax4[1].set_ylabel(r"$[Na^+]_i$ (mM)")
    if graph_fr:
        ax4[1].set_xlabel("Temps (s)")
    else:
        ax4[1].set_xlabel("Time (s)")

    if arg_max_na > arg_min_na:
        ax4[0].legend(loc='lower right')
        ax4[1].legend(loc='lower right')
    else:
        ax4[0].legend(loc='upper right')
        ax4[1].legend(loc='upper right')
    
    # If you want to manually choose the location of the legend
    manual = False
    if manual:
        ax4[0].legend(loc='lower right')
        ax4[1].legend(loc='upper right')

del dend_nai_begin
del dend_nai_end
del soma_nai


if gab == 1:
    plt.figure(f"GABA external concentration {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm")
    for i, val in enumerate(dend_gabo_begin):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':')
            else:
                plt.plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':')
        elif i == len(dend_gabo_begin)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--')
            else:
                plt.plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--')
        else:
            plt.plot(t, val, color=colormap_begin[i])
    if graph_fr:
        plt.plot(t, soma_gabo, color='red', label=f"dans le soma")
    else:
        plt.plot(t, soma_gabo, color='red', label=f"in soma")

    plt.xlim((puff_time-decal)/1000-0.01, (puff_time-decal)/1000+0.1)
    plt.ylabel(r"$[GABA]_o$ (mM)")
    if graph_fr:
        plt.xlabel("Temps (s)")
    else:
        plt.xlabel("Time (s)")
    plt.legend(loc='upper right')
    plt.tight_layout()


    plt.figure(f"GABA external concentration {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm")
    for i, val in enumerate(dend_gabo_end):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--')
            else:
                plt.plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--')
        elif i == len(dend_gabo_end)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':')
            else:
                plt.plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':')
        else:
            plt.plot(t, val, color=colormap_end[i])

    plt.xlim((puff_time-decal)/1000-0.01, (puff_time-decal)/1000+0.1)
    plt.ylabel(r"$[GABA]_o$ (mM)")
    if graph_fr:
        plt.xlabel("Temps (s)")
    else:
        plt.xlabel("Time (s)")
    plt.legend(loc='upper right')
    plt.tight_layout()


    fig3, ax3 = plt.subplots(2,1)
    fig3.canvas.manager.set_window_title('gab')

    if graph_fr:
        ax3[0].set_title(f"Vers le soma (de {record_pos_begin[0]*lenght:.1f} à {record_pos_begin[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax3[0].set_title(f"GABA extracellular concentration from {record_pos_begin[0]*lenght:.1f} to {record_pos_begin[-1]*lenght:.1f} µm")

    for i, val in enumerate(dend_gabo_begin):
        if i == 0:
            if graph_fr:
                ax3[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle=':')
            else:
                ax3[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle=':')
        elif i == len(dend_gabo_begin)-1:
            if graph_fr:
                ax3[0].plot(t, val, color='black', label=f"à {record_pos_begin[i]*lenght:.1f} µm dans la dendrite", linestyle='--')
            else:
                ax3[0].plot(t, val, color='black', label=f"at {record_pos_begin[i]*lenght:.1f} µm in dendrite", linestyle='--')
        else:
            ax3[0].plot(t, val, color=colormap_begin[i])

    if graph_fr:
        ax3[0].plot(t, soma_gabo, color='red', label=f"dans le soma")
    else:
        ax3[0].plot(t, soma_gabo, color='red', label=f"in soma")


    if graph_fr:
        ax3[1].set_title(f"Vers l'opposé du soma (de {record_pos_end[0]*lenght:.1f} à {record_pos_end[-1]*lenght:.1f} µm dans la dendrite)")
    else:
        ax3[1].set_title(f"GABA extracellular concentration from {record_pos_end[0]*lenght:.1f} to {record_pos_end[-1]*lenght:.1f} µm")
    for i, val in enumerate(dend_gabo_end):
        if i == 0:
            if graph_fr:
                ax3[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle='--')
            else:
                ax3[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle='--')
        elif i == len(dend_gabo_end)-1:
            if graph_fr:
                ax3[1].plot(t, val, color='black', label=f"à {record_pos_end[i]*lenght:.1f} µm dans la dendrite", linestyle=':')
            else:
                ax3[1].plot(t, val, color='black', label=f"at {record_pos_end[i]*lenght:.1f} µm in dendrite", linestyle=':')
        else:
            ax3[1].plot(t, val, color=colormap_end[i])

    ax3[0].set_xlim((puff_time-decal)/1000-0.01, (puff_time-decal)/1000+0.1)
    ax3[1].set_xlim((puff_time-decal)/1000-0.01, (puff_time-decal)/1000+0.1)
    if graph_fr:
        ax3[1].set_xlabel("Temps (s)")
    else:
        ax3[1].set_xlabel("Time (s)")
    ax3[0].set_ylabel(r"$[GABA]_o$ (mM)")
    ax3[1].set_ylabel(r"$[GABA]_o$ (mM)")
    ax3[0].legend(loc='upper right')
    ax3[1].legend(loc='upper right')

del dend_gabo_begin
del dend_gabo_end
del soma_gabo


if mp_soma == 1:
    plt.figure('MP_all_curves')
    for i in range(len(dend_v)):
        if i == 0:
            if graph_fr:
                plt.plot(t, dend_v[i], label=f'à {syn_pos[i]*lenght:.2f} µm dans la dendrite', linestyle='--', color='black', zorder=3)
            else:
                plt.plot(t, dend_v[i], label=f'at {syn_pos[i]*lenght:.2f} µm in dendrite', linestyle='--', color='black', zorder=3)
        elif i == len(dend_v)-1:
            if graph_fr:
                plt.plot(t, dend_v[i], label=f'à {syn_pos[i]*lenght:.2f} µm dans la dendrite', linestyle=':', color='black', zorder=3)
            else:
                plt.plot(t, dend_v[i], label=f'at {syn_pos[i]*lenght:.2f} µm in dendrite', linestyle=':', color='black', zorder=3)
        else:
            plt.plot(t, dend_v[i], color=colormap_syn_pos[i])
    if graph_fr:
        plt.plot(t, soma_v, label=f'dans le soma', color='red')
        plt.xlabel('Temps [s]')
    else:
        plt.plot(t, soma_v, label=f'in soma', color='red')
        plt.xlabel('Time [s]')

    plt.ylabel('Membrane potential [mV]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend()


    plt.figure('MP_one_curve')
    if graph_fr:
        plt.plot(t, dend_v[center_ind], label=f'à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite')
        plt.plot(t, soma_v, label=f'dans le soma')
        plt.xlabel('Temps [ms]')
        plt.ylabel('Potential de membrane [mV]')
    else:
        plt.plot(t, dend_v[center_ind], label=f'at {syn_pos[center_ind]*lenght:.2f} µm in dendrite')
        plt.plot(t, soma_v, label=f'in soma')
        plt.xlabel('Time [ms]')
        plt.ylabel('Membrane potential [mV]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend()

if current_soma == 1:
    plt.figure('Ionic_currents_dend_one_point')
    if graph_fr:
        plt.plot(t, dend_icl[center_ind], label=r'$I_{cl}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="orange")
        plt.plot(t, dend_ik[center_ind], label=r'$I_{k}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="darkgreen")
        plt.plot(t, dend_ina[center_ind], label=r'$I_{na}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="deepskyblue")
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.plot(t, dend_icl[center_ind], label=r'$I_{cl}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="orange")
        plt.plot(t, dend_ik[center_ind], label=r'$I_{k}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="darkgreen")
        plt.plot(t, dend_ina[center_ind], label=r'$I_{na}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="deepskyblue")
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend()

    plt.figure('Ionic_currents_soma')
    if graph_fr:
        plt.plot(t, soma_icl, label=r'$I_{cl}$ dans le soma', color="orange")
        plt.plot(t, soma_ik, label=r'$I_{k}$ dans le soma', color="darkgreen")
        plt.plot(t, soma_ina, label=r'$I_{na}$ dans le soma', color="deepskyblue")
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.plot(t, soma_icl, label=r'$I_{cl}$ in soma', color="orange")
        plt.plot(t, soma_ik, label=r'$I_{k}$ in soma', color="darkgreen")
        plt.plot(t, soma_ina, label=r'$I_{na}$ in soma', color="deepskyblue")
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    #plt.plot(t, soma_icl + soma_ik + soma_ina, label=r"$I_{cl}+I_{k}+I_{na}$", color='black', linestyle='--')
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend()

if rev_pot_soma == 1:
    if graph_fr:
        plt.figure('Reversal_potentials')
        plt.plot(t, dend_ecl[center_ind], label=r'$E_{cl}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="red")
        plt.plot(t, dend_ek[center_ind], label=r'$E_{k}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="darkgreen")
        plt.plot(t, dend_ena[center_ind], label=r'$E_{na}$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="deepskyblue")
        plt.plot(t, dend_v[center_ind], label=r'$MP$' + f' à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite', color="orange")

        plt.plot(t, soma_ecl, label=r'$E_{cl}$ dans le soma', color="red", linestyle='--')
        plt.plot(t, soma_ek, label=r'$E_{k}$ dans le soma', color="darkgreen", linestyle='--')
        plt.plot(t, soma_ena, label=r'$E_{na}$ dans le soma', color="deepskyblue", linestyle='--')
        plt.plot(t, soma_v, label=r'$MP$ dans le soma', color="orange", linestyle='--')

        plt.xlabel('Temps [s]')
        plt.ylabel('Potentiel réversible [mV]')
    else:
        plt.figure('Reversal_potentials')
        plt.plot(t, dend_ecl[center_ind], label=r'$E_{cl}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="red")
        plt.plot(t, dend_ek[center_ind], label=r'$E_{k}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="darkgreen")
        plt.plot(t, dend_ena[center_ind], label=r'$E_{na}$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="deepskyblue")
        plt.plot(t, dend_v[center_ind], label=r'$MP$' + f' at {syn_pos[center_ind]*lenght:.2f} µm in dendrite', color="orange")

        plt.plot(t, soma_ecl, label=r'$E_{cl}$ in soma', color="red", linestyle='--')
        plt.plot(t, soma_ek, label=r'$E_{k}$ in soma', color="darkgreen", linestyle='--')
        plt.plot(t, soma_ena, label=r'$E_{na}$ in soma', color="deepskyblue", linestyle='--')
        plt.plot(t, soma_v, label=r'$MP$ in soma', color="orange", linestyle='--')

        plt.xlabel('Time [s]')
        plt.ylabel('Reversal potential [mV]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')


if rev_pot_dend_all == 1:
    max_e = max([np.max(dend_ecl), np.max(dend_ek), np.max(dend_ena), np.max(dend_v)])
    min_e = min([np.min(dend_ecl), np.min(dend_ek), np.min(dend_ena), np.min(dend_v)])
    max_e_without_ena = max([np.max(dend_ecl), np.max(dend_ek), np.max(dend_v)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('Reversal potentials at synapses')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, dend_v[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_v[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_v[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_v[j+12], color='orange', alpha=0.7)

        ax_[0, j].plot(t, dend_ek[j], color='deepskyblue', alpha=0.7)
        ax_[1, j].plot(t, dend_ek[j+4], color='deepskyblue', alpha=0.7)
        ax_[2, j].plot(t, dend_ek[j+8], color='deepskyblue', alpha=0.7)
        ax_[3, j].plot(t, dend_ek[j+12], color='deepskyblue', alpha=0.7)

        ax_[0, j].plot(t, dend_ena[j], color='orchid', alpha=0.7)
        ax_[1, j].plot(t, dend_ena[j+4], color='orchid', alpha=0.7)
        ax_[2, j].plot(t, dend_ena[j+8], color='orchid', alpha=0.7)
        ax_[3, j].plot(t, dend_ena[j+12], color='orchid', alpha=0.7)

        ax_[0, j].plot(t, dend_ecl[j], color='red', alpha=0.7)
        ax_[1, j].plot(t, dend_ecl[j+4], color='red', alpha=0.7)
        ax_[2, j].plot(t, dend_ecl[j+8], color='red', alpha=0.7)
        ax_[3, j].plot(t, dend_ecl[j+12], color='red', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Potentiel\n[mV]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_e-10, max_e+10)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Potential\n[mV]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_e-10, max_e+10)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)



    # Même graphique sans Ena
    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('Reversal potentials at synapses (without E_na)')
    for j in range(4):
        ax_[0, j].plot(t, dend_v[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_v[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_v[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_v[j+12], color='orange', alpha=0.7)

        ax_[0, j].plot(t, dend_ek[j], color='deepskyblue', alpha=0.7)
        ax_[1, j].plot(t, dend_ek[j+4], color='deepskyblue', alpha=0.7)
        ax_[2, j].plot(t, dend_ek[j+8], color='deepskyblue', alpha=0.7)
        ax_[3, j].plot(t, dend_ek[j+12], color='deepskyblue', alpha=0.7)

        ax_[0, j].plot(t, dend_ecl[j], color='red', alpha=0.7)
        ax_[1, j].plot(t, dend_ecl[j+4], color='red', alpha=0.7)
        ax_[2, j].plot(t, dend_ecl[j+8], color='red', alpha=0.7)
        ax_[3, j].plot(t, dend_ecl[j+12], color='red', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Potentiel\n[mV]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_e-10, max_e_without_ena+10)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Potential\n[mV]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_e-10, max_e_without_ena+10)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)


    b = np.zeros_like(dend_ecl[0])
    plt.figure("Legend_of_graphs")
    plt.plot(t, b, color='orange', label=r'$MP$')
    plt.plot(t, b, color='deepskyblue', label=r'$E_{K}$')
    plt.plot(t, b, color='orchid', label=r'$E_{Na}$')
    plt.plot(t, b, color='red', label=r'$E_{Cl}$')
    plt.ylim(1,2)
    plt.legend(loc='center')

del soma_v
del dend_v
del dend_ecl
del soma_ecl
del dend_ek
del soma_ek
del dend_ena
del soma_ena


if icl_dend_all == 1:
    max_icl_big_currents = max([np.max(dend_icl), np.max(dend_icl_synapses), max(soma_icl)])
    min_icl_big_currents = min([np.min(dend_icl), np.min(dend_icl_synapses), min(soma_icl)])
    max_icl_smaller_currents = max([np.max(dend_icl_kcc2), np.max(dend_icl_nkcc1), np.max(dend_icl_leak)])
    min_icl_smaller_currents = min([np.min(dend_icl_kcc2), np.min(dend_icl_nkcc1), np.min(dend_icl_leak)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('I_cl at synapses (bigger currents)')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, dend_icl[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_icl[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_icl[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_icl[j+12], color='orange', alpha=0.7)

        ax_[0, j].plot(t, dend_icl_synapses[j], color='chocolate', alpha=0.7)
        ax_[1, j].plot(t, dend_icl_synapses[j+4], color='chocolate', alpha=0.7)
        ax_[2, j].plot(t, dend_icl_synapses[j+8], color='chocolate', alpha=0.7)
        ax_[3, j].plot(t, dend_icl_synapses[j+12], color='chocolate', alpha=0.7)
    ax_[0,0].plot(t, soma_icl, linestyle='--', color='red')

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Current\n[pA]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_icl_big_currents, max_icl_big_currents+1)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Courant\n[pA]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_icl_big_currents, max_icl_big_currents+1)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)


    a = np.zeros_like(dend_icl_synapses[0])
    plt.figure('I_cl at synapses (bigger currents) label')
    plt.plot(t, a, color='orange', label=r'$I_{cl,tot}$')
    plt.plot(t, a, color='chocolate', label=r'$I_{cl,synapse}$')
    plt.plot(t, a, color='red', label=r'$I_{cl, soma}$')
    plt.ylim(1,2)
    plt.legend(loc='center')


    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('I_cl at synapses (smaller currents)')
    for j in range(4):
        ax_[0, j].plot(t, dend_icl_nkcc1[j], color='deepskyblue', alpha=0.7)
        ax_[1, j].plot(t, dend_icl_nkcc1[j+4], color='deepskyblue', alpha=0.7)
        ax_[2, j].plot(t, dend_icl_nkcc1[j+8], color='deepskyblue', alpha=0.7)
        ax_[3, j].plot(t, dend_icl_nkcc1[j+12], color='deepskyblue', alpha=0.7)

        ax_[0, j].plot(t, dend_icl_kcc2[j], color='darkgreen', alpha=0.7)
        ax_[1, j].plot(t, dend_icl_kcc2[j+4], color='darkgreen', alpha=0.7)
        ax_[2, j].plot(t, dend_icl_kcc2[j+8], color='darkgreen', alpha=0.7)
        ax_[3, j].plot(t, dend_icl_kcc2[j+12], color='darkgreen', alpha=0.7)

        ax_[0, j].plot(t, dend_icl_leak[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_icl_leak[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_icl_leak[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_icl_leak[j+12], color='orange', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Current\n[pA]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_icl_smaller_currents, max_icl_smaller_currents)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Courant\n[pA]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_icl_smaller_currents, max_icl_smaller_currents)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    

    a = np.zeros_like(dend_icl_synapses[0])
    plt.figure('I_cl at synapses (smaller currents) label')
    plt.plot(t, a, color='orange', label=r'$I_{cl,leak}$')
    plt.plot(t, a, color='deepskyblue', label=r'$I_{cl,nkcc1}$')
    plt.plot(t, a, color='darkgreen', label=r'$I_{cl,kcc2}$')
    plt.ylim(1,2)
    plt.legend(loc='center')

if icl_dend_soma_add_check == 1:
    plt.figure(f'Chloride_currents_dend')
    if graph_fr:
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
        plt.title(f'Courants de chlore à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite')
    else:
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
        plt.title(f'Chloride currents at {syn_pos[center_ind]*lenght:.2f} µm in dendrite')
    plt.plot(t, dend_icl[center_ind], label=r'$I_{cl}$', color="darkorange", lw=2)
    plt.plot(t, dend_icl_kcc2[center_ind], label=r'$I_{cl,kcc2}$', color="mediumblue", lw=2)
    plt.plot(t, dend_icl_nkcc1[center_ind], label=r'$I_{cl,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, dend_icl_leak[center_ind], label=r'$I_{cl,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, dend_icl_synapses[center_ind], label=r'$I_{cl,synapses}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, dend_icl_leak[center_ind] + dend_icl_synapses[center_ind] + dend_icl_nkcc1[center_ind] + dend_icl_kcc2[center_ind],
            label=r'$I_{cl,kcc2}+I_{cl,nkcc1}+$' + '\n' + r'$I_{cl,leak}+I_{cl,synapse}$', color="black",
            linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')


    plt.figure(f'Chloride_currents_dend (zoom)')
    if graph_fr:
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
        plt.title(f'Courants de chlore à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite')
    else:
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
        plt.title(f'Chloride currents at {syn_pos[center_ind]*lenght:.2f} µm in dendrite')
    plt.plot(t, dend_icl[center_ind], label=r'$I_{cl}$', color="darkorange", lw=2)
    plt.plot(t, dend_icl_kcc2[center_ind], label=r'$I_{cl,kcc2}$', color="mediumblue", lw=2)
    plt.plot(t, dend_icl_nkcc1[center_ind], label=r'$I_{cl,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, dend_icl_leak[center_ind], label=r'$I_{cl,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, dend_icl_synapses[center_ind], label=r'$I_{cl,synapses}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, dend_icl_leak[center_ind] + dend_icl_synapses[center_ind] + dend_icl_nkcc1[center_ind] + dend_icl_kcc2[center_ind],
            label=r'$I_{cl,kcc2}+I_{cl,nkcc1}+$' + '\n' + r'$I_{cl,leak}+I_{cl,synapse}$', color="black",
            linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(0, 8)
    plt.ylim(-0.02, 0.2)
    plt.legend(loc='center right')


    plt.figure('Chloride_currents_soma')
    if graph_fr:
        plt.title(f'Courants de chlore dans le soma')
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.title(f'Chloride currents in soma')
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.plot(t, soma_icl, label=r'$I_{cl}$', color="darkorange", lw=2)
    plt.plot(t, soma_icl_kcc2, label=r'$I_{cl,kcc2}$', color="mediumblue", lw=2)
    plt.plot(t, soma_icl_nkcc1, label=r'$I_{cl,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, soma_icl_leak, label=r'$I_{cl,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, soma_icl_leak + soma_icl_nkcc1 + soma_icl_kcc2,
            label=r'$I_{cl,kcc2}+I_{cl,nkcc1}+I_{cl,leak}$', color="black", linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')


if ik_dend_all == 1:
    max_ik = max([np.max(dend_ik), np.max(dend_ik_kcc2), np.max(dend_ik_nkcc1), np.max(dend_ik_leak), np.max(dend_ik_nak), np.max(dend_ik_hh)])
    min_ik = min([np.min(dend_ik), np.min(dend_ik_kcc2), np.min(dend_ik_nkcc1), np.min(dend_ik_leak), np.min(dend_ik_nak), np.min(dend_ik_hh)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('I_k at synapses')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, dend_ik[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_ik[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_ik[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_ik[j+12], color='orange', alpha=0.7)

        ax_[0, j].plot(t, dend_ik_nkcc1[j], color='deepskyblue', alpha=0.7)
        ax_[1, j].plot(t, dend_ik_nkcc1[j+4], color='deepskyblue', alpha=0.7)
        ax_[2, j].plot(t, dend_ik_nkcc1[j+8], color='deepskyblue', alpha=0.7)
        ax_[3, j].plot(t, dend_ik_nkcc1[j+12], color='deepskyblue', alpha=0.7)

        ax_[0, j].plot(t, dend_ik_kcc2[j], color='darkgreen', alpha=0.7)
        ax_[1, j].plot(t, dend_ik_kcc2[j+4], color='darkgreen', alpha=0.7)
        ax_[2, j].plot(t, dend_ik_kcc2[j+8], color='darkgreen', alpha=0.7)
        ax_[3, j].plot(t, dend_ik_kcc2[j+12], color='darkgreen', alpha=0.7)

        ax_[0, j].plot(t, dend_ik_leak[j], color='orchid', alpha=0.7)
        ax_[1, j].plot(t, dend_ik_leak[j+4], color='orchid', alpha=0.7)
        ax_[2, j].plot(t, dend_ik_leak[j+8], color='orchid', alpha=0.7)
        ax_[3, j].plot(t, dend_ik_leak[j+12], color='orchid', alpha=0.7)

        ax_[0, j].plot(t, dend_ik_nak[j], color='red', alpha=0.7)
        ax_[1, j].plot(t, dend_ik_nak[j+4], color='red', alpha=0.7)
        ax_[2, j].plot(t, dend_ik_nak[j+8], color='red', alpha=0.7)
        ax_[3, j].plot(t, dend_ik_nak[j+12], color='red', alpha=0.7)
        
        ax_[0, j].plot(t, dend_ik_hh[j], color='chocolate', alpha=0.7)
        ax_[1, j].plot(t, dend_ik_hh[j+4], color='chocolate', alpha=0.7)
        ax_[2, j].plot(t, dend_ik_hh[j+8], color='chocolate', alpha=0.7)
        ax_[3, j].plot(t, dend_ik_hh[j+12], color='chocolate', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Courant\n[pA]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_ik, max_ik)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Current\n[pA]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_ik, max_ik)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)


    c = np.zeros_like(dend_ik[0])
    plt.figure("Ik_label")
    plt.plot(t, c, color='orange', label=r'$I_{K,tot}$')
    plt.plot(t, c, color='deepskyblue', label=r'$I_{K,nkc1}$')
    plt.plot(t, c, color='darkgreen', label=r'$I_{K,kcc2}$')
    plt.plot(t, c, color='orchid', label=r'$I_{K,leak}$')
    plt.plot(t, c, color='red', label=r'$I_{K,NaK}$')
    plt.plot(t, c, color='chocolate', label=r'$I_{K,hh}$')
    plt.ylim(1,2)
    plt.legend(loc='center')

if ik_dend_soma_add_check == 1:
    plt.figure(f'Potassium_currents_dend')
    if graph_fr:
        plt.title(f'Courants de potassium à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite')
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.title(f'Potassium currents at {syn_pos[center_ind]*lenght:.2f} µm in dendrite')
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.plot(t, dend_ik[center_ind], label=r'$I_{k}$', color="darkorange", lw=2)
    plt.plot(t, dend_ik_kcc2[center_ind], label=r'$I_{k,kcc2}$', color="mediumblue", lw=2)
    plt.plot(t, dend_ik_nkcc1[center_ind], label=r'$I_{k,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, dend_ik_leak[center_ind], label=r'$I_{k,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, dend_ik_nak[center_ind], label=r'$I_{k,nak}$', color="forestgreen", lw=2, alpha=1)
    plt.plot(t, dend_ik_hh[center_ind], label=r'$I_{k,hh}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, dend_ik_leak[center_ind] + dend_ik_nkcc1[center_ind] + dend_ik_kcc2[center_ind] + dend_ik_nak[center_ind] + dend_ik_hh[center_ind],
            label=r'$I_{k,kcc2}+I_{k,nkcc1}+$' + '\n' + r'$I_{k,leak}+I_{k,nak}+I_{k,hh}$', color="black",
            linestyle='--', alpha=0.3, zorder=3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')

    plt.figure(f'Potassium_currents_soma')
    if graph_fr:
        plt.title(f'Courants de potassium dans le soma')
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.title(f'Potassium currents in soma')
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.plot(t, soma_ik, label=r'$I_{k}$', color="darkorange", lw=2)
    plt.plot(t, soma_ik_kcc2, label=r'$I_{k,kcc2}$', color="mediumblue", lw=2)
    plt.plot(t, soma_ik_nkcc1, label=r'$I_{k,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, soma_ik_leak, label=r'$I_{k,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, soma_ik_nak, label=r'$I_{k,nak}$', color="forestgreen", lw=2, alpha=1)
    plt.plot(t, soma_ik_hh, label=r'$I_{k,hh}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, soma_ik_leak + soma_ik_nkcc1 + soma_ik_kcc2 + soma_ik_nak + soma_ik_hh,
            label=r'$I_{k,kcc2}+I_{k,nkcc1}+$' + '\n' + r'$I_{k,leak}+I_{k,nak}+I_{k,hh}$', color="black",
            linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')


if ina_dend_all == 1:
    max_ina = max([np.max(dend_ina), np.max(dend_ina_nkcc1), np.max(dend_ina_leak), np.max(dend_ina_nak), np.max(dend_ina_hh)])
    min_ina = min([np.min(dend_ina), np.min(dend_ina_nkcc1), np.min(dend_ina_leak), np.min(dend_ina_nak), np.min(dend_ina_hh)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('I_na at synapses')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, dend_ina[j], color='orange', alpha=0.7)
        ax_[1, j].plot(t, dend_ina[j+4], color='orange', alpha=0.7)
        ax_[2, j].plot(t, dend_ina[j+8], color='orange', alpha=0.7)
        ax_[3, j].plot(t, dend_ina[j+12], color='orange', alpha=0.7)

        ax_[0, j].plot(t, dend_ina_nkcc1[j], color='deepskyblue', alpha=0.7)
        ax_[1, j].plot(t, dend_ina_nkcc1[j+4], color='deepskyblue', alpha=0.7)
        ax_[2, j].plot(t, dend_ina_nkcc1[j+8], color='deepskyblue', alpha=0.7)
        ax_[3, j].plot(t, dend_ina_nkcc1[j+12], color='deepskyblue', alpha=0.7)

        ax_[0, j].plot(t, dend_ina_leak[j], color='orchid', alpha=0.7)
        ax_[1, j].plot(t, dend_ina_leak[j+4], color='orchid', alpha=0.7)
        ax_[2, j].plot(t, dend_ina_leak[j+8], color='orchid', alpha=0.7)
        ax_[3, j].plot(t, dend_ina_leak[j+12], color='orchid', alpha=0.7)

        ax_[0, j].plot(t, dend_ina_nak[j], color='red', alpha=0.7)
        ax_[1, j].plot(t, dend_ina_nak[j+4], color='red', alpha=0.7)
        ax_[2, j].plot(t, dend_ina_nak[j+8], color='red', alpha=0.7)
        ax_[3, j].plot(t, dend_ina_nak[j+12], color='red', alpha=0.7)

        ax_[0, j].plot(t, dend_ina_hh[j], color='chocolate', alpha=0.7)
        ax_[1, j].plot(t, dend_ina_hh[j+4], color='chocolate', alpha=0.7)
        ax_[2, j].plot(t, dend_ina_hh[j+8], color='chocolate', alpha=0.7)
        ax_[3, j].plot(t, dend_ina_hh[j+12], color='chocolate', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Courant\n[pA]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_ina, max_ina)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Current\n[pA]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_ina, max_ina)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)


    d = np.zeros_like(dend_ina[0])
    plt.figure("Ina_label")
    plt.plot(t, d, color='orange', label=r'$I_{Na,tot}$')
    plt.plot(t, d, color='deepskyblue', label=r'$I_{Na,nkcc1}$')
    plt.plot(t, d, color='orchid', label=r'$I_{Na,leak}$')
    plt.plot(t, d, color='red', label=r'$I_{Na,Nak}$')
    plt.plot(t, d, color='chocolate', label=r'$I_{Na,hh}$')
    plt.ylim(1,2)
    plt.legend(loc='center')

if ina_dend_soma_add_check == 1:
    plt.figure(f'Sodium_currents_dend')
    if graph_fr:
        plt.title(f'Courants de sodium à {syn_pos[center_ind]*lenght:.2f} µm dans la dendrite')
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.title(f'Sodium currents at {syn_pos[center_ind]*lenght:.2f} µm in dendrite')
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.plot(t, dend_ina[center_ind], label=r'$I_{na}$', color="darkorange", lw=2)
    plt.plot(t, dend_ina_nkcc1[center_ind], label=r'$I_{na,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, dend_ina_leak[center_ind], label=r'$I_{na,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, dend_ina_nak[center_ind], label=r'$I_{na,nak}$', color="forestgreen", lw=2, alpha=1)
    plt.plot(t, dend_ina_hh[center_ind], label=r'$I_{na,hh}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, dend_ina_leak[center_ind] + dend_ina_nkcc1[center_ind] + dend_ina_nak[center_ind] + dend_ina_hh[center_ind],
            label=r'$I_{na,nkcc1}+I_{na,leak}+$' + "\n" + r"$I_{na,nak}+I_{na,hh}$", color="black",
            linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='lower right')

    plt.figure(f'Sodium_currents_soma')
    if graph_fr:
        plt.title(f'Courants de sodium dans le soma')
        plt.xlabel('Temps [s]')
        plt.ylabel('Courant ionique [pA]')
    else:
        plt.title(f'Sodium currents in soma')
        plt.xlabel('Time [s]')
        plt.ylabel('Ionic current [pA]')
    plt.plot(t, soma_ina, label=r'$I_{na}$', color="darkorange", lw=2)
    plt.plot(t, soma_ina_nkcc1, label=r'$I_{na,nkcc1}$', color="darkkhaki", lw=2)
    plt.plot(t, soma_ina_leak, label=r'$I_{na,leak}$', color="darkorchid", lw=2, alpha=1)
    plt.plot(t, soma_ina_nak, label=r'$I_{na,nak}$', color="forestgreen", lw=2, alpha=1)
    plt.plot(t, soma_ina_hh, label=r'$I_{na,hh}$', color="deepskyblue", lw=2, alpha=1)
    plt.plot(t, soma_ina_leak + soma_ina_nkcc1 + soma_ina_nak + soma_ina_hh,
            label=r'$I_{na,nkcc1}+I_{na,leak}+$' + '\n' + r'$I_{na,nak}+I_{na,hh}$', color="black",
            linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='lower right')


if isynapse == 1:
    max_isyn = max([np.max(dend_igaba), np.max(dend_ihco3), np.max(dend_icl_synapses)])
    min_isyn = min([np.min(dend_igaba), np.min(dend_ihco3), np.min(dend_icl_synapses)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('I at synapses')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, dend_igaba[j], color='black', alpha=0.4, linestyle='--', zorder=3)
        ax_[1, j].plot(t, dend_igaba[j+4], color='black', alpha=0.4, linestyle='--', zorder=3)
        ax_[2, j].plot(t, dend_igaba[j+8], color='black', alpha=0.4, linestyle='--', zorder=3)
        ax_[3, j].plot(t, dend_igaba[j+12], color='black', alpha=0.4, linestyle='--', zorder=3)

        ax_[0, j].plot(t, dend_ihco3[j], color='deepskyblue', alpha=1)
        ax_[1, j].plot(t, dend_ihco3[j+4], color='deepskyblue', alpha=1)
        ax_[2, j].plot(t, dend_ihco3[j+8], color='deepskyblue', alpha=1)
        ax_[3, j].plot(t, dend_ihco3[j+12], color='deepskyblue', alpha=1)

        ax_[0, j].plot(t, dend_icl_synapses[j], color='orange', alpha=1)
        ax_[1, j].plot(t, dend_icl_synapses[j+4], color='orange', alpha=1)
        ax_[2, j].plot(t, dend_icl_synapses[j+8], color='orange', alpha=1)
        ax_[3, j].plot(t, dend_icl_synapses[j+12], color='orange', alpha=1)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Courant\n[pA]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_isyn, max_isyn+0.5)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Current\n[pA]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_isyn, max_isyn+0.5)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)


    e = np.zeros_like(dend_igaba[0])
    plt.figure("Isynapse_label")
    plt.plot(t, e, color='black', linestyle='--', label=r'$I_{tot,synapse}$')
    plt.plot(t, e, color='deepskyblue', label=r'$I_{hco3}$')
    plt.plot(t, e, color='orange', label=r'$I_{cl,synapse}$')
    plt.ylim(1,2)
    plt.legend(loc='center')

del dend_icl
del soma_icl
del dend_ik
del soma_ik
del dend_ina
del soma_ina

if gcl_ghco3_grel == 1:
    max_g = max([np.max(gcl*1000), np.max(ghco3*1000)])
    min_g = min([np.min(gcl*1000), np.min(ghco3*1000)])

    fig_, ax_ = plt.subplots(4, 4)
    fig_.canvas.manager.set_window_title('g_cl and ghco3 at synapses')
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, gcl[j]*1000, color='orange', label=r'$g_{cl}$')
        ax_[1, j].plot(t, gcl[j+4]*1000, color='orange', label=r'$g_{cl}$')
        ax_[2, j].plot(t, gcl[j+8]*1000, color='orange', label=r'$g_{cl}$')
        ax_[3, j].plot(t, gcl[j+12]*1000, color='orange', label=r'$g_{cl}$')

        ax_[0, j].plot(t, ghco3[j]*1000, color='deepskyblue', alpha=0.7, label=r'$g_{hco3}$')
        ax_[1, j].plot(t, ghco3[j+4]*1000, color='deepskyblue', alpha=0.7, label=r'$g_{hco3}$')
        ax_[2, j].plot(t, ghco3[j+8]*1000, color='deepskyblue', alpha=0.7, label=r'$g_{hco3}$')
        ax_[3, j].plot(t, ghco3[j+12]*1000, color='deepskyblue', alpha=0.7, label=r'$g_{hco3}$')

        #ax_[0, j].plot(t, grel[j], color='darkgreen', alpha=0.7)
        #ax_[1, j].plot(t, grel[j+4], color='darkgreen', alpha=0.7)
        #ax_[2, j].plot(t, grel[j+8], color='darkgreen', alpha=0.7)
        #ax_[3, j].plot(t, grel[j+12], color='darkgreen', alpha=0.7)

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel(r'$g$ [nS]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_g, max_g)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
                ax_[i,j].legend(fontsize=12)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel(r'$g$ [nS]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(min_g, max_g)
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
                ax_[i,j].legend(fontsize=12)

del gcl
del ghco3
del grel

if open_states == 1:
    plt.figure('open_states_o1')
    for i, val in enumerate(o1*rnum):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(o1)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            plt.plot(t, val, color=colormap_syn_pos[i])
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')
    if graph_fr:
        plt.ylabel("Nombre de canaux\ndans l'état ouvert 1 [-]")
        plt.xlabel("Temps [s]")
    else:
        plt.ylabel("Number of channels\nin the open state 1 [-]")
        plt.xlabel("Time [s]")


    plt.figure('open_states_o2')
    for i, val in enumerate(o2*rnum):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(o2)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            plt.plot(t, val, color=colormap_syn_pos[i])
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')
    if graph_fr:
        plt.ylabel("Nombre de canaux\ndans l'état ouvert 2 [-]")
        plt.xlabel("Temps [s]")
    else:
        plt.ylabel("Number of channels\nin the open state 2 [-]")
        plt.xlabel("Time [s]")


    plt.figure('open_states_o3')
    for i, val in enumerate(o3*rnum):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(o3)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            plt.plot(t, val, color=colormap_syn_pos[i])
    plt.xlim(xlimit_min, xlimit_max)
    plt.legend(loc='upper right')
    if graph_fr:
        plt.ylabel("Nombre de canaux\ndans l'état ouvert 3 [-]")
        plt.xlabel("Temps [s]")
    else:
        plt.ylabel("Number of channels\nin the open state 3 [-]")
        plt.xlabel("Time [s]")


    plt.figure('open_states_total')
    for i, val in enumerate((o1+o2+o3)*rnum):
        if i == 0:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle=':', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle=':', zorder=3)
        elif i == len(o1)-1:
            if graph_fr:
                plt.plot(t, val, color='black', label=f"à {syn_pos[i]*lenght:.1f} µm dans la dendrite", linestyle='--', zorder=3)
            else:
                plt.plot(t, val, color='black', label=f"at {syn_pos[i]*lenght:.1f} µm in dendrite", linestyle='--', zorder=3)
        else:
            plt.plot(t, val, color=colormap_syn_pos[i])
    plt.xlim(xlimit_min, xlimit_max)
    if graph_fr:
        plt.ylabel("Nombre de canaux ouverts [-]")
        plt.xlabel("Temps [s]")
    else:
        plt.ylabel("Number of opened channels [-]")
        plt.xlabel("Time [s]")
    plt.legend(loc='upper right')

if open_states_all == 1:
    op = (o1+o2+o3)*rnum

    fig_, ax_ = plt.subplots(4, 4)
    for j in range(4):
        ax_[0,j].set_title(f'{syn_pos[j]*lenght:.1f} µm', fontsize=12)
        ax_[1,j].set_title(f'{syn_pos[j+4]*lenght:.1f} µm', fontsize=12)
        ax_[2,j].set_title(f'{syn_pos[j+8]*lenght:.1f} µm', fontsize=12)
        ax_[3,j].set_title(f'{syn_pos[j+12]*lenght:.1f} µm', fontsize=12)

        ax_[0, j].plot(t, op[j], color='black')
        ax_[1, j].plot(t, op[j+4], color='black')
        ax_[2, j].plot(t, op[j+8], color='black')
        ax_[3, j].plot(t, op[j+12], color='black')

    if graph_fr:
        for i in range(4):
            ax_[i, 0].set_ylabel('Canaux\nouverts\n[-]')
            ax_[3, i].set_xlabel('Temps [s]')
            for j in range(4):
                ax_[i,j].set_ylim(op.min(), op.max())
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    else:
        for i in range(4):
            ax_[i, 0].set_ylabel('Opened\nchannels\n[-]')
            ax_[3, i].set_xlabel('Time [s]')
            for j in range(4):
                ax_[i,j].set_ylim(op.min(), op.max())
                ax_[i,j].set_xlim(xlimit_min, xlimit_max)
                ax_[i,j].tick_params(axis='both', labelsize=18)
    
    del op

del o1
del o2
del o3

plt.show()