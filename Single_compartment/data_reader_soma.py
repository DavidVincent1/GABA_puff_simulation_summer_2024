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
path = r"Test_soma\unclamped_somasim_simlen=650000_dt=(5,0.025,0.025)_test.hdf5"
f = h5py.File(path, 'r')


# To see the keys of the datatset
#print(list(f.keys()))


decal = 590000   # Offset to skip the initial stabilization of the simulation
graph_fr = False # If True, the graphs axis, titles and legends will be in french


# Graphs choices. Put 1 if you want the graph and 0 if not. -----------------------------------------------------------------------
show_info = 0 # Print information on the simulation

important = 1 # Intracellular concentrations in soma, extracellular potassium concentration, membrane potential, reversal potentials

concentration = 0  # Intracellular concentrations in soma

mp_soma = 0          # 2 graphs : Membrane potential in soma and at one point in dendrite. Membrane potentials at each synapse.
current_soma = 0     # All currents in soma and at one point in dendrite
rev_pot_soma = 0     # Reversal potentials in soma and at one point in dendrite

icl_dend_soma_add_check = 1 # Chloride currents in the soma and at one point in the dendrite
ik_dend_soma_add_check = 1  # Potassium currents in the soma and at one point in the dendrite
ina_dend_soma_add_check = 1 # Sodium currents in the soma and at one point in the dendrite


# Separation of the dataset in arrays ----------------------------------------------------
t_and_soma_v = f["mp_soma"] # Time and MP in soma
#print(len(t_and_soma_v[0]))

soma_conc = f["conc_soma"]             # Concentrations in soma
soma_e = f["e_soma"]                   # Reversal potentials in soma
soma_current_cl = f["soma_current_cl"] # Icl in soma
soma_current_k = f["soma_current_k"]   # Ik in soma
soma_current_na = f["soma_current_na"] # Ina in soma


# Other values used later -----------------------------------------------------
sim_time = soma_conc.attrs["simulation_length"]        # Simulation lenght [ms]


# Everything separated -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
t_ = t_and_soma_v[0]  # Time array [ms]
t = (t_ - decal)/1000 # Time array with shifted zero [s]

soma_v = t_and_soma_v[1] # Membrane potential in soma

# Concentration arrays in dendrite and in soma
soma_cli, soma_ki, soma_nai, soma_ko = soma_conc[0], soma_conc[1], soma_conc[2], soma_conc[3]
print(soma_cli[-1])
print(soma_ki[-1])
print(soma_nai[-1])

# Chloride currents arrays in dendrite and in soma
soma_icl, soma_icl_kcc2, soma_icl_nkcc1, soma_icl_leak, soma_icl_clc2 = soma_current_cl[0], soma_current_cl[1], soma_current_cl[2], soma_current_cl[3], soma_current_cl[4]
print(soma_icl_leak[-1])
print(soma_icl_clc2[-1])

# Potassium currents arrays in dendrite and in soma
soma_ik, soma_ik_kcc2, soma_ik_nkcc1, soma_ik_leak, soma_ik_nak, soma_ik_hh = soma_current_k[0], soma_current_k[1], soma_current_k[2], soma_current_k[3], soma_current_k[4], soma_current_k[5]

# Sodium currents arrays in dendrite and in soma
soma_ina, soma_ina_nkcc1, soma_ina_leak, soma_ina_nak, soma_ina_hh = soma_current_na[0], soma_current_na[1], soma_current_na[2], soma_current_na[3], soma_current_na[4]

# Reversal potentials arrays in dendrite and in soma
soma_ecl, soma_ek, soma_ena = soma_e[0], soma_e[1], soma_e[2]


# Creation of the colormap and x-axis limits for the graphs------
xlimit_min = 0
xlimit_max = (sim_time-decal)/1000


# Deletion of the arrays that won't be used again
del t_and_soma_v
del soma_conc
del soma_e
del soma_current_cl
del soma_current_k
del soma_current_na



if important == 1:
    fig, ax = plt.subplots(2,3)
    ax[0,0].plot(t, soma_ko, color='green')
    ax[0,0].set_xlabel("Time (s)")
    ax[0,0].set_ylabel(r"$[K^+]_o$ (mM)")
    ax[0,0].set_xlim(xlimit_min, xlimit_max)
    ax[0,0].set_ylim(min(soma_ko[int(decal/5):])-1, max(soma_ko[int(decal/5):])+1)

    ax[1,0].plot(t, soma_ki, color="darkgreen")
    ax[1,0].set_xlabel("Time (s)")
    ax[1,0].set_ylabel(r"$[K^+]_i$ (mM)")
    ax[1,0].set_xlim(xlimit_min, xlimit_max)
    ax[1,0].set_ylim(min(soma_ki[int(decal/5):])-1, max(soma_ki[int(decal/5):])+1)

    ax[0,1].plot(t, soma_cli, color='orange')
    ax[0,1].set_xlabel("Time (s)")
    ax[0,1].set_ylabel(r"$[Cl^-]_i$ (mM)")
    ax[0,1].set_xlim(xlimit_min, xlimit_max)
    ax[0,1].set_ylim(min(soma_cli[int(decal/5):])-1, max(soma_cli[int(decal/5):])+1)

    ax[1,1].plot(t, soma_nai, color='blue')
    ax[1,1].set_xlabel("Time (s)")
    ax[1,1].set_ylabel(r"$[Na^+]_i$ (mM)")
    ax[1,1].set_xlim(xlimit_min, xlimit_max)
    ax[1,1].set_ylim(min(soma_nai[int(decal/5):])-1, max(soma_nai[int(decal/5):])+1)

    ax[0,2].plot(t, soma_v, color="black")
    ax[0,2].set_xlabel("Time (s)")
    ax[0,2].set_ylabel(r"MP (mV)")
    ax[0,2].set_xlim(xlimit_min, xlimit_max)
    ax[0,2].set_ylim(min(soma_v[int(decal/5):])-1, max(soma_v[int(decal/5):])+1)

    ax[1,2].plot(t, soma_v, color="black", label='MP')
    ax[1,2].plot(t, soma_ecl, color="orange", label=r'$E_{cl}$')
    ax[1,2].plot(t, soma_ek, color='darkgreen', label=r'$E_{k}$')
    ax[1,2].plot(t, soma_ena, color='blue', label=r'$E_{na}$')
    ax[1,2].set_xlabel("Time (s)")
    ax[1,2].set_ylabel(r"Potential (mV)")
    ax[1,2].set_xlim(xlimit_min, xlimit_max)
    ax[1,2].set_ylim(min([min(soma_v[int(decal/5):]), min(soma_ecl[int(decal/5):]), min(soma_ek[int(decal/5):]), min(soma_ena[int(decal/5):])])-1,
                    max([max(soma_v[int(decal/5):]), max(soma_ecl[int(decal/5):]), max(soma_ek[int(decal/5):]), max(soma_ena[int(decal/5):])])+1)
    ax[1,2].legend(fontsize=16)




if concentration == 1:
    fig, ax = plt.subplots(2,2)
    ax[0,0].plot(t, soma_ko, color='black')
    ax[0,1].plot(t, soma_cli, color='orange')
    ax[1,0].plot(t, soma_nai, color='blue')
    ax[1,1].plot(t, soma_ki, color='green')

    ax[0,0].set_xlabel("Time (s)")
    ax[0,1].set_xlabel("Time (s)")
    ax[1,0].set_xlabel("Time (s)")
    ax[1,1].set_xlabel("Time (s)")

    ax[0,0].set_ylabel(r"$[K^+]_o$ (mM)")
    ax[0,1].set_ylabel(r"$[Cl^-]_i$ (mM)")
    ax[1,0].set_ylabel(r"$[Na^+]_i$ (mM)")
    ax[1,1].set_ylabel(r"$[K^+]_i$ (mM)")

    ax[0,0].set_xlim(xlimit_min, xlimit_max)
    ax[0,1].set_xlim(xlimit_min, xlimit_max)
    ax[1,0].set_xlim(xlimit_min, xlimit_max)
    ax[1,1].set_xlim(xlimit_min, xlimit_max)

    ax[0,0].set_ylim(min(soma_ko[int(decal/5):])-1, max(soma_ko[int(decal/5):])+1)
    ax[0,1].set_ylim(min(soma_cli[int(decal/5):])-1, max(soma_cli[int(decal/5):])+1)
    ax[1,0].set_ylim(min(soma_nai[int(decal/5):])-1, max(soma_nai[int(decal/5):])+1)
    ax[1,1].set_ylim(min(soma_ki[int(decal/5):])-1, max(soma_ki[int(decal/5):])+1)

del soma_cli, soma_ki, soma_ko, soma_nai


if mp_soma == 1:
    plt.figure('MP_one_curve')
    if graph_fr:
        plt.plot(t, soma_v)
        plt.xlabel('Temps [ms]')
        plt.ylabel('Potential de membrane [mV]')
    else:
        plt.plot(t, soma_v)
        plt.xlabel('Time [ms]')
        plt.ylabel('Membrane potential [mV]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.ylim(min(soma_v[int(decal/5):])-1, max(soma_v[int(decal/5):])+1)

if current_soma == 1:
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
    plt.ylim(min([min(soma_icl[int(decal/5):]), min(soma_ik[int(decal/5):]), min(soma_ina[int(decal/5):])])-1,
            max([max(soma_icl[int(decal/5):]), max(soma_ik[int(decal/5):]), max(soma_ina[int(decal/5):])])+1)
    plt.legend()

if rev_pot_soma == 1:
    if graph_fr:
        plt.figure('Reversal_potentials')
        plt.plot(t, soma_ecl, label=r'$E_{cl}$ dans le soma', color="red", linestyle='--')
        plt.plot(t, soma_ek, label=r'$E_{k}$ dans le soma', color="darkgreen", linestyle='--')
        plt.plot(t, soma_ena, label=r'$E_{na}$ dans le soma', color="deepskyblue", linestyle='--')
        plt.plot(t, soma_v, label=r'$MP$ dans le soma', color="orange", linestyle='--')

        plt.xlabel('Temps [s]')
        plt.ylabel('Potentiel réversible [mV]')
    else:
        plt.figure('Reversal_potentials')
        plt.plot(t, soma_ecl, label=r'$E_{cl}$ in soma', color="red", linestyle='--')
        plt.plot(t, soma_ek, label=r'$E_{k}$ in soma', color="darkgreen", linestyle='--')
        plt.plot(t, soma_ena, label=r'$E_{na}$ in soma', color="deepskyblue", linestyle='--')
        plt.plot(t, soma_v, label=r'$MP$ in soma', color="orange", linestyle='--')

        plt.xlabel('Time [s]')
        plt.ylabel('Reversal potential [mV]')
    plt.xlim(xlimit_min, xlimit_max)
    plt.ylim(min([min(soma_ecl[int(decal/5):]), min(soma_ek[int(decal/5):]), min(soma_ena[int(decal/5):]), min(soma_v[int(decal/5):])])-1,
            max([max(soma_ecl[int(decal/5):]), max(soma_ek[int(decal/5):]), max(soma_ena[int(decal/5):]), max(soma_v[int(decal/5):])])+1)
    plt.legend(loc='upper right')

del soma_v
del soma_ecl
del soma_ek
del soma_ena


if icl_dend_soma_add_check == 1:
    max_icl_soma = max([max(soma_icl[int(decal/5):]),
                        max(soma_icl_clc2[int(decal/5):]), max(soma_icl_kcc2[int(decal/5):]),
                        max(soma_icl_nkcc1[int(decal/5):]), max(soma_icl_leak[int(decal/5):])])+1
    min_icl_soma = min([min(soma_icl[int(decal/5):]),
                        min(soma_icl_clc2[int(decal/5):]), min(soma_icl_kcc2[int(decal/5):]),
                        min(soma_icl_nkcc1[int(decal/5):]), min(soma_icl_leak[int(decal/5):])])-1

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
    plt.plot(t, soma_icl_clc2, label=r'$I_{cl,clc2}$', color="darkorchid", linestyle=':', lw=2, alpha=1)
    plt.plot(t, soma_icl_leak + soma_icl_nkcc1 + soma_icl_kcc2 + soma_icl_clc2,
            label=r'$I_{cl,kcc2}+I_{cl,nkcc1}+I_{cl,leak}+I_{cl,clc2}$', color="black", linestyle='--', alpha=0.3, lw=1.5)
    plt.xlim(xlimit_min, xlimit_max)
    plt.ylim(min_icl_soma, max_icl_soma)
    plt.legend(loc='upper right')

if ik_dend_soma_add_check == 1:
    max_ik_soma = max([max(soma_ik[int(decal/5):]), max(soma_ik_hh[int(decal/5):]),
                        max(soma_ik_nak[int(decal/5):]), max(soma_ik_kcc2[int(decal/5):]),
                        max(soma_ik_nkcc1[int(decal/5):]), max(soma_ik_leak[int(decal/5):])])+1
    min_ik_soma = min([min(soma_ik[int(decal/5):]), min(soma_ik_hh[int(decal/5):]),
                        min(soma_ik_nak[int(decal/5):]), min(soma_ik_kcc2[int(decal/5):]),
                        min(soma_ik_nkcc1[int(decal/5):]), min(soma_ik_leak[int(decal/5):])])-1
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
    plt.ylim(min_ik_soma, max_ik_soma)
    plt.legend(loc='upper right')

if ina_dend_soma_add_check == 1:
    max_ina_soma = max([max(soma_ina[int(decal/5):]), max(soma_ina_hh[int(decal/5):]),
                        max(soma_ina_nak[int(decal/5):]),
                        max(soma_ina_nkcc1[int(decal/5):]), max(soma_ina_leak[int(decal/5):])])+1
    min_ina_soma = min([min(soma_ina[int(decal/5):]), min(soma_ina_hh[int(decal/5):]),
                        min(soma_ina_nak[int(decal/5):]),
                        min(soma_ina_nkcc1[int(decal/5):]), min(soma_ina_leak[int(decal/5):])])-1
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

del soma_icl
del soma_ik
del soma_ina

plt.show()