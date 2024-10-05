import time
texp = time.time()

import numpy as np
import matplotlib.pyplot as plt
from function import NeuronCell, syna_kcc2_nkcc1
from tqdm import tqdm
from neuron import h
from matplotlib import rcParams
from multiprocessing import Pool
import os
import pickle
import parameters as p
rcParams.update({'font.size': 22}) # Graph parameters


# Loading mod files and for better simulation --
h.nrn_load_dll(r"mod_files\nrnmech.dll")
h.load_file("stdrun.hoc")


# INPUT 1
# path where to save the pickle files. It makes it possible to edit the graphs later throught the
# 'data_reader_pickle.py' file.
filepath = r"dataset" 


# temperature in celsius -
h.celsius = p.temperature


# KCC2 and NKCC1 strenght lists ------------------------------------------------------
# INPUT 2
nb_kcc2 = 5  # the 2D grid will be nb_kcc2 x nb_nkcc1
nb_nkcc1 = 5 # ^

# INPUT 3
range_kcc2 = 1e-6, 5e-4  # Range will begin at range_kcc2[0] and end at range_kcc2
range_nkcc1 = 1e-6, 5e-4 # Range will begin at range_nkcc1[0] and end at range_nkcc1

Ukcc2 = np.linspace(range_kcc2[1], range_kcc2[0], nb_kcc2)   # Range of U_KCC2 values
Unkcc1 = np.linspace(range_nkcc1[0], range_nkcc1[1], nb_nkcc1) # Range of U_nkc1 values


# Arrays for results ---------------------------------
cells_delta_chloride = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_chloride = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_potassium = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_sodium = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_mp = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_ecl = np.zeros((nb_kcc2, nb_nkcc1))


# Number of synapses and synapses positions -------------------------------------------------------------
syn_nb = int(p.dend_lenght_with_synapses * p.syn_per_micron)                       # number of synapses
synapse_pos = np.linspace(0.01, p.dend_lenght_with_synapses/p.dend_lenght, syn_nb) # position of synapses


# Recording vectors positions for concentration (K+, Na+, Cl-) ----------
record_pos1 = np.linspace(0.01, (p.position_of_puff*2)/p.dend_lenght, 40)


# Assure that there is a recording vector at the GABA puff position
mini, mini_pos = record_pos1[0], 0
for i,j in enumerate(record_pos1):
    if abs(j-p.position_of_puff/p.dend_lenght) < abs(mini-p.position_of_puff/p.dend_lenght):
        mini = j
        mini_pos = i
record_pos1[mini_pos] = p.position_of_puff/p.dend_lenght


# GABA puff position and index ----------------------------------------------------------------------
puff_pos = p.position_of_puff/p.dend_lenght # position of the puff event on the dendrite (first part)
puff_indice = int((p.time_of_puff-100)/p.dt1 + (p.time_of_puff-p.time_for_stabilization)/p.dt2)


# Function for the simulations -----------------------------------------------------------------
def run_simulation(i, j):
    """This function create a cell and run a simulation. It is essential to use multiprocessing.

    Args:
        i (int): First index of the grid. Corresponds to a certain U_kcc2 value.
        j (int): Seconde index of the grid. Corresponds to a certain U_nkcc1 value.

    Returns:
        Tuple : Contain 13 elements.
                Tuple[0] : i
                Tuple[1] : j
                Tuple[2] : Maximum variation of intracellular chloride due to GABA stimulation
                Tuple[3] : Intracellular chloride concentration just before GABA stimulation
                Tuple[4] : Intracellular potassium concentration just before GABA stimulation
                Tuple[5] : Intracellular sodium concentration just before GABA stimulation
                Tuple[6] : Membrane potential just before GABA stimulation
                Tuple[7] : Chloride reversal potential just before GABA stimulation
                Tuple[8] : Sodium
    """
    # Initial membrane potential
    MP = -100

    # Creation of the cell
    my_cell = NeuronCell(0,
                            number_of_dendrite_segments=p.dend_nseg,
                            number_of_soma_segments=p.soma_nseg,
                            number_of_dendrite_segments2=p.dend2_nseg,
                            dendrite_length_um=p.dend_lenght,
                            dendrite_length_um2=p.dend2_lenght,
                            cli_0=p.intial_cli,
                            nai_0=p.initial_nai,
                            ki_0=p.initial_ki,
                            clo_0=p.clo,
                            nao_0=p.nao,
                            ko_0=p.ko,
                            ukcc2=Ukcc2[i],
                            unkcc1=Unkcc1[j])

    # GABA stimulation event
    GABA_puff_event = syna_kcc2_nkcc1(first=puff_indice,
                            cell=my_cell,
                            pos=synapse_pos,
                            puff_pos=puff_pos,
                            puff_time=p.time_of_puff,
                            puff_conc=p.concentration_of_puff,
                            tau=p.tau_GABA,
                            dgab=p.Dgaba,
                            rnum=[p.rnum],
                            clamp=p.clamp,
                            clamp_amp=p.clamp_amp,
                            rmp_initial=MP,
                            sim_time=p.simulation_lenght,
                            record_pos=record_pos1,
                            skip=p.time_for_stabilization,
                            pipette=p.pipett)

    # Tuple returned by the function
    results = (i, j, GABA_puff_event[0], GABA_puff_event[1],
                GABA_puff_event[2], GABA_puff_event[3], GABA_puff_event[4],
                GABA_puff_event[5], GABA_puff_event[6], GABA_puff_event[7],
                GABA_puff_event[8], GABA_puff_event[9], GABA_puff_event[10])
    del my_cell
    del GABA_puff_event
    return results

# Function to collect the results --------------------------------------------------------------------------
def collect_result(result):
    """This function is there to collect the resultas of run_simulation function when using multiprocessing. 

    Args:
        result (Tuple): Contain 13 elements (from run_simulation function).
                Tuple[0] : i
                Tuple[1] : j
                Tuple[2] : Maximum variation of intracellular chloride due to GABA stimulation
                Tuple[3] : Intracellular chloride concentration just before GABA stimulation
                Tuple[4] : Intracellular potassium concentration just before GABA stimulation
                Tuple[5] : Intracellular sodium concentration just before GABA stimulation
                Tuple[6] : Membrane potential just before GABA stimulation
                Tuple[7] : Chloride reversal potential just before GABA stimulation
                Tuple[8] : Sodium
    """
    i, j, delta_cl, stable_cl, stable_k, stable_na, stable_mp, stable_ecl, nai, ki, cli, v, t = result
    cells_delta_chloride[i, j] = delta_cl
    cells_stable_chloride[i, j] = stable_cl
    cells_stable_potassium[i, j] = stable_k
    cells_stable_sodium[i, j] = stable_na
    cells_stable_mp[i, j] = stable_mp
    cells_stable_ecl[i, j] = stable_ecl

    # Graphs of all the curves from all the simulations
    plt.figure("Sodium steady state")
    plt.plot(t, nai)
    plt.figure("Potassium steady state")
    plt.plot(t, ki)
    plt.figure("Chloride steady state")
    plt.plot(t, cli)
    plt.figure("Resting mp")
    plt.plot(t, v)


# if __name__ bloc necessary for multiprocessing
if __name__ == '__main__':
    # Graphs of all the curves from all the simulations
    fig1 = plt.figure("Sodium steady state")
    plt.ylabel(r"$[Na^+]_i$ in soma [mM]")
    plt.xlabel("Time [ms]")
    fig2 = plt.figure("Potassium steady state")
    plt.ylabel(r"$[K^+]_i$ in soma [mM]")
    plt.xlabel("Time [ms]")
    fig3 = plt.figure("Chloride steady state")
    plt.ylabel(r"$[Cl^-]_i$ in soma [mM]")
    plt.xlabel("Time [ms]")
    fig4 = plt.figure("Resting mp")
    plt.ylabel(r"$MP$ in soma [mV]")
    plt.xlabel("Time [ms]")

    # Time of the beginning of the multiprocessing
    # Used later to calculate the computational time of the code
    texp = time.time()

    # number of usable cpu 
    num_cores = os.cpu_count()

    # Multiprocessing
    with Pool(num_cores) as pool:
        for i in tqdm(range(len(Ukcc2))):
            for j in range(len(Unkcc1)):
                pool.apply_async(run_simulation, args=(i, j), callback=collect_result)
        pool.close()
        pool.join()

    tcomp = time.time() - texp
    print("Simulation done")
    print("Computational time [s] : ", tcomp)
    print("Computational time [min, s] : ", tcomp//60, tcomp-(tcomp//60)*60)


    # Graphs -------------------------------------------------------------------------------------------
    # Labels lists for the grid axes
    Unkcc1_label = [f"{Unkcc1[i]:2g}" for i in range(len(Unkcc1))]
    Ukcc2_label = [f"{Ukcc2[i]:2g}" for i in range(len(Ukcc2))]

    # Maximum chloride variation due to the GABA puff
    fig5 = plt.figure("Delta chloride")
    im = plt.imshow(cells_delta_chloride, cmap='cividis')
    cbar = plt.colorbar(im, extend='both', spacing='proportional', label=r'$(\Delta [Cl^-]_i)_{max}$ [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_delta_chloride[i, j]:.3f}', ha="center", va="center", color="w")
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # Chloride concentration just before the GABA puff
    fig6 = plt.figure("Stable chloride")
    im2 = plt.imshow(cells_stable_chloride, cmap='cividis')
    cbar2 = plt.colorbar(im2, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_chloride[i, j]:.3f}', ha="center", va="center", color="w")
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # Potassium concentration just before the GABA puff
    fig7 = plt.figure("Stable potassium")
    im3 = plt.imshow(cells_stable_potassium, cmap='cividis')
    cbar3 = plt.colorbar(im3, extend='both', spacing='proportional', label=r'$[K^+]_i$ before puff [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_potassium[i, j]:.2f}', ha="center", va="center", color="w")
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # Sodium concentration just before the GABA puff
    fig8 = plt.figure("Stable sodium")
    im4 = plt.imshow(cells_stable_sodium, cmap='cividis')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_sodium[i, j]:.3f}', ha="center", va="center", color="w")
    cbar4 = plt.colorbar(im4, extend='both', spacing='proportional', label=r'$[Na^+]_i$ before puff [mM]')
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # Membrane potential just before the GABA puff
    fig9 = plt.figure("Stable mp")
    im5 = plt.imshow(cells_stable_mp, cmap='cividis')
    cbar5 = plt.colorbar(im5, extend='both', spacing='proportional', label=r'$MP$ before puff [mV]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_mp[i, j]:.3f}', ha="center", va="center", color="w")
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # Chloride reversal potential just before the GABA puff
    fig10 = plt.figure("Stable Ecl")
    im5 = plt.imshow(cells_stable_ecl, cmap='cividis')
    cbar5 = plt.colorbar(im5, extend='both', spacing='proportional', label=r'$E_{cl}$ before puff [mV]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_ecl[i, j]:.3f}', ha="center", va="center", color="w")
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [mM/ms]")
    plt.ylabel(r"$U_{kcc2}$ [mM/ms]")

    # To save the graphs as pickle files.
    if p.clamp:
        folder_path = filepath + f'\pickle_kcc2={range_kcc2[0]}-{range_kcc2[1]}_nkcc1={range_nkcc1[0]}-{range_nkcc1[1]}_gclc2={p.soma_gclc2}_clamp={p.clamp_amp}'
    else:
        folder_path = filepath + f'\pickle_kcc2={range_kcc2[0]}-{range_kcc2[1]}_nkcc1={range_nkcc1[0]}-{range_nkcc1[1]}_gclc2={p.soma_gclc2}_unclamp_gnaother={p.dend_gnaother}'
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    else:
        input = input('The folder for the pickle files already exists. Proceed anyway ? Yes/No : ')
        if input == 'No' or input == 'no' or input == 'Non' or input == 'non':
            raise Exception('Manually stoped')

    pickle.dump(fig1, open(folder_path + '\Sodium_steady_state.pickle', 'wb'))
    pickle.dump(fig2, open(folder_path + '\Potassium_steady_state.pickle', 'wb'))
    pickle.dump(fig3, open(folder_path + '\Chlore_steady_state.pickle', 'wb'))
    pickle.dump(fig4, open(folder_path + '\Resting_mp.pickle', 'wb'))
    pickle.dump(fig5, open(folder_path + '\Delta_chloride.pickle', 'wb'))
    pickle.dump(fig6, open(folder_path + '\Stable_chloride.pickle', 'wb'))
    pickle.dump(fig7, open(folder_path + '\Stable_potassium.pickle', 'wb'))
    pickle.dump(fig8, open(folder_path + '\Stable_sodium.pickle', 'wb'))
    pickle.dump(fig9, open(folder_path + '\Stable_mp.pickle', 'wb'))
    pickle.dump(fig10, open(folder_path + '\Stable_ecl.pickle', 'wb'))


plt.show()