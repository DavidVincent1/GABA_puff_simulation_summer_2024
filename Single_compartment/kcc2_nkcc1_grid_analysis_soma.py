import time
texp = time.time()

import numpy as np
import matplotlib.pyplot as plt
from function import Soma, syna_kcc2_nkcc1_soma
from tqdm import tqdm
from neuron import h
from matplotlib import rcParams
from multiprocessing import Pool
import os
import pickle
import parameters_soma as p
rcParams.update({'font.size': 22}) # Graph parameters


# Loading mod files and for better simulation --
h.nrn_load_dll(r"mod_files\nrnmech.dll")
h.load_file("stdrun.hoc")


# INPUT 1
# path where to save the pickle files. It makes it possible to edit the graphs later throught the
# 'data_reader_pickle.py' file.
filepath = r"Single_compartment\dataset"
save_ = False # True or False


# temperature in celsius -
h.celsius = p.temperature


# KCC2 and NKCC1 strenght lists --------------------------------------------------------
# INPUT 2
nb_kcc2 = 10  # the 2D grid will be nb_kcc2 x nb_nkcc1
nb_nkcc1 = 10 # ^

# INPUT 3
range_kcc2 = 0.5e-6, 1e-4  # Range will begin at range_kcc2[0] and end at range_kcc2
range_nkcc1 = 0.5e-6, 1e-4 # Range will begin at range_nkcc1[0] and end at range_nkcc1

Ukcc2 = np.linspace(range_kcc2[1], range_kcc2[0], nb_kcc2)     # Range of U_KCC2 values
Unkcc1 = np.linspace(range_nkcc1[0], range_nkcc1[1], nb_nkcc1) # Range of U_nkc1 values
print(Ukcc2)
print(Unkcc1)


# Arrays for results ---------------------------------
cells_delta_chloride   = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_chloride  = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_potassium = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_sodium    = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_mp        = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_ecl       = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_icl       = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_icl_kcc2  = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_icl_nkcc1 = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_icl_clc2  = np.zeros((nb_kcc2, nb_nkcc1))
cells_stable_icl_leak  = np.zeros((nb_kcc2, nb_nkcc1))


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
    MP = -70

    # Creation of the cell
    my_cell = Soma(0,
                    number_of_soma_segments=p.soma_nseg,
                    cli_0=p.intial_cli,
                    nai_0=p.initial_nai,
                    ki_0=p.initial_ki,
                    clo_0=p.clo,
                    nao_0=p.nao,
                    ko_0=p.ko,
                    ukcc2=Ukcc2[i],
                    unkcc1=Unkcc1[j])

    # GABA stimulation event
    Potassic_event = syna_kcc2_nkcc1_soma(cell=my_cell,
                                        choc_time=p.tchoc,
                                        kchoc=p.kchoc,
                                        tauchoc=p.tauchoc,
                                        clamp=p.clamp,
                                        clamp_amp=p.clamp_amp,
                                        rmp_initial=MP,
                                        sim_time=p.simulation_lenght,
                                        skip=p.time_for_stabilization,
                                        pipette=p.pipett)

    # Tuple returned by the function
    results = (i, j, Potassic_event[0], Potassic_event[1],
                Potassic_event[2], Potassic_event[3], Potassic_event[4],
                Potassic_event[5], Potassic_event[6], Potassic_event[7],
                Potassic_event[8], Potassic_event[9], Potassic_event[10],
                Potassic_event[11], Potassic_event[12], Potassic_event[13],
                Potassic_event[14], Potassic_event[15])
    del my_cell
    del Potassic_event
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
    (i, j, delta_cl, stable_cl, stable_k, stable_na, stable_mp, stable_ecl, cli, ki, nai, v, t,
    stable_icl, stable_icl_kcc2, stable_icl_nkcc1, stable_icl_clc2, stable_icl_leak) = result

    cells_delta_chloride[i, j] = delta_cl
    cells_stable_chloride[i, j] = stable_cl
    cells_stable_potassium[i, j] = stable_k
    cells_stable_sodium[i, j] = stable_na
    cells_stable_mp[i, j] = stable_mp
    cells_stable_ecl[i, j] = stable_ecl
    cells_stable_icl[i, j] = stable_icl
    cells_stable_icl_kcc2[i, j] = stable_icl_kcc2
    cells_stable_icl_nkcc1[i, j] = stable_icl_nkcc1
    cells_stable_icl_clc2[i, j] = stable_icl_clc2
    cells_stable_icl_leak[i, j] = stable_icl_leak

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
    Unkcc1_label = [f"{(Unkcc1[i]*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}" for i in range(len(Unkcc1))]
    Ukcc2_label  = [f"{(Ukcc2[i]*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}" for i in range(len(Ukcc2))]

    # Maximum chloride variation due to the potassic event
    fig5 = plt.figure("Delta chloride")
    im = plt.imshow(cells_delta_chloride, cmap='cividis')
    cbar = plt.colorbar(im, extend='both', spacing='proportional', label=r'$(\Delta [Cl^-]_i)_{max}$ [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_delta_chloride[i, j]:.3f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride concentration just before the potassic event
    fig6 = plt.figure("Stable chloride")
    im2 = plt.imshow(cells_stable_chloride, cmap='cividis')
    cbar2 = plt.colorbar(im2, extend='both', spacing='proportional', label=r'$[Cl^-]_i$ before puff [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_chloride[i, j]:.3f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Potassium concentration just before the potassic event
    fig7 = plt.figure("Stable potassium")
    im3 = plt.imshow(cells_stable_potassium, cmap='cividis')
    cbar3 = plt.colorbar(im3, extend='both', spacing='proportional', label=r'$[K^+]_i$ before puff [mM]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_potassium[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Sodium concentration just before the potassic event
    fig8 = plt.figure("Stable sodium")
    im4 = plt.imshow(cells_stable_sodium, cmap='cividis')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_sodium[i, j]:.3f}', ha="center", va="center", color="w", fontsize=10)
    cbar4 = plt.colorbar(im4, extend='both', spacing='proportional', label=r'$[Na^+]_i$ before puff [mM]')
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Membrane potential just before the potassic event
    fig9 = plt.figure("Stable mp")
    im5 = plt.imshow(cells_stable_mp, cmap='cividis')
    cbar5 = plt.colorbar(im5, extend='both', spacing='proportional', label=r'$MP$ before puff [mV]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_mp[i, j]:.3f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride reversal potential just before the potassic event
    fig10 = plt.figure("Stable Ecl")
    im5 = plt.imshow(cells_stable_ecl, cmap='cividis')
    cbar5 = plt.colorbar(im5, extend='both', spacing='proportional', label=r'$E_{cl}$ before puff [mV]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_ecl[i, j]:.3f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride current just before the potassic event
    fig11 = plt.figure("Stable Icl")
    im6 = plt.imshow(cells_stable_icl, cmap='cividis')
    cbar6 = plt.colorbar(im6, extend='both', spacing='proportional', label=r'$I_{cl}$ before choc [mA]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_icl[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride KCC2 current just before the potassic event
    fig12 = plt.figure("Stable Icl kcc2")
    im7 = plt.imshow(cells_stable_icl_kcc2, cmap='cividis')
    cbar7 = plt.colorbar(im7, extend='both', spacing='proportional', label=r'$I_{cl,kcc2}$ before choc [mA]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_icl_kcc2[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride NKCC1 current just before the potassic event
    fig13 = plt.figure("Stable Icl NKCC1")
    im8 = plt.imshow(cells_stable_icl_nkcc1, cmap='cividis')
    cbar8 = plt.colorbar(im8, extend='both', spacing='proportional', label=r'$I_{cl,nkcc1}$ before choc [mA]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_icl_nkcc1[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride CLC2 current just before the potassic event
    fig14 = plt.figure("Stable Icl CLC2")
    im9 = plt.imshow(cells_stable_icl_clc2, cmap='cividis')
    cbar9 = plt.colorbar(im9, extend='both', spacing='proportional', label=r'$I_{cl,clc2}$ before choc [mA]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_icl_clc2[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    # Chloride leak current just before the potassic event
    fig15 = plt.figure("Stable Icl leak")
    im10 = plt.imshow(cells_stable_icl_leak, cmap='cividis')
    cbar10 = plt.colorbar(im10, extend='both', spacing='proportional', label=r'$I_{cl,leak}$ before choc [mA]')
    for i in range(len(Ukcc2)):
        for j in range(len(Unkcc1)):
            text = plt.text(j, i, f'{cells_stable_icl_leak[i, j]:.2f}', ha="center", va="center", color="w", fontsize=10)
    plt.xticks(np.arange(len(Unkcc1)), labels=Unkcc1_label, rotation=45, ha="right", rotation_mode="anchor")
    plt.yticks(np.arange(len(Ukcc2)), labels=Ukcc2_label)
    plt.xlabel(r"$U_{nkcc1}$ [$mA/cm^2$]")
    plt.ylabel(r"$U_{kcc2}$ [$mA/cm^2$]")

    if save_:
        # To save the graphs as pickle files.
        if p.clamp:
            folder_path = filepath + f'\pickle_{nb_kcc2}_x_{nb_nkcc1}_kcc2={range_kcc2[0]}-{range_kcc2[1]}_nkcc1={range_nkcc1[0]}-{range_nkcc1[1]}_gclc2={p.soma_gclc2:.2g}_leak={p.soma_gcl:.2g}_clamp={p.clamp_amp}'
        else:
            folder_path = filepath + f'\pickle_{nb_kcc2}_x_{nb_nkcc1}_kcc2={range_kcc2[0]}-{range_kcc2[1]}_nkcc1={range_nkcc1[0]}-{range_nkcc1[1]}_gclc2={p.soma_gclc2:.2g}_leak={p.soma_gcl:.2g}_unclamp'
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
        pickle.dump(fig11, open(folder_path + '\Stable_icl.pickle', 'wb'))
        pickle.dump(fig12, open(folder_path + '\Stable_icl_kcc2.pickle', 'wb'))
        pickle.dump(fig13, open(folder_path + '\Stable_icl_nkcc1.pickle', 'wb'))
        pickle.dump(fig14, open(folder_path + '\Stable_icl_clc2.pickle', 'wb'))
        pickle.dump(fig15, open(folder_path + '\Stable_icl_leak.pickle', 'wb'))


plt.show()