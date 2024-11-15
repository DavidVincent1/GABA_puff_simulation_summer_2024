import time
texp = time.time()

import numpy as np
import matplotlib.pyplot as plt
from function import Soma_leak_on_clc2, syna_kcc2_nkcc1_soma
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
save_ = True # True or False


# temperature in celsius -
h.celsius = p.temperature


# KCC2 and NKCC1 strenght lists --------------------------------------------------------
# INPUT 2
nb_ratio = 50  # the 2D grid will be nb_kcc2 x nb_nkcc1

# INPUT 3
max_conductance = 1.5e-5

gclc2 = np.linspace(1.35e-6, max_conductance, nb_ratio)
gcl_leak = max_conductance - gclc2
ratio = gcl_leak/gclc2


# Arrays for results ---------------------------------
cells_delta_chloride   = np.zeros(nb_ratio)
cells_stable_chloride  = np.zeros(nb_ratio)
cells_stable_potassium = np.zeros(nb_ratio)
cells_stable_sodium    = np.zeros(nb_ratio)
cells_stable_mp        = np.zeros(nb_ratio)
cells_stable_ecl       = np.zeros(nb_ratio)
cells_stable_icl       = np.zeros(nb_ratio)
cells_stable_icl_kcc2  = np.zeros(nb_ratio)
cells_stable_icl_nkcc1 = np.zeros(nb_ratio)
cells_stable_icl_clc2  = np.zeros(nb_ratio)
cells_stable_icl_leak  = np.zeros(nb_ratio)


# Function for the simulations -----------------------------------------------------------------
def run_simulation(i):
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
    MP = -65

    # Creation of the cell
    my_cell = Soma_leak_on_clc2(0,
                    number_of_soma_segments=p.soma_nseg,
                    cli_0=p.intial_cli,
                    nai_0=p.initial_nai,
                    ki_0=p.initial_ki,
                    clo_0=p.clo,
                    nao_0=p.nao,
                    ko_0=p.ko,
                    ukcc2=p.U_kcc2,
                    unkcc1=p.U_nkcc1,
                    gcl_leak=gcl_leak[i],
                    gclc2=gclc2[i])

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
    results = (i, Potassic_event[0], Potassic_event[1],
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
    (i, delta_cl, stable_cl, stable_k, stable_na, stable_mp, stable_ecl, cli, ki, nai, v, t,
    stable_icl, stable_icl_kcc2, stable_icl_nkcc1, stable_icl_clc2, stable_icl_leak) = result

    cells_delta_chloride[i] = delta_cl
    cells_stable_chloride[i] = stable_cl
    cells_stable_potassium[i] = stable_k
    cells_stable_sodium[i] = stable_na
    cells_stable_mp[i] = stable_mp
    cells_stable_ecl[i] = stable_ecl
    cells_stable_icl[i] = stable_icl
    cells_stable_icl_kcc2[i] = stable_icl_kcc2
    cells_stable_icl_nkcc1[i] = stable_icl_nkcc1
    cells_stable_icl_clc2[i] = stable_icl_clc2
    cells_stable_icl_leak[i] = stable_icl_leak

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
        for i in tqdm(range(nb_ratio)):
            pool.apply_async(run_simulation, args=(i,), callback=collect_result)
        pool.close()
        pool.join()

    tcomp = time.time() - texp
    print("Simulation done")
    print("Computational time [s] : ", tcomp)
    print("Computational time [min, s] : ", tcomp//60, tcomp-(tcomp//60)*60)


    # Graphs -------------------------------------------------------------------------------------------
    fig1 = plt.figure("Chloride_rest_concentration")
    plt.plot(ratio, cells_stable_chloride, color='black')
    plt.xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    plt.ylabel(r"$[Cl^-]_i$ (mM)")


    fig2 = plt.figure("Membrane_rest_potential")
    plt.plot(ratio, cells_stable_mp, color='blue')
    plt.xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    plt.ylabel(r"Membrane potential (mV)")


    fig3 = plt.figure("Chloride_current")
    plt.plot(ratio, cells_stable_icl, color='orange', label=r'$I_{cl,tot}$')
    plt.plot(ratio, cells_stable_icl_kcc2, color='mediumblue', label=r'$I_{cl,kcc2}$')
    plt.plot(ratio, cells_stable_icl_nkcc1, color='darkkhaki', label=r'$I_{cl,nkcc1}$')
    plt.plot(ratio, cells_stable_icl_clc2, color='darkorchid', label=r'$I_{cl,clc2}$')
    plt.plot(ratio, cells_stable_icl_leak, color='chocolate', label=r'$I_{cl,leak}$')
    plt.plot(ratio, cells_stable_icl_kcc2+cells_stable_icl_nkcc1+cells_stable_icl_clc2+cells_stable_icl_leak,
            color='black', linestyle='--', label=r'$I_{cl,add check}$')
    plt.xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    plt.ylabel(r"Current (mA)")
    plt.legend(fontsize=18)

    fig4, ax = plt.subplots(1,3)
    ax[0].plot(ratio, cells_stable_chloride, color='black')
    ax[0].set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax[0].set_ylabel(r"$[Cl^-]_i$ (mM)")

    ax[1].plot(ratio, cells_stable_mp, color='blue')
    ax[1].set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax[1].set_ylabel(r"Membrane potential (mV)")

    ax[2].plot(ratio, cells_stable_icl, color='orange', label=r'$I_{cl,tot}$')
    ax[2].plot(ratio, cells_stable_icl_kcc2, color='mediumblue', label=r'$I_{cl,kcc2}$')
    ax[2].plot(ratio, cells_stable_icl_nkcc1, color='darkkhaki', label=r'$I_{cl,nkcc1}$')
    ax[2].plot(ratio, cells_stable_icl_clc2, color='darkorchid', label=r'$I_{cl,clc2}$')
    ax[2].plot(ratio, cells_stable_icl_leak, color='chocolate', label=r'$I_{cl,leak}$')
    ax[2].plot(ratio, cells_stable_icl_kcc2+cells_stable_icl_nkcc1+cells_stable_icl_clc2+cells_stable_icl_leak,
            color='black', linestyle='--', label=r'$I_{cl,add check}$')
    ax[2].set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax[2].set_ylabel(r"Current (mA)")
    ax[2].legend(fontsize=18)

    # Création de la figure et de la grille
    fig4 = plt.figure(figsize=(10, 8))
    grid = fig4.add_gridspec(2, 2)

    # Premier panneau : première ligne, première colonne
    ax1 = fig4.add_subplot(grid[0, 0])
    ax1.plot(ratio, cells_stable_chloride, color='black')
    ax1.set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax1.set_ylabel(r"$[Cl^-]_i$ (mM)")

    # Deuxième panneau : première ligne, deuxième colonne
    ax2 = fig4.add_subplot(grid[0, 1])
    ax2.plot(ratio, cells_stable_mp, color='blue')
    ax2.set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax2.set_ylabel(r"MP (mV)")

    # Troisième panneau : deuxième ligne, fusion des deux colonnes
    ax3 = fig4.add_subplot(grid[1, :])  # Fusionne les deux colonnes
    ax3.plot(ratio, cells_stable_icl, color='orange', label=r'$I_{cl,tot}$')
    ax3.plot(ratio, cells_stable_icl_kcc2, color='mediumblue', label=r'$I_{cl,kcc2}$')
    ax3.plot(ratio, cells_stable_icl_nkcc1, color='darkkhaki', label=r'$I_{cl,nkcc1}$')
    ax3.plot(ratio, cells_stable_icl_clc2, color='darkorchid', label=r'$I_{cl,clc2}$')
    ax3.plot(ratio, cells_stable_icl_leak, color='chocolate', label=r'$I_{cl,leak}$')
    ax3.plot(ratio, cells_stable_icl_kcc2+cells_stable_icl_nkcc1+cells_stable_icl_clc2+cells_stable_icl_leak,
            color='black', linestyle='--', label=r'$I_{cl,add check}$')
    ax3.set_xlabel(r"$g_{cl,leak}/g_{clc2}$ (-)")
    ax3.set_ylabel(r"Current (mA)")
    ax3.legend(fontsize=18)


    if save_:
        # To save the graphs as pickle files.
        if p.clamp:
            folder_path = filepath + f'\pickle_{nb_ratio}_gclc2={gclc2[0]:.2g}-{gclc2[-1]:.2g}_gleak={gcl_leak[0]:.2g}-{gcl_leak[-1]:.2g}_Ukcc2={(p.U_kcc2*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}_Unkcc1={(p.U_nkcc1*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}_clamp={p.clamp_amp}'
        else:
            folder_path = filepath + f'\pickle_{nb_ratio}_gclc2={gclc2[0]:.2g}-{gclc2[-1]:.2g}_gleak={gcl_leak[0]:.2g}-{gcl_leak[-1]:.2g}_Ukcc2={(p.U_kcc2*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}_Unkcc1={(p.U_nkcc1*p.F*p.soma_vol)/(p.soma_surf*1e4):.2g}_unclamp'
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
        else:
            input = input('The folder for the pickle files already exists. Proceed anyway ? Yes/No : ')
            if input == 'No' or input == 'no' or input == 'Non' or input == 'non':
                raise Exception('Manually stoped')

        pickle.dump(fig1, open(folder_path + '\Chloride_rest_concentration.pickle', 'wb'))
        pickle.dump(fig2, open(folder_path + '\Membrane_rest_potential.pickle', 'wb'))
        pickle.dump(fig3, open(folder_path + '\Chloride_current.pickle', 'wb'))
        pickle.dump(fig4, open(folder_path + '\Three_panels.pickle', 'wb'))


plt.show()