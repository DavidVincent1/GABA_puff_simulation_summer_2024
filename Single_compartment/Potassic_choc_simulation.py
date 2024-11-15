import time
texp = time.time()

from neuron import h
from neuron.units import mV
from matplotlib import rcParams
import numpy as np
import parameters_soma as p
from function import Soma, save_sim_choc_potass, choc_potass
from pathlib import Path
rcParams.update({'font.size': 22}) # Graph parameters


# Loading mod files and for better simulation --------------------------
h.nrn_load_dll(r"mod_files\nrnmech.dll")
h.load_file("stdrun.hoc")


# Temperature in celsius -
h.celsius = p.temperature


# Creation of the cell ----------------------------------
if p.clamp:
    mp = p.clamp_amp # 1500 ms pip, -90 mV
else:
    mp = -71*mV

my_cell = Soma(0,
                number_of_soma_segments=p.soma_nseg,
                cli_0=p.intial_cli,
                nai_0=p.initial_nai,
                ki_0=p.initial_ki,
                clo_0=p.clo,
                nao_0=p.nao,
                ko_0=p.ko,
                ukcc2=p.U_kcc2,
                unkcc1=p.U_nkcc1)


# Potassic choc event and simulation ------------------------------
print("Potassic choc simulation\n")
Potassic_choc_event = choc_potass(my_cell,
                                clamp=p.clamp,
                                clamp_amp=p.clamp_amp * mV,
                                rmp_initial=mp,
                                sim_time=p.simulation_lenght,
                                skip=p.time_for_stabilization,
                                choc_time=p.tchoc,
                                kchoc=p.kchoc,
                                tauchoc=p.tauchoc,
                                tdur=p.tdur)


tcomp = time.time() - texp
print("Simulation done")
print("Computational time [s] : ", tcomp)
print("Computational time [min, s] : ", int(tcomp//60), tcomp-(tcomp//60)*60)


# Arrays for save function ---------------------------------------------------------------------------------------------------
dataset_t_mp_soma = Potassic_choc_event[0]
dataset_conc_soma = Potassic_choc_event[1]
dataset_e_soma = Potassic_choc_event[2]
dataset_icl_soma, dataset_ik_soma, dataset_ina_soma = Potassic_choc_event[3], Potassic_choc_event[4], Potassic_choc_event[5]


# path for dataset and name of the file-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
name = f"unclamped_somasim_simlen={p.simulation_lenght:.2g}_dt=({p.dt1},{p.dt2},{p.dt3})_kcc2={p.U_kcc2:.2g}_nkcc1={p.U_nkcc1:.2g}_gclc2={p.soma_gclc2:.2g}_.hdf5"
filepath_h5 = Path.cwd()/"Single_compartment\dataset"/name


# Saving -------------------------------------------------
save_sim_choc_potass(time_mp_soma=dataset_t_mp_soma,
                    soma_concentration=dataset_conc_soma,
                    reversal_pot_soma=dataset_e_soma,
                    soma_currents_cl=dataset_icl_soma,
                    soma_currents_k=dataset_ik_soma,
                    soma_currents_na=dataset_ina_soma,
                    filepath=filepath_h5,
                    dt=(p.dt1, p.dt2, p.dt3),
                    sim_lenght=p.simulation_lenght,
                    rmp=mp,
                    skcc2=p.U_kcc2,
                    snkcc1=p.U_nkcc1,
                    cl_i=my_cell.cli_0,
                    cl_o=my_cell.clo_0,
                    na_i=my_cell.nai_0,
                    na_o=my_cell.nao_0,
                    k_i=my_cell.ki_0,
                    k_o=my_cell.ko_0, 
                    volt_clamp=p.clamp,
                    clamp_v=p.clamp_amp)