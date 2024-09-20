import time
texp = time.time()

from neuron import h
from neuron.units import ms, mV, mM
from matplotlib import rcParams
import numpy as np
import parameters as p
from function import NeuronCell, save_sim, syna
from pathlib import Path
rcParams.update({'font.size': 22}) # Graph parameters


# Loading mod files and for better simulation --------------------------
h.nrn_load_dll(r"mod_files\nrnmech.dll")
h.load_file("stdrun.hoc")


# Temperature in celsius -
h.celsius = p.temperature


# Number of synapses and synapses positions -------------------------------------------------------------
syn_nb = int(p.dend_lenght_with_synapses * p.syn_per_micron) # number of synapses
synapse_pos = np.linspace(0.01, p.dend_lenght_with_synapses/p.dend_lenght, syn_nb) # position of synapses


# Recording vectors positions for concentration (K+, Na+, Cl-) -----------------------
record_pos1 = np.linspace(0.01, (p.position_of_puff*2)/p.dend_lenght, p.number_of_rec)


# Assure that there is a recording vector at the GABA puff position ----------------------------
mini, mini_pos = record_pos1[0], 0
for i,j in enumerate(record_pos1):
    if abs(j-p.position_of_puff/p.dend_lenght) < abs(mini-p.position_of_puff/p.dend_lenght):
        mini = j
        mini_pos = i
record_pos1[mini_pos] = p.position_of_puff/p.dend_lenght


# GABA puff position --------------------------------------------------------------------------------
puff_pos = p.position_of_puff/p.dend_lenght # position of the puff event on the dendrite (first part)


# Creation of the cell ----------------------------------
mp = p.clamp_amp # 1500 ms pip, -90 mV
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
                    ukcc2=p.U_kcc2,
                    unkcc1=p.U_nkcc1)


# GABA puff event and simulation ------------------------------
GABA_puff_event = syna(my_cell,
                        pos=synapse_pos,
                        puff_pos=puff_pos,
                        puff_time=p.time_of_puff * ms,
                        puff_conc=p.concentration_of_puff * mM,
                        tau=p.tau_GABA * ms,
                        dgab=p.Dgaba,
                        rnum=[p.rnum],
                        clamp=p.clamp,
                        clamp_amp=p.clamp_amp * mV,
                        rmp_initial=mp,
                        sim_time=p.simulation_lenght,
                        record_pos=record_pos1,
                        pipette=p.pipett,
                        skip=p.time_for_stabilization)


tcomp = time.time() - texp
print("Simulation done")
print("Computational time [s] : ", tcomp)
print("Computational time [min, s] : ", int(tcomp//60), tcomp-(tcomp//60)*60)


# Arrays for save function ---------------------------------------------------------------------------------------
dataset_t_mp_soma, dataset_mp_dend = GABA_puff_event[0], GABA_puff_event[1]
dataset_conc_dend, dataset_conc_soma = GABA_puff_event[2], GABA_puff_event[3]
dataset_e_soma, dataset_e_dend = GABA_puff_event[4], GABA_puff_event[5]
dataset_icl_dend, dataset_ik_dend, dataset_ina_dend = GABA_puff_event[6], GABA_puff_event[7], GABA_puff_event[8]
dataset_icl_soma, dataset_ik_soma, dataset_ina_soma = GABA_puff_event[9], GABA_puff_event[10], GABA_puff_event[11]
dataset_iother_dend = GABA_puff_event[12]
dataset_g_and_o = GABA_puff_event[13]


# path for dataset and name of the file--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if p.clamp:
    name = f"voltclamped_{p.clamp_amp}mV_syn_nb_{len(synapse_pos)}_sim_lenght_{p.simulation_lenght}_dt_({p.dt1},{p.dt2},{p.dt3})_L_{p.dend_lenght}_kcc2_{p.U_kcc2}_nkcc1_{p.U_nkcc1}_rnum={p.rnum}_puffconc={p.concentration_of_puff}_{my_cell.dend.DCl_iondifus}_gclc2=_{p.soma_gclc2}.hdf5"
    filepath_h5 = Path.cwd()/"dataset"/name
else:
    name = f"unclamped_{mp}mV_syn_nb_{len(synapse_pos)}_sim_lenght_{p.simulation_lenght}_dt_({p.dt1},{p.dt2},{p.dt3})_L_{p.dend_lenght}_kcc2_{p.U_kcc2}_nkcc1_{p.U_nkcc1}_rnum={p.rnum}_puffconc={p.concentration_of_puff}_{my_cell.dend.DCl_iondifus}.hdf5"
    filepath_h5 = Path.cwd()/"dataset"/name


# Saving ----------------------------------------------
save_sim(time_mp_soma=dataset_t_mp_soma,
        mp_dend=dataset_mp_dend,
        dend_concentration=dataset_conc_dend,
        soma_concentration=dataset_conc_soma,
        reversal_pot_soma=dataset_e_soma,
        reversal_pot_dend=dataset_e_dend,
        dend_currents_cl=dataset_icl_dend,
        dend_currents_k=dataset_ik_dend,
        dend_currents_na=dataset_ina_dend,
        soma_currents_cl=dataset_icl_soma,
        soma_currents_k=dataset_ik_soma,
        soma_currents_na=dataset_ina_soma,
        dend_currents_other=dataset_iother_dend,
        dend_g_and_o=dataset_g_and_o,
        filepath=filepath_h5,
        dt=(p.dt1, p.dt2, p.dt3),
        sim_lenght=p.simulation_lenght,
        rnum=p.rnum,
        dend_lenght=(my_cell.dend.L, my_cell.dend2.L),
        record_location=record_pos1,
        rmp=mp,
        nb_dend_seg=(my_cell.dend.nseg, my_cell.dend2.nseg),
        skcc2=p.U_kcc2,
        snkcc1=p.U_nkcc1,
        cl_i=my_cell.cli_0,
        cl_o=my_cell.clo_0,
        na_i=my_cell.nai_0,
        na_o=my_cell.nao_0,
        k_i=my_cell.ki_0,
        k_o=my_cell.ko_0, 
        tau=p.tau_GABA,
        time=p.time_of_puff,
        dgab=p.Dgaba,
        conc=p.concentration_of_puff,
        syn_pos=synapse_pos,
        pos=puff_pos,
        volt_clamp=p.clamp,
        clamp_v=p.clamp_amp)