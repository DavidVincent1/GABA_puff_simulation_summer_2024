from neuron import h
from neuron.units import ms, µm, mV, mM
import matplotlib.pyplot as plt
import math
import h5py
import numpy as np
import parameters_soma as p


# Class 1 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# To create a neuron composed of a soma and a dendrite
class Soma:
    def __init__(self, gid, number_of_soma_segments=1, cli_0=3.763, nai_0=9.814, ki_0=128.347, clo_0=130.5, nao_0=147.5, ko_0=3.5, ukcc2=0.003, unkcc1=2e-5):
        self.number_of_soma_segments = number_of_soma_segments
        self.cli_0 = cli_0
        self.nai_0 = nai_0
        self.ki_0 = ki_0
        self.clo_0 = clo_0
        self.nao_0 = nao_0
        self.ko_0 = ko_0
        self._gid = gid
        self.gabo_0 = 0
        self.kcc2 = ukcc2
        self.nkcc1 = unkcc1
        self._setup_morphology()
        self._setup_biophysics()

    def _setup_morphology(self):
        # Creation of the section (soma)
        self.soma = h.Section(name="soma", cell=self)
    
        # Number of segments in each section
        self.soma.nseg = self.number_of_soma_segments

        # lenght and diameter of each section
        self.soma.L = self.soma.diam = p.soma_diam * µm

    def _setup_biophysics(self):
        # Mecanims insertion in each section
        self.soma.insert("hhrat")
        self.soma.insert("iondifus")
        self.soma.insert("kcc2")
        self.soma.insert("nkcc1")
        self.soma.insert("nakpump")
        self.soma.insert("leak")
        self.soma.insert("clc2")

        # Chloride diffusion coefficient
        self.soma.DCl_iondifus = p.soma_DCl

        # set volume and surface for KCC2 and NKCC1 for soma and dendrite
        self.soma.Vi_kcc2 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_kcc2 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2
        self.soma.Vi_nkcc1 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_nkcc1 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2

        # Na-K pump parameters
        self.soma.imax_nakpump = p.soma_imax
        self.soma.km_k_nakpump = p.soma_kmk
        self.soma.km_na_nakpump = p.soma_kmna

        # Leak paramters
        self.soma.gk_leak = p.soma_gk
        self.soma.gna_leak = p.soma_gna
        self.soma.gnaother_leak = p.soma_gnaother
        self.soma.gcl_leak = p.soma_gcl

        # There is non specific leak channels in the leak mecanism
        # and the HH mecanism, but they are not used.
        self.soma.gfix_leak = 0 
        self.soma.gl_hhrat = 0

        # HH parameters
        self.soma.gnabar_hhrat = p.soma_gnabar   # valeur originale : .12
        self.soma.gkbar_hhrat = p.soma_gkbar # valeur originale : .036

        self.soma.U_kcc2 = self.kcc2
        self.soma.U_nkcc1 = self.nkcc1

        # CLC-2 parameters
        self.soma.gclc2_clc2 = p.soma_gclc2

        # General cell parameters
        self.soma.Ra = p.axial_resistance
        self.soma.cm = p.membrane_capacitance
        self.soma.clamp_iondifus = 0 # It means there is no voltage clamp initially

        # Initial intracellular concentrations
        self.soma.nai = self.nai_0 # Realistics values : 5 to 30 [mM]
        self.soma.ki = self.ki_0   # Realistics values : 60 to 170 [mM]
        self.soma.cli = self.cli_0 # Realistics values : 8 to 30 [mM]
        self.soma.cli0_iondifus = self.cli_0
        self.soma.nai0_iondifus = self.nai_0
        self.soma.ki0_iondifus = self.ki_0
        self.soma.hco3i0_iondifus = 15 # test 15

        # Extracellular concentrations (fix values)
        self.soma.hco3o0_iondifus = 26 # test 25
        self.soma.nao = self.nao_0
        self.soma.ko = self.ko_0 # 3.5
        self.soma.clo = self.clo_0
        self.soma.clo0_iondifus = self.clo_0
        self.soma.nao0_iondifus = self.nao_0
        self.soma.ko0_iondifus = self.ko_0

        # CLC-2 parameters
        self.soma.vhalf_clc2 = p.vhalf
        self.soma.vslope_clc2 = p.vslope
        self.soma.ptau_clc2 = p.ptau

    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)

class Soma_leak_on_clc2:
    def __init__(self, gid, number_of_soma_segments=1, cli_0=3.763, nai_0=9.814, ki_0=128.347, clo_0=130.5,
                nao_0=147.5, ko_0=3.5, ukcc2=0.003, unkcc1=2e-5, gcl_leak=1e-5, gclc2=1e-5):
        self.number_of_soma_segments = number_of_soma_segments
        self.cli_0 = cli_0
        self.nai_0 = nai_0
        self.ki_0 = ki_0
        self.clo_0 = clo_0
        self.nao_0 = nao_0
        self.ko_0 = ko_0
        self._gid = gid
        self.gabo_0 = 0
        self.kcc2 = ukcc2
        self.nkcc1 = unkcc1
        self.gcl_leak = gcl_leak
        self.gclc2 = gclc2
        self._setup_morphology()
        self._setup_biophysics()

    def _setup_morphology(self):
        # Creation of the section (soma)
        self.soma = h.Section(name="soma", cell=self)
    
        # Number of segments in each section
        self.soma.nseg = self.number_of_soma_segments

        # lenght and diameter of each section
        self.soma.L = self.soma.diam = p.soma_diam * µm

    def _setup_biophysics(self):
        # Mecanims insertion in each section
        self.soma.insert("hhrat")
        self.soma.insert("iondifus")
        self.soma.insert("kcc2")
        self.soma.insert("nkcc1")
        self.soma.insert("nakpump")
        self.soma.insert("leak")
        self.soma.insert("clc2")

        # Chloride diffusion coefficient
        self.soma.DCl_iondifus = p.soma_DCl

        # set volume and surface for KCC2 and NKCC1 for soma and dendrite
        self.soma.Vi_kcc2 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_kcc2 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2
        self.soma.Vi_nkcc1 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_nkcc1 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2

        # Na-K pump parameters
        self.soma.imax_nakpump = p.soma_imax
        self.soma.km_k_nakpump = p.soma_kmk
        self.soma.km_na_nakpump = p.soma_kmna

        # Leak paramters
        self.soma.gk_leak = p.soma_gk
        self.soma.gna_leak = p.soma_gna
        self.soma.gnaother_leak = p.soma_gnaother
        self.soma.gcl_leak = self.gcl_leak

        # There is non specific leak channels in the leak mecanism
        # and the HH mecanism, but they are not used.
        self.soma.gfix_leak = 0 
        self.soma.gl_hhrat = 0

        # HH parameters
        self.soma.gnabar_hhrat = p.soma_gnabar   # valeur originale : .12
        self.soma.gkbar_hhrat = p.soma_gkbar # valeur originale : .036

        self.soma.U_kcc2 = self.kcc2
        self.soma.U_nkcc1 = self.nkcc1

        # CLC-2 parameters
        self.soma.gclc2_clc2 = self.gclc2

        # General cell parameters
        self.soma.Ra = p.axial_resistance
        self.soma.cm = p.membrane_capacitance
        self.soma.clamp_iondifus = 0 # It means there is no voltage clamp initially

        # Initial intracellular concentrations
        self.soma.nai = self.nai_0 # Realistics values : 5 to 30 [mM]
        self.soma.ki = self.ki_0   # Realistics values : 60 to 170 [mM]
        self.soma.cli = self.cli_0 # Realistics values : 8 to 30 [mM]
        self.soma.cli0_iondifus = self.cli_0
        self.soma.nai0_iondifus = self.nai_0
        self.soma.ki0_iondifus = self.ki_0
        self.soma.hco3i0_iondifus = 15 # test 15

        # Extracellular concentrations (fix values)
        self.soma.hco3o0_iondifus = 26 # test 25
        self.soma.nao = self.nao_0
        self.soma.ko = self.ko_0 # 3.5
        self.soma.clo = self.clo_0
        self.soma.clo0_iondifus = self.clo_0
        self.soma.nao0_iondifus = self.nao_0
        self.soma.ko0_iondifus = self.ko_0

        # CLC-2 parameters
        self.soma.vhalf_clc2 = p.vhalf
        self.soma.vslope_clc2 = p.vslope
        self.soma.ptau_clc2 = p.ptau

    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)


# Class 1 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 2 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# Main function to do the GABA puff simulation
def choc_potass(cell, clamp=False, clamp_amp=-70, rmp_initial=-72.38, sim_time=10000, skip=0, choc_time=10000, kchoc=13.5, tauchoc=10000, tdur=10000):
    # Initialization at 0 of the 'messenger' concentration ----------
    # This 'messenger' makes the link between the puff (puff.mod) and
    # the extracellular GABA concentration (iondiffus.mod)
    for seg in cell.soma:
        seg.messi = 0
    
    # Setting the potassic choc ------------------------
    choc = h.CHOCpot(cell.soma(0.5))
    choc.tchoc = choc_time * ms
    choc.kchoc = kchoc * mM
    choc.tauchoc = tauchoc * ms
    choc.ko0 = p.ko
    choc.tdur = tdur * ms

    # Recording vectors for concentrations --------------
    soma_cli = h.Vector().record(cell.soma(0.5)._ref_cli)
    soma_ko = h.Vector().record(cell.soma(0.5)._ref_ko)
    soma_nai = h.Vector().record(cell.soma(0.5)._ref_nai)
    soma_ki = h.Vector().record(cell.soma(0.5)._ref_ki)

    # Chloride currents, soma ---------------------------------------
    soma_icl = h.Vector().record(cell.soma(0.5)._ref_icl)
    soma_icl_kcc2 = h.Vector().record(cell.soma(0.5)._ref_icl_kcc2)
    soma_icl_nkcc1 = h.Vector().record(cell.soma(0.5)._ref_icl_nkcc1)
    soma_icl_leak = h.Vector().record(cell.soma(0.5)._ref_icl_leak)
    soma_icl_clc2 = h.Vector().record(cell.soma(0.5)._ref_icl_clc2)

    # Potassium currents, soma ------------------------------------
    soma_ik = h.Vector().record(cell.soma(0.5)._ref_ik)
    soma_ik_kcc2 = h.Vector().record(cell.soma(0.5)._ref_ik_kcc2)
    soma_ik_nkcc1 = h.Vector().record(cell.soma(0.5)._ref_ik_nkcc1)
    soma_ik_leak = h.Vector().record(cell.soma(0.5)._ref_ik_leak)
    soma_ik_nak = h.Vector().record(cell.soma(0.5)._ref_ik_nakpump)
    soma_ik_hh = h.Vector().record(cell.soma(0.5)._ref_ik_hhrat)

    # Sodium currents, soma -----------------------------------------
    soma_ina = h.Vector().record(cell.soma(0.5)._ref_ina)
    soma_ina_nkcc1 = h.Vector().record(cell.soma(0.5)._ref_ina_nkcc1)
    soma_ina_leak = h.Vector().record(cell.soma(0.5)._ref_ina_leak)
    soma_ina_nak = h.Vector().record(cell.soma(0.5)._ref_ina_nakpump)
    soma_ina_hh = h.Vector().record(cell.soma(0.5)._ref_ina_hhrat)

    # Recording vectors for reversal potential, soma ----
    soma_ecl = h.Vector().record(cell.soma(0.5)._ref_ecl)
    soma_ena = h.Vector().record(cell.soma(0.5)._ref_ena)
    soma_ek = h.Vector().record(cell.soma(0.5)._ref_ek)

    # Recording vector for membrane potential, soma -
    soma_v = h.Vector().record(cell.soma(0.5)._ref_v)

    # Recording vector for the time and initialization of the membrane potential
    t = h.Vector().record(h._ref_t)
    if clamp:
        h.finitialize(clamp_amp*mV)
    else:
        h.finitialize(rmp_initial)


    # Simulation -------------------------------------------
    # Separated in 3 temporal window
    # First at dt1 for the initial stabilization of the cell
    # Second at dt2 for the time juste before the puff
    # Third at dt3 for the rest of the simulation 
    h.tstop = sim_time
    if sim_time > skip:
        h.dt = p.dt1*ms
        h.continuerun(skip*ms)
    h.dt = p.dt2*ms
    h.continuerun((choc_time-100)*ms)
    h.dt = p.dt3*ms
    h.continuerun(sim_time)


    # Surface of the soma ----------------------------
    # Correspond to (segment lenght) * (circonference)
    soma_surface_area = cell.soma(0.5).area()
    soma_surface_area *= (1e-8) # in cm2

    # Current density to total current conversion in soma
    soma_icl *= soma_surface_area * (1e9)
    soma_icl_kcc2 *= soma_surface_area * (1e9)
    soma_icl_nkcc1 *= soma_surface_area * (1e9)
    soma_icl_leak *= soma_surface_area * (1e9)
    soma_icl_clc2 *= soma_surface_area * (1e9)

    soma_ik *= soma_surface_area * (1e9)
    soma_ik_kcc2 *= soma_surface_area * (1e9)
    soma_ik_nkcc1 *= soma_surface_area * (1e9)
    soma_ik_leak *= soma_surface_area * (1e9)
    soma_ik_nak *= soma_surface_area * (1e9)
    soma_ik_hh *= soma_surface_area * (1e9)

    soma_ina *= soma_surface_area * (1e9)
    soma_ina_nkcc1 *= soma_surface_area * (1e9)
    soma_ina_leak *= soma_surface_area * (1e9)
    soma_ina_nak *= soma_surface_area * (1e9)
    soma_ina_hh *= soma_surface_area * (1e9)


    # Results tuple that will be return by the function ------------------------------------------
    result_mp = (t, soma_v)
    result_con_soma = (soma_cli, soma_ki, soma_nai, soma_ko)
    result_e_soma = (soma_ecl, soma_ek, soma_ena)
    result_icl_soma = (soma_icl, soma_icl_kcc2, soma_icl_nkcc1, soma_icl_leak, soma_icl_clc2)
    result_ik_soma = (soma_ik, soma_ik_kcc2, soma_ik_nkcc1, soma_ik_leak, soma_ik_nak, soma_ik_hh)
    result_ina_soma = (soma_ina, soma_ina_nkcc1, soma_ina_leak, soma_ina_nak, soma_ina_hh)

    return result_mp, result_con_soma, result_e_soma, result_icl_soma, result_ik_soma, result_ina_soma


# Function 2 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 3 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# Fonction pour sauvegarder les résultats d'une simulation
def save_sim_choc_potass(time_mp_soma, soma_concentration,
                        reversal_pot_soma, soma_currents_cl,
                        soma_currents_k, soma_currents_na, filepath,
                        dt, sim_lenght, rmp, skcc2, snkcc1, cl_i, cl_o, na_i,
                        na_o, k_i, k_o, volt_clamp=False, clamp_v=False):
    h5_file = h5py.File(filepath, "w")

    dataset1 = h5_file.create_dataset("mp_soma", data=time_mp_soma)
    dataset1.attrs["Shape"] = f"Time array mp_soma[0]\nMP at 0.5 in soma mp_soma[1]"

    dataset2 = h5_file.create_dataset("conc_soma", data=soma_concentration)
    dataset2.attrs["temporal_resolution"] = dt
    dataset2.attrs["simulation_length"] = sim_lenght
    dataset2.attrs["voltage_clamped"] = volt_clamp
    dataset2.attrs["voltage_clamp_mV"] = clamp_v
    dataset2.attrs["initial_membrane_potential_mV"] = rmp
    dataset2.attrs["kcc2_strength"] = skcc2
    dataset2.attrs["nkcc1_strength"] = snkcc1
    dataset2.attrs["cli_0"] = cl_i
    dataset2.attrs["clo_0"] = cl_o
    dataset2.attrs["nai_0"] = na_i
    dataset2.attrs["nao_0"] = na_o
    dataset2.attrs["ki_0"] = k_i
    dataset2.attrs["ko_0"] = k_o
    dataset2.attrs["Shape"] = f"""Chloride concentration at 0.5 in soma conc_soma[0]
    Potassium concentration at 0.5 in soma conc_soma[1]
    Sodium concentration at 0.5 in soma conc_soma[2]
    External GABA concentration at 0.5 in soma conc_soma[3]
    Potassium extracellular concentration at 0.5 in soma conc_soma[4]"""

    dataset3 = h5_file.create_dataset("e_soma", data=reversal_pot_soma)
    dataset3.attrs["Shape"] = f"""Ecl at 0.5 in soma e_soma[0]
    Ek at 0.5 in soma e_soma[1]
    Ena at 0.5 in soma e_soma[2]""" 

    dataset4 = h5_file.create_dataset("soma_current_cl", data=soma_currents_cl)
    dataset4.attrs["Shape"] = f"""Icl at 0.5 in soma soma_current_cl[0]
    Icl_kcc2 at 0.5 in soma soma_current_cl[1]
    Icl_nkcc1 at 0.5 in soma soma_current_cl[2]
    Icl_leak at 0.5 in soma soma_current_cl[3]
    Icl_clc2 at 0.5 in soma soma_current_cl[4]"""   

    dataset5 = h5_file.create_dataset("soma_current_k", data=soma_currents_k)
    dataset5.attrs["Shape"] = f"""Ik at 0.5 in soma soma_current_k[0]
    Ik_kcc2 at 0.5 in soma soma_current_k[1]
    Ik_nkcc1 at 0.5 in soma soma_current_k[2]
    Ik_leak at 0.5 in soma soma_current_k[3]
    Ik_hh at 0.5 in soma soma_current_k[4]"""   

    dataset6 = h5_file.create_dataset("soma_current_na", data=soma_currents_na)
    dataset6.attrs["Shape"] = f"""Ina at 0.5 in soma soma_current_na[0]
    Ina_nkcc1 at 0.5 in soma soma_current_na[1]
    Ina_leak at 0.5 in soma soma_current_na[2]
    Ina_hh at 0.5 in soma soma_current_na[3]"""    


    h5_file.close()


# Function 3 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 4 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# To print information about a dataset
def show_info_sim(time_mp_soma, mp_dend,
                        dend_concentration, soma_concentration,
                        reversal_pot_soma, reversal_pot_dend,
                        dend_currents_cl, dend_currents_k, dend_currents_na,
                        soma_currents_cl, soma_currents_k, soma_currents_na,
                        dend_currents_other, dend_g):
    """Print to the terminal information about a dataset.

    Args:
        time_mp_soma (array): Time and membrane potential in soma array (from a h5py dataset)
        mp_dend (array): Membrane potential in dendrite array (from a h5py dataset)
        dend_concentration (array): Concentrations in dendrite array (from a h5py dataset)
        soma_concentration (array): Concentrations in soma array (from a h5py dataset)
        reversal_pot_soma (array): Reversal potentials in soma array (from a h5py dataset)
        reversal_pot_dend (array): Reversal potentials in dendrite array (from a h5py dataset)
        dend_currents_cl (array): Chloride currents in dendrite array (from a h5py dataset)
        dend_currents_k (array): Potassium currents in dendrite array (from a h5py dataset)
        dend_currents_na (array): Sodium currents in dendrite array (from a h5py dataset)
        soma_currents_cl (array): Chloride currents in soma array (from a h5py dataset)
        soma_currents_k (array): Potassium currents in soma array (from a h5py dataset)
        soma_currents_na (array): Sodium currents in soma array (from a h5py dataset)
        dend_currents_other (array): Synapse other currents in dendrite array (from a h5py dataset)
        dend_g (array): Conductances at synapses in dendrite array (from a h5py dataset)
    """
    print(r'Temporal resolution (ms) : ', dend_concentration.attrs["temporal_resolution"])
    print(r'Simulation lenght (ms)   : ', dend_concentration.attrs["simulation_length"])
    print('Recording location (fraction of the dendrite lenght) : \n', dend_concentration.attrs["recording_locations"])
    print(r'Number of dendrite segments in first and second parts (-) : ', dend_concentration.attrs["number_dendrite_segments"])
    print(r'Voltage clamp (True or False) : ', dend_concentration.attrs["voltage_clamped"])
    print(r'Voltage clamp value (mV)      : ', dend_concentration.attrs["voltage_clamp_mV"])
    print("")

    # Information about the neuron
    print('Synapses positions (fraction of dendrite lenght) : \n', dend_concentration.attrs["Synapses positions"])
    print(r'Number of GABA receptors on synapse (-) : ', dend_concentration.attrs["number_of_receptors"])
    print(r'Dendrite lenght (um)   : ', dend_concentration.attrs["dendrite_full_length_um"])
    print(r'KCC2 strenght (mM/ms)  : ', dend_concentration.attrs["kcc2_strength"])
    print(r'nKCC1 strenght (mM/ms) : ', dend_concentration.attrs["nkcc1_strength"])
    print(r'Initial [Cl^-]_i (mM)  : ', dend_concentration.attrs["cli_0"])
    print(r'[Cl^-]_o (mM)          : ', dend_concentration.attrs["clo_0"])
    print(r'Initial [K^+]_i (mM)   : ', dend_concentration.attrs["ki_0"])
    print(r'[K^+]_o (mM)           : ', dend_concentration.attrs["ko_0"])
    print(r'Initial [Na^+]_i (mM)  : ', dend_concentration.attrs["nai_0"])
    print(r'[Na^+]_o (mM)          : ', dend_concentration.attrs["nao_0"])
    print(r'Initial MP (mV)        : ', dend_concentration.attrs["initial_membrane_potential_mV"])
    print("")

    # Information on the GABA puff
    print(r'Puff position (fraction of dendrite lenght)        : ', dend_concentration.attrs["puff_position"])
    print(r'Puff concentration (mM)                            : ', dend_concentration.attrs["puff_concentration"])
    print(r'Time at wich the puff event occurs (ms)            : ', dend_concentration.attrs["puff_time"])
    print(r'Longitunal diffusion coefficient for GABA (um2/ms) : ', dend_concentration.attrs["puff_longitudinal_diffusion_coefficient"])
    print(r'Exchange with bath time constant (ms)              : ', dend_concentration.attrs["puff_taubath"])
    print("")

    # Shape of data to use them
    print(time_mp_soma.attrs["Shape"])
    print('')
    print(mp_dend.attrs["Shape"])
    print('')
    print(dend_concentration.attrs["Shape"])
    print('')
    print(soma_concentration.attrs["Shape"])
    print('')
    print(reversal_pot_soma.attrs["Shape"])
    print('')
    print(reversal_pot_dend.attrs["Shape"])
    print('')
    print(dend_currents_cl.attrs["Shape"])
    print('')
    print(dend_currents_k.attrs["Shape"])
    print('')
    print(dend_currents_na.attrs["Shape"])
    print('')
    print(soma_currents_cl.attrs["Shape"])
    print('')
    print(soma_currents_k.attrs["Shape"])
    print('')
    print(soma_currents_na.attrs["Shape"])
    print('')
    print(dend_currents_other.attrs["Shape"])
    print('')
    print(dend_g.attrs["Shape"])
    print('')

# Function 4 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 5 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# Simulation for 2D matix
def syna_kcc2_nkcc1_soma(cell, choc_time, kchoc, tauchoc, clamp=False, clamp_amp=-70,
                    rmp_initial=-72.38, sim_time=10000, skip=0, pipette=(1500*ms, 8*mM, 140*mM, 12*mM)):
    """_summary_

    Args:
        first (bool or int): If True, the function calculate the GABA puff index. Else, it is the puff index.
        cell (NeuronCell class object): The cell of the simulation.
        puff_pos (float): Position of the GABA puff on the dendrite (first part). Number between 0 and 1.
        puff_time (float): Time of the GABA puff in the simulation [ms].
        puff_conc (float): Concentration of the GABA puff [mM]
        pos (array): Array of numbers between 0 and 1 corresponding to synapse positions.
        tau (float): Exchange with bath constant of the GABA [ms]
        dgab (float): Diffusion coefficient of the GABA [µm2/ms]
        rnum (float or list): Number of receptors per synapse.
                              If the input is a list, it corresponds to the number of receptors for each synapses individually.
                              If it's a float, it corresponds to the number of receptors for each synapses.
        clamp (bool, optional): True if there is a voltage clamp on the soma. False if not. Defaults to False.
        clamp_amp (int, optional): Amplitude of the voltage clamp if there is one. Defaults to -70. [mV]
        rmp_initial (float, optional): Initial membrane potential. Defaults to -72.38. [mV]
        sim_time (int, optional): Lenght of the simulation. Defaults to 10000. [ms]
        record_pos (array, optional): Recording vectors positions. Defaults to np.linspace(0.01,0.99,100).
        skip (float, optional): Time for the initial stabilization of the cell. Defaults to 0.
        pipette (tuple, optional): Caracteristic of the patch clamp pipette if there is a voltage clamp.
                                   Defaults to (1500*ms, 8*mM, 140*mM, 12*mM).

    Returns:
        _type_: _description_
    """

    # If the soma is voltage clamped -----------------------------------------------------------
    if clamp:
        # Simulate an electrode
        # See https://github.com/neuronsimulator/nrn/blob/master/src/nrnoc/svclmp.mod
        cell.soma.clamp_iondifus = 1       # Activates the exchange with the patch clamp pipette
        vclamp = h.SEClamp(cell.soma(0.5)) # The clamp is added to the soma here
        vclamp.amp1 = clamp_amp*mV         # Amplitude of the voltage clamp [mV]
        vclamp.dur1 = sim_time*ms          # Duration of the clamp [ms]
        vclamp.rs = 4                      # Resistance of the electrode [megaohm]
    
    # Parameters of the pipette if there's a clamp
    h.tau_iondifus = pipette[0]    # Exchange with pipette constant [ms]
    h.clipip_iondifus = pipette[1] # Chloride concentration in piette solution [mM]
    h.kipip_iondifus = pipette[2]  # Potassium concentration in piette solution [mM]
    h.naipip_iondifus = pipette[3] # Sodium concentration in piette solution [mM]

    # Setting the potassic choc ------------------------
    choc = h.CHOCpot(cell.soma(0.5))
    choc.tchoc = choc_time * ms
    choc.kchoc = kchoc * mM
    choc.tauchoc = tauchoc * ms
    choc.ko0 = p.ko

    # Initialization at 0 of the 'messenger' concentration ----------
    # This 'messenger' makes the link between the puff (puff.mod) and
    # the extracellular GABA concentration (iondiffus.mod)
    for seg in cell.soma:
        seg.messi = 0

    # Recording vectors for concentrations --------------
    soma_cli = h.Vector().record(cell.soma(0.5)._ref_cli)
    soma_nai = h.Vector().record(cell.soma(0.5)._ref_nai)
    soma_ki = h.Vector().record(cell.soma(0.5)._ref_ki)

    # Chloride currents, soma ---------------------------------------
    soma_icl = h.Vector().record(cell.soma(0.5)._ref_icl)
    soma_icl_kcc2 = h.Vector().record(cell.soma(0.5)._ref_icl_kcc2)
    soma_icl_nkcc1 = h.Vector().record(cell.soma(0.5)._ref_icl_nkcc1)
    soma_icl_leak = h.Vector().record(cell.soma(0.5)._ref_icl_leak)
    soma_icl_clc2 = h.Vector().record(cell.soma(0.5)._ref_icl_clc2)

    # Recording vectors for reversal potential, soma ----
    soma_ecl = h.Vector().record(cell.soma(0.5)._ref_ecl)

    # Recording vector for membrane potential, soma -
    soma_v = h.Vector().record(cell.soma(0.5)._ref_v)

    # Recording vector for the time and initialization of the membrane potential
    t = h.Vector().record(h._ref_t)
    if clamp:
        h.finitialize(clamp_amp*mV)
    else:
        h.finitialize(rmp_initial)

    # Simulation -------------------------------------------
    # Separated in 3 temporal window
    # First at dt1 for the initial stabilization of the cell
    # Second at dt2 for the time juste before the puff
    # Third at dt3 for the rest of the simulation 
    h.tstop = sim_time
    if sim_time > skip:
        h.dt = p.dt1*ms
        h.continuerun(skip*ms)
    h.dt = p.dt2*ms
    h.continuerun((choc_time-100)*ms)
    h.dt = p.dt3*ms
    h.continuerun(sim_time)

    # Surface of the soma ----------------------------
    # Correspond to (segment lenght) * (circonference)
    soma_surface_area = cell.soma(0.5).area()
    soma_surface_area *= (1e-8) # in cm2

    # Current density to total current conversion in soma
    soma_icl *= soma_surface_area * (1e9)
    soma_icl_kcc2 *= soma_surface_area * (1e9)
    soma_icl_nkcc1 *= soma_surface_area * (1e9)
    soma_icl_leak *= soma_surface_area * (1e9)
    soma_icl_clc2 *= soma_surface_area * (1e9)

    # Index of the time just before the GABA puff event
    potassic_ind = int(skip/p.dt1 + (choc_time-skip)/p.dt2)
    stable_conc = int(potassic_ind - 2)

    # Maximum delta chloride -----------------------------------
    soma_cli_numpy = np.array(soma_cli)

    max_chloride = max(soma_cli_numpy[potassic_ind:])
    min_chloride = min(soma_cli_numpy[potassic_ind:])

    delta_chloride_max = max_chloride - soma_cli[stable_conc]
    delta_chloride_min = min_chloride - soma_cli[stable_conc]
    delta_chloride = 0
    if abs(delta_chloride_max) > abs(delta_chloride_min):
        delta_chloride = delta_chloride_max
    else:
        delta_chloride = delta_chloride_min

    return (delta_chloride, soma_cli[stable_conc], soma_ki[stable_conc], soma_nai[stable_conc], soma_v[stable_conc], 
            soma_ecl[stable_conc], soma_cli, soma_ki, soma_nai, soma_v, t, soma_icl[stable_conc], soma_icl_kcc2[stable_conc],
            soma_icl_nkcc1[stable_conc], soma_icl_clc2[stable_conc], soma_icl_leak[stable_conc])

# Function 5 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------