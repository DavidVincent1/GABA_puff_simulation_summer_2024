from neuron import h
from neuron.units import ms, µm, mV, mM
import matplotlib.pyplot as plt
import math
import h5py
import numpy as np
import parameters as p


# Class 1 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# To create a neuron composed of a soma and a dendrite
class NeuronCell:
    def __init__(self, gid, number_of_dendrite_segments=1, number_of_dendrite_segments2=1, number_of_soma_segments=1, dendrite_length_um=100, dendrite_length_um2=100, cli_0=3.763, nai_0=9.814, ki_0=128.347, clo_0=130.5, nao_0=147.5, ko_0=3.5, ukcc2=0.003, unkcc1=2e-5):
        self.number_of_dendrite_segments = number_of_dendrite_segments
        self.number_of_dendrite_segments2 = number_of_dendrite_segments2
        self.number_of_soma_segments = number_of_soma_segments
        self.dendrite_length_um = dendrite_length_um
        self.dendrite_length_um2 = dendrite_length_um2
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
        # Creation of the sections (soma and dendrite)
        self.soma = h.Section(name="soma", cell=self)
        self.dend = h.Section(name="dend", cell=self)
        self.dend2 = h.Section(name="dend2", cell=self)
    
        # Number of segments in each section
        self.soma.nseg = self.number_of_soma_segments
        self.dend.nseg = self.number_of_dendrite_segments
        self.dend2.nseg = self.number_of_dendrite_segments2
    
        # Connection of the soma and the dendrite
        self.dend.connect(self.soma)
        self.dend2.connect(self.dend)
        self.all = self.soma.wholetree()

        # lenght and diameter of each section
        self.soma.L = self.soma.diam = p.soma_diam * µm
        self.dend.L = self.dendrite_length_um * µm
        self.dend2.L = self.dendrite_length_um2 * µm
        self.dend.diam = p.dend_diam * µm#1 * µm
        self.dend2.diam = p.dend2_diam * µm#1 * µm

    def _setup_biophysics(self):
        # Mecanims insertion in each section
        self.soma.insert("hhrat")
        self.soma.insert("iondifus")
        self.soma.insert("kcc2")
        self.soma.insert("nkcc1")
        self.soma.insert("nakpump")
        self.soma.insert("leak")
        self.soma.insert("clc2")

        self.dend.insert("hhrat")
        self.dend.insert("iondifus")
        self.dend.insert("kcc2")
        self.dend.insert("nkcc1")
        self.dend.insert("nakpump")
        self.dend.insert("leak")
        self.dend.insert("clc2")

        self.dend2.insert("hhrat")
        self.dend2.insert("iondifus")
        self.dend2.insert("kcc2")
        self.dend2.insert("nkcc1")
        self.dend2.insert("nakpump")
        self.dend2.insert("leak")
        self.dend2.insert("clc2")

        # Chloride diffusion coefficient
        self.soma.DCl_iondifus = p.soma_DCl
        self.dend.DCl_iondifus = p.dend_DCl
        self.dend2.DCl_iondifus = p.dend2_DCl

        # set volume and surface for KCC2 and NKCC1 for soma and dendrite
        self.soma.Vi_kcc2 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_kcc2 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2
        self.soma.Vi_nkcc1 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_nkcc1 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2

        self.dend.Vi_kcc2 = self.dend.L*math.pi*(self.dend.diam/2)**2
        self.dend.S_kcc2 = 2*math.pi*self.dend.diam/2*self.dend.L+2*math.pi*(self.dend.diam/2)**2
        self.dend.Vi_nkcc1 = self.dend.L*math.pi*(self.dend.diam/2)**2
        self.dend.S_nkcc1 = 2*math.pi*self.dend.diam/2*self.dend.L+2*math.pi*(self.dend.diam/2)**2

        self.dend2.Vi_kcc2 = self.dend2.L*math.pi*(self.dend2.diam/2)**2
        self.dend2.S_kcc2 = 2*math.pi*self.dend2.diam/2*self.dend2.L+2*math.pi*(self.dend2.diam/2)**2
        self.dend2.Vi_nkcc1 = self.dend2.L*math.pi*(self.dend2.diam/2)**2
        self.dend2.S_nkcc1 = 2*math.pi*self.dend2.diam/2*self.dend2.L+2*math.pi*(self.dend2.diam/2)**2

        # Na-K pump parameters
        self.soma.imax_nakpump = p.soma_imax
        self.dend.imax_nakpump = p.dend_imax
        self.dend2.imax_nakpump = p.dend2_imax

        self.soma.km_k_nakpump = p.soma_kmk
        self.dend.km_k_nakpump = p.dend_kmk
        self.dend2.km_k_nakpump = p.dend2_kmk

        self.soma.km_na_nakpump = p.soma_kmna
        self.dend.km_na_nakpump = p.dend_kmna
        self.dend2.km_na_nakpump = p.dend2_kmna

        # Leak paramters
        self.soma.gk_leak = p.soma_gk
        self.dend.gk_leak = p.dend_gk
        self.dend2.gk_leak = p.dend2_gk

        self.soma.gna_leak = p.soma_gna
        self.dend.gna_leak = p.dend_gna
        self.dend2.gna_leak = p.dend2_gna

        self.soma.gnaother_leak = p.soma_gnaother
        self.dend.gnaother_leak = p.dend_gnaother
        self.dend2.gnaother_leak = p.dend2_gnaother

        self.soma.gcl_leak = p.soma_gcl
        self.dend.gcl_leak = p.dend_gcl
        self.dend2.gcl_leak = p.dend2_gcl

        # There is non specific leak channels in the leak mecanism
        # and the HH mecanism, but they are not used.
        self.soma.gfix_leak = 0 
        self.dend.gfix_leak = 0
        self.dend2.gfix_leak = 0
        self.soma.gl_hhrat = 0
        self.dend.gl_hhrat = 0
        self.dend2.gl_hhrat = 0

        # HH parameters
        self.soma.gnabar_hhrat = p.soma_gnabar   # valeur originale : .12
        self.dend.gnabar_hhrat = p.dend_gnabar
        self.dend2.gnabar_hhrat = p.dend2_gnabar

        self.soma.gkbar_hhrat = p.soma_gkbar # valeur originale : .036
        self.dend.gkbar_hhrat = p.dend_gkbar
        self.dend2.gkbar_hhrat = p.dend2_gkbar

        self.soma.U_kcc2 = self.kcc2
        self.dend.U_kcc2 = self.kcc2 * (self.soma.Vi_kcc2 * self.dend.S_kcc2)/(self.soma.S_kcc2 * self.dend.Vi_kcc2)
        self.dend2.U_kcc2 = self.kcc2 * (self.soma.Vi_kcc2 * self.dend2.S_kcc2)/(self.soma.S_kcc2 * self.dend2.Vi_kcc2)

        self.soma.U_nkcc1 = self.nkcc1
        self.dend.U_nkcc1 = self.nkcc1 * (self.soma.Vi_nkcc1 * self.dend.S_nkcc1)/(self.soma.S_nkcc1 * self.dend.Vi_nkcc1)
        self.dend2.U_nkcc1 = self.nkcc1 * (self.soma.Vi_nkcc1 * self.dend2.S_nkcc1)/(self.soma.S_nkcc1 * self.dend2.Vi_nkcc1)

        # CLC-2 parameters
        self.soma.gclc2_clc2 = p.soma_gclc2
        self.dend.gclc2_clc2 = p.dend_gclc2
        self.dend2.gclc2_clc2 = p.dend2_gclc2



        for sec1 in self.all:
            # General cell parameters
            sec1.Ra = p.axial_resistance
            sec1.cm = p.membrane_capacitance
            sec1.clamp_iondifus = 0 # It means there is no voltage clamp initially

            # Initial intracellular concentrations
            sec1.nai = self.nai_0 # Realistics values : 5 to 30 [mM]
            sec1.ki = self.ki_0   # Realistics values : 60 to 170 [mM]
            sec1.cli = self.cli_0 # Realistics values : 8 to 30 [mM]
            sec1.cli0_iondifus = self.cli_0
            sec1.nai0_iondifus = self.nai_0
            sec1.ki0_iondifus = self.ki_0
            sec1.hco3i0_iondifus = 15 # test 15

            # Extracellular concentrations (fix values)
            sec1.hco3o0_iondifus = 26 # test 25
            sec1.nao = self.nao_0
            sec1.ko = self.ko_0 # 3.5
            sec1.clo = self.clo_0
            sec1.clo0_iondifus = self.clo_0
            sec1.nao0_iondifus = self.nao_0
            sec1.ko0_iondifus = self.ko_0

            # CLC-2 parameters
            sec1.vhalf_clc2 = p.vhalf
            sec1.vslope_clc2 = p.vslope
            sec1.ptau_clc2 = p.ptau

    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)

class NeuronCellOneFork:
    def __init__(self, gid, number_of_dendrite_segments=1, number_of_dendrite_segments2=1,
                number_of_fork_segments=1, number_of_soma_segments=1, dendrite_length_um=100,
                dendrite_length_um2=100, fork_length_um=100, cli_0=3.763, nai_0=9.814,
                ki_0=128.347, clo_0=130.5, nao_0=147.5, ko_0=3.5, ukcc2=0.003, unkcc1=2e-5,
                fork_position=100/2):
        self.number_of_dendrite_segments = number_of_dendrite_segments
        self.number_of_dendrite_segments2 = number_of_dendrite_segments2
        self.number_of_soma_segments = number_of_soma_segments
        self.number_of_fork_segments = number_of_fork_segments

        self.dendrite_length_um = dendrite_length_um
        self.dendrite_length_um2 = dendrite_length_um2
        self.fork_length_um = fork_length_um
        self.fork_position = fork_position

        self.cli_0 = cli_0
        self.nai_0 = nai_0
        self.ki_0 = ki_0
        self.clo_0 = clo_0
        self.nao_0 = nao_0
        self.ko_0 = ko_0
        self.gabo_0 = 0

        self._gid = gid
        self.kcc2 = ukcc2
        self.nkcc1 = unkcc1
        self._setup_morphology()
        self._setup_biophysics()

    def _setup_morphology(self):
        # Creation of the sections (soma and dendrite)
        self.soma = h.Section(name="soma", cell=self)
        self.dend = h.Section(name="dend", cell=self)
        self.dend2 = h.Section(name="dend2", cell=self)
        self.fork = h.Section(name="fork", cell=self)
    
        # Number of segments in each section
        self.soma.nseg = self.number_of_soma_segments
        self.dend.nseg = self.number_of_dendrite_segments
        self.dend2.nseg = self.number_of_dendrite_segments2
        self.fork.nseg = self.number_of_fork_segments
    
        # Connection of the soma and the dendrite
        self.dend.connect(self.soma)
        self.dend2.connect(self.dend)
        self.fork.connect(self.dend(self.fork_position), 0)
        self.all = self.soma.wholetree()

        # lenght and diameter of each section
        self.soma.L = self.soma.diam = p.soma_diam * µm
        self.dend.L = self.dendrite_length_um * µm
        self.dend2.L = self.dendrite_length_um2 * µm
        self.dend.diam = p.dend_diam * µm#1 * µm
        self.dend2.diam = p.dend2_diam * µm#1 * µm
        self.fork.L = self.fork_length_um * µm

    def _setup_biophysics(self):
        # Mecanims insertion in each section
        self.soma.insert("hhrat")
        self.soma.insert("iondifus")
        self.soma.insert("kcc2")
        self.soma.insert("nkcc1")
        self.soma.insert("nakpump")
        self.soma.insert("leak")
        self.soma.insert("clc2")

        self.dend.insert("hhrat")
        self.dend.insert("iondifus")
        self.dend.insert("kcc2")
        self.dend.insert("nkcc1")
        self.dend.insert("nakpump")
        self.dend.insert("leak")
        self.dend.insert("clc2")

        self.dend2.insert("hhrat")
        self.dend2.insert("iondifus")
        self.dend2.insert("kcc2")
        self.dend2.insert("nkcc1")
        self.dend2.insert("nakpump")
        self.dend2.insert("leak")
        self.dend2.insert("clc2")

        self.fork.insert("hhrat")
        self.fork.insert("iondifus")
        self.fork.insert("kcc2")
        self.fork.insert("nkcc1")
        self.fork.insert("nakpump")
        self.fork.insert("leak")
        self.fork.insert("clc2")

        # Chloride diffusion coefficient
        self.soma.DCl_iondifus = p.soma_DCl
        self.dend.DCl_iondifus = p.dend_DCl
        self.dend2.DCl_iondifus = p.dend2_DCl
        self.fork.DCl_iondifus = p.fork_DCl

        # set volume and surface for KCC2 and NKCC1 for soma and dendrite
        self.soma.Vi_kcc2 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_kcc2 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2
        self.soma.Vi_nkcc1 = self.soma.L*math.pi*(self.soma.diam/2)**2
        self.soma.S_nkcc1 = 2*math.pi*self.soma.diam/2*self.soma.L+2*math.pi*(self.soma.diam/2)**2

        self.dend.Vi_kcc2 = self.dend.L*math.pi*(self.dend.diam/2)**2
        self.dend.S_kcc2 = 2*math.pi*self.dend.diam/2*self.dend.L+2*math.pi*(self.dend.diam/2)**2
        self.dend.Vi_nkcc1 = self.dend.L*math.pi*(self.dend.diam/2)**2
        self.dend.S_nkcc1 = 2*math.pi*self.dend.diam/2*self.dend.L+2*math.pi*(self.dend.diam/2)**2

        self.dend2.Vi_kcc2 = self.dend2.L*math.pi*(self.dend2.diam/2)**2
        self.dend2.S_kcc2 = 2*math.pi*self.dend2.diam/2*self.dend2.L+2*math.pi*(self.dend2.diam/2)**2
        self.dend2.Vi_nkcc1 = self.dend2.L*math.pi*(self.dend2.diam/2)**2
        self.dend2.S_nkcc1 = 2*math.pi*self.dend2.diam/2*self.dend2.L+2*math.pi*(self.dend2.diam/2)**2

        self.fork.Vi_kcc2 = self.fork.L*math.pi*(self.fork.diam/2)**2
        self.fork.S_kcc2 = 2*math.pi*self.fork.diam/2*self.fork.L+2*math.pi*(self.fork.diam/2)**2
        self.fork.Vi_nkcc1 = self.fork.L*math.pi*(self.fork.diam/2)**2
        self.fork.S_nkcc1 = 2*math.pi*self.fork.diam/2*self.fork.L+2*math.pi*(self.fork.diam/2)**2

        # Na-K pump parameters
        self.soma.imax_nakpump = p.soma_imax
        self.dend.imax_nakpump = p.dend_imax
        self.dend2.imax_nakpump = p.dend2_imax
        self.fork.imax_nakpump = p.fork_imax

        self.soma.km_k_nakpump = p.soma_kmk
        self.dend.km_k_nakpump = p.dend_kmk
        self.dend2.km_k_nakpump = p.dend2_kmk
        self.fork.km_k_nakpump = p.fork_kmk

        self.soma.km_na_nakpump = p.soma_kmna
        self.dend.km_na_nakpump = p.dend_kmna
        self.dend2.km_na_nakpump = p.dend2_kmna
        self.fork.km_na_nakpump = p.fork_kmna

        # Leak paramters
        self.soma.gk_leak = p.soma_gk
        self.dend.gk_leak = p.dend_gk
        self.dend2.gk_leak = p.dend2_gk
        self.fork.gk_leak = p.fork_gk

        self.soma.gna_leak = p.soma_gna
        self.dend.gna_leak = p.dend_gna
        self.dend2.gna_leak = p.dend2_gna
        self.fork.gna_leak = p.fork_gna

        self.soma.gnaother_leak = p.soma_gnaother
        self.dend.gnaother_leak = p.dend_gnaother
        self.dend2.gnaother_leak = p.dend2_gnaother
        self.fork.gnaother_leak = p.fork_gnaother

        self.soma.gcl_leak = p.soma_gcl
        self.dend.gcl_leak = p.dend_gcl
        self.dend2.gcl_leak = p.dend2_gcl
        self.fork.gcl_leak = p.fork_gcl

        # There is non specific leak channels in the leak mecanism
        # and the HH mecanism, but they are not used.
        self.soma.gfix_leak = 0 
        self.dend.gfix_leak = 0
        self.dend2.gfix_leak = 0
        self.fork.gfix_leak = 0
        self.soma.gl_hhrat = 0
        self.dend.gl_hhrat = 0
        self.dend2.gl_hhrat = 0
        self.fork.gl_hhrat = 0

        # HH parameters
        self.soma.gnabar_hhrat = p.soma_gnabar   # valeur originale : .12
        self.dend.gnabar_hhrat = p.dend_gnabar
        self.dend2.gnabar_hhrat = p.dend2_gnabar
        self.fork.gnabar_hhrat = p.fork_gnabar

        self.soma.gkbar_hhrat = p.soma_gkbar # valeur originale : .036
        self.dend.gkbar_hhrat = p.dend_gkbar
        self.dend2.gkbar_hhrat = p.dend2_gkbar
        self.fork.gkbar_hhrat = p.fork_gkbar

        self.soma.U_kcc2 = self.kcc2
        self.dend.U_kcc2 = self.kcc2 * (self.soma.Vi_kcc2 * self.dend.S_kcc2)/(self.soma.S_kcc2 * self.dend.Vi_kcc2)
        self.dend2.U_kcc2 = self.kcc2 * (self.soma.Vi_kcc2 * self.dend2.S_kcc2)/(self.soma.S_kcc2 * self.dend2.Vi_kcc2)
        self.fork.U_kcc2 = self.kcc2 * (self.soma.Vi_kcc2 * self.fork.S_kcc2)/(self.soma.S_kcc2 * self.fork.Vi_kcc2)

        self.soma.U_nkcc1 = self.nkcc1
        self.dend.U_nkcc1 = self.nkcc1 * (self.soma.Vi_nkcc1 * self.dend.S_nkcc1)/(self.soma.S_nkcc1 * self.dend.Vi_nkcc1)
        self.dend2.U_nkcc1 = self.nkcc1 * (self.soma.Vi_nkcc1 * self.dend2.S_nkcc1)/(self.soma.S_nkcc1 * self.dend2.Vi_nkcc1)
        self.fork.U_nkcc1 = self.nkcc1 * (self.soma.Vi_nkcc1 * self.fork.S_nkcc1)/(self.soma.S_nkcc1 * self.fork.Vi_nkcc1)

        # CLC-2 parameters
        self.soma.gclc2_clc2 = p.soma_gclc2
        self.dend.gclc2_clc2 = p.dend_gclc2
        self.dend2.gclc2_clc2 = p.dend2_gclc2
        self.fork.gclc2_clc2 = p.fork_gclc2



        for sec1 in self.all:
            # General cell parameters
            sec1.Ra = p.axial_resistance
            sec1.cm = p.membrane_capacitance
            sec1.clamp_iondifus = 0 # It means there is no voltage clamp initially

            # Initial intracellular concentrations
            sec1.nai = self.nai_0 # Realistics values : 5 to 30 [mM]
            sec1.ki = self.ki_0   # Realistics values : 60 to 170 [mM]
            sec1.cli = self.cli_0 # Realistics values : 8 to 30 [mM]
            sec1.cli0_iondifus = self.cli_0
            sec1.nai0_iondifus = self.nai_0
            sec1.ki0_iondifus = self.ki_0
            sec1.hco3i0_iondifus = 15 # test 15

            # Extracellular concentrations (fix values)
            sec1.hco3o0_iondifus = 26 # test 25
            sec1.nao = self.nao_0
            sec1.ko = self.ko_0 # 3.5
            sec1.clo = self.clo_0
            sec1.clo0_iondifus = self.clo_0
            sec1.nao0_iondifus = self.nao_0
            sec1.ko0_iondifus = self.ko_0

            # CLC-2 parameters
            sec1.vhalf_clc2 = p.vhalf
            sec1.vslope_clc2 = p.vslope
            sec1.ptau_clc2 = p.ptau

    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)


# Class 1 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 2 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# Main function to do the GABA puff simulation
def syna(cell, puff_pos, puff_time, puff_conc, pos, tau, dgab, rnum, clamp=False,
        clamp_amp=-70, rmp_initial=-72.38, sim_time=10000, record_pos=np.linspace(0.01,0.99,100),
        skip=0, pipette=(1500*ms, 8*mM, 140*mM, 12*mM)):
    """The function simulates a GABA puff event and records ionic concentrations, currents, conductances, reversal potential
        and other.

    Args:
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
        Tuple : Each element is another tuple containing arrays of the different recorded values.
                    Tuple[0][0-1] : Time array and membrane potential in soma array.
                    Tuple[1] : Membrane potential at each synapse positions in dendrite.
                    Tuple[2][0-3] : Arrays for the ionic concentrations at each recording point in dendrite (Cl-, K+, Na+, GABA)
                    Tuple[3][0-3] : Arrays for the ionic concentrations in soma (Cl-, K+, Na+, GABA)
                    Tuple[4][0-2] : Arrays for the reversal potentials in soma (Ecl, Ek, Ena)
                    Tuple[5][0-2] : Arrays for the reversal potentials at each synapse positions in dendrite (Ecl, Ek, Ena)
                    Tuple[6][0-4] : Arrays for the chloride currents at each synapse positions in dendrite (icl, icl_kcc2, icl_nkcc1, icl_leak, icl_gaba)
                    Tuple[7][0-5] : Arrays for the potassium currents at each synapse positions in dendrite (ik, ik_kcc2, ik_nkcc1, ik_leak, ik_nak, ik_hh)
                    Tuple[8][0-4] : Arrays for the sodium currents at each synapse positions in dendrite (ina, ina_nkcc1, ina_leak, ina_nak, ina_hh)
                    Tuple[9][0-3] : Arrays for the chloride currents in soma (icl, icl_kcc2, icl_nkcc1, icl_leak)
                    Tuple[10][0-5] : Arrays for the potassium currents in soma (ik, ik_kcc2, ik_nkcc1, ik_leak, ik_nak, ik_hh)
                    Tuple[11][0-4] : Arrays for the sodium currents in soma (ina, ina_nkcc1, ina_leak, ina_nak, ina_hh)
                    Tuple[12][0-1] : Arrays for the HCO3- current and for the total current at each synapse position in dendrite (ihco3, isynapse)
                    Tuple[13][0-5] : Arrays for the conductances at each synapse position in dendrite and 
                                     for the number of channels in the 3 different open states
                                     (gcl, ghco3, gtotal, open state 1, open state 2, open state 3) 
    """

    print(f"-------- Cellule {cell._gid} --------")
    # Creation of the necessary arrays ------------------------------
    if len(pos) == 1:
        if pos[0] < 1: position = pos
        else: position = np.linspace(0.01, 0.99, pos[0])
    else:
        position = pos # entre 0 et 1
    
    if len(rnum) == 1: rnum_ = [rnum[0] for _ in range(len(position))]
    else: rnum_ = rnum
    print("Position des synapses (0 à 1, fraction de la longueur) :\n", position)
    print("Position de la puff : ", puff_pos)
    print("Temps de la puff (ms) : ", puff_time)
    print("Concentration de la puff (mM) : ", puff_conc)
    print("Nombre de récepteurs gaba à chaque synapse : ", rnum_[0])
    print("Recording position : \n", record_pos)


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


    # GABA diffusion parameters ------------------------------------------
    h.DGab_iondifus = dgab       # Diffusion coefficient of GABA [µm2/ms]
    h.taugaba_iondifus = tau*ms  # Exchange with bath constant ([GABA]=0) [ms]
    h.fhspace_iondifus = 0.3* µm # width of the anulus in wich GABA diffuse [µm]
    print("Tau (échange avec le bain, ms) : ", tau)
    print("Coefficient de diffusion longitudinale (um2/ms) : ", dgab)


    # GABA puff event -----------------------------------------------------------------
    stim = h.NetStim()                    # Creation of the structure of the GABA puff
    stim.number = 1                       # Number of events
    stim.start = puff_time                # Time at wich the event occurs
    puff = h.gabpuff(cell.dend(puff_pos)) # Creation of the puff 

    # Makes the connection between the event and the simulation (stim and puff mecanism) 
    netcon = h.NetCon(stim, puff, 0, 0, 0, sec=cell.dend)
    netcon.weight[0] = puff_conc # Define the GABA puff concentration


    # Synapses --------------------------------------------------------------------
    gaba_R = [0]*len(position) # To stock the synapses
    for i in range(len(position)):
        gaba_R[i] = h.gaghk(cell.dend(position[i])) # Creation of the point process
        gaba_R[i].Rnumber = rnum_[i]                # Number of GABA receptors


    # Initialization at 0 of the 'messenger' concentration ----------
    # This 'messenger' makes the link between the puff (puff.mod) and
    # the extracellular GABA concentration (iondiffus.mod)
    for seg in cell.all:
        seg.messi = 0
    print(f"-------- Fin Cellule {cell._gid} --------")


    # Recording vectors for concentrations -----------------------
    dend_cli, dend_ki, dend_nai, dend_gab = [], [], [], []
    for j in record_pos:
        dend_cli.append(h.Vector().record(cell.dend(j)._ref_cli))
        dend_ki.append(h.Vector().record(cell.dend(j)._ref_ki))
        dend_gab.append(h.Vector().record(cell.dend(j)._ref_gabo))
        dend_nai.append(h.Vector().record(cell.dend(j)._ref_nai))
    soma_cli = h.Vector().record(cell.soma(0.5)._ref_cli)
    soma_nai = h.Vector().record(cell.soma(0.5)._ref_nai)
    soma_ki = h.Vector().record(cell.soma(0.5)._ref_ki)
    soma_gab = h.Vector().record(cell.soma(0.5)._ref_gabo)

    # Recording verctors for currents, reversal potentials, MP, conductance, open states in dend ---------------
    dend_icl, dend_icl_kcc2, dend_icl_nkcc1, dend_icl_leak, dend_icl_gag, dend_icl_clc2 = [], [], [], [], [], []
    dend_ik, dend_ik_kcc2, dend_ik_nkcc1, dend_ik_leak, dend_ik_nak, dend_ik_hh = [], [], [], [], [], []
    dend_ina, dend_ina_nkcc1, dend_ina_leak, dend_ina_nak, dend_ina_hh = [], [], [], [], []
    dend_ihco3, dend_igaba = [], []
    dend_surface_area = []
    dend_ecl, dend_ena, dend_ek = [], [], []
    dend_v = []
    dend_gcl, dend_ghco3, dend_ggab = [], [], []
    dend_o1, dend_o2, dend_o3 = [], [], []
    for k in range(len(position)):
        # Cloride currents, dend
        dend_icl.append(h.Vector().record(cell.dend(position[k])._ref_icl))
        dend_icl_kcc2.append(h.Vector().record(cell.dend(position[k])._ref_icl_kcc2))
        dend_icl_nkcc1.append(h.Vector().record(cell.dend(position[k])._ref_icl_nkcc1))
        dend_icl_leak.append(h.Vector().record(cell.dend(position[k])._ref_icl_leak))
        dend_icl_clc2.append(h.Vector().record(cell.dend(position[k])._ref_icl_clc2))
        dend_icl_gag.append(h.Vector().record(gaba_R[k]._ref_icl))

        # Potassium currents, dend
        dend_ik.append(h.Vector().record(cell.dend(position[k])._ref_ik))
        dend_ik_kcc2.append(h.Vector().record(cell.dend(position[k])._ref_ik_kcc2))
        dend_ik_nkcc1.append(h.Vector().record(cell.dend(position[k])._ref_ik_nkcc1))
        dend_ik_leak.append(h.Vector().record(cell.dend(position[k])._ref_ik_leak))
        dend_ik_nak.append(h.Vector().record(cell.dend(position[k])._ref_ik_nakpump))
        dend_ik_hh.append(h.Vector().record(cell.dend(position[k])._ref_ik_hhrat))

        # Sodium currents, dend
        dend_ina.append(h.Vector().record(cell.dend(position[k])._ref_ina))
        dend_ina_nkcc1.append(h.Vector().record(cell.dend(position[k])._ref_ina_nkcc1))
        dend_ina_leak.append(h.Vector().record(cell.dend(position[k])._ref_ina_leak))
        dend_ina_nak.append(h.Vector().record(cell.dend(position[k])._ref_ina_nakpump))
        dend_ina_hh.append(h.Vector().record(cell.dend(position[k])._ref_ina_hhrat))

        # HCO3 currents, dend
        dend_ihco3.append(h.Vector().record(gaba_R[k]._ref_ihco3))
        dend_igaba.append(h.Vector().record(gaba_R[k]._ref_igaba)) # igaba = ihco3 + icl

        # Area of the dendrite segment.
        # Corresponds to (segment lenght) * (circonference)
        dend_surface_area.append(cell.dend(position[k]).area()*(1e-8)) # 1e-8 -> µm2 to cm2

        # Reversal potential, dend
        dend_ecl.append(h.Vector().record(cell.dend(position[k])._ref_ecl))
        dend_ena.append(h.Vector().record(cell.dend(position[k])._ref_ena))
        dend_ek.append(h.Vector().record(cell.dend(position[k])._ref_ek))

        # Membrane potential, dend
        dend_v.append(h.Vector().record(cell.dend(position[k])._ref_v))

        # Conductance and open states, dend
        dend_gcl.append(h.Vector().record(gaba_R[k]._ref_gcl)) # µS
        dend_ghco3.append(h.Vector().record(gaba_R[k]._ref_ghco3)) # µS
        dend_ggab.append(h.Vector().record(gaba_R[k]._ref_grel)) # -
        dend_o1.append(h.Vector().record(gaba_R[k]._ref_O1))
        dend_o2.append(h.Vector().record(gaba_R[k]._ref_O2))
        dend_o3.append(h.Vector().record(gaba_R[k]._ref_O3))

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
    h.continuerun((puff_time-100)*ms)
    h.dt = p.dt3*ms
    h.continuerun(sim_time)


    # Surface of the soma ----------------------------
    # Correspond to (segment lenght) * (circonference)
    soma_surface_area = cell.soma(0.5).area()
    soma_surface_area *= (1e-8) # in cm2

    # Chloride currents in pA (1e9 for mA -> pA) ---------
    for l in range(len(dend_surface_area)):
        dend_icl[l] *= dend_surface_area[l] * (1e9)
        dend_icl_gag[l] *= 1000 # nA -> pA 
        dend_icl_kcc2[l] *= dend_surface_area[l] * (1e9)
        dend_icl_nkcc1[l] *= dend_surface_area[l] * (1e9)
        dend_icl_leak[l] *= dend_surface_area[l] * (1e9)
        dend_icl_clc2[l] *= dend_surface_area[l] * (1e9)

        dend_ik[l] *= dend_surface_area[l] * (1e9)
        dend_ik_kcc2[l] *= dend_surface_area[l] * (1e9)
        dend_ik_nkcc1[l] *= dend_surface_area[l] * (1e9)
        dend_ik_leak[l] *= dend_surface_area[l] * (1e9)
        dend_ik_nak[l] *= dend_surface_area[l] * (1e9)
        dend_ik_hh[l] *= dend_surface_area[l] * (1e9)

        dend_ina[l] *= dend_surface_area[l] * (1e9)
        dend_ina_nkcc1[l] *= dend_surface_area[l] * (1e9)
        dend_ina_leak[l] *= dend_surface_area[l] * (1e9)
        dend_ina_nak[l] *= dend_surface_area[l] * (1e9)
        dend_ina_hh[l] *= dend_surface_area[l] * (1e9)

        dend_ihco3[l] *= 1000 # nA -> pA 
        dend_igaba[l] *= 1000 # nA -> pA 


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
    result_mp_dend = dend_v
    result_con_dend = (dend_cli, dend_ki, dend_nai, dend_gab)
    result_con_soma = (soma_cli, soma_ki, soma_nai, soma_gab)
    result_e_soma = (soma_ecl, soma_ek, soma_ena)
    result_e_dend = (dend_ecl, dend_ek, dend_ena)
    result_icl_dend = (dend_icl, dend_icl_kcc2, dend_icl_nkcc1, dend_icl_leak, dend_icl_gag, dend_icl_clc2)
    result_ik_dend = (dend_ik, dend_ik_kcc2, dend_ik_nkcc1, dend_ik_leak, dend_ik_nak, dend_ik_hh)
    result_ina_dend = (dend_ina, dend_ina_nkcc1, dend_ina_leak, dend_ina_nak, dend_ina_hh)
    result_icl_soma = (soma_icl, soma_icl_kcc2, soma_icl_nkcc1, soma_icl_leak, soma_icl_clc2)
    result_ik_soma = (soma_ik, soma_ik_kcc2, soma_ik_nkcc1, soma_ik_leak, soma_ik_nak, soma_ik_hh)
    result_ina_soma = (soma_ina, soma_ina_nkcc1, soma_ina_leak, soma_ina_nak, soma_ina_hh)
    dend_iother = (dend_ihco3, dend_igaba)
    result_g = (dend_gcl, dend_ghco3, dend_ggab, dend_o1, dend_o2, dend_o3)


    return result_mp, result_mp_dend, result_con_dend, result_con_soma, result_e_soma, result_e_dend, result_icl_dend, result_ik_dend, result_ina_dend, result_icl_soma, result_ik_soma, result_ina_soma, dend_iother, result_g


# Function 2 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------


# Function 3 (begin) ------------------------------------------------------------------------------------------------------------------------------------------------
# Fonction pour sauvegarder les résultats d'une simulation
def save_sim(time_mp_soma, mp_dend,
                    dend_concentration, soma_concentration,
                    reversal_pot_soma, reversal_pot_dend,
                    dend_currents_cl, dend_currents_k, dend_currents_na,
                    soma_currents_cl, soma_currents_k, soma_currents_na,
                    dend_currents_other, dend_g_and_o, filepath,
                    dt, sim_lenght, rnum, dend_lenght, record_location,
                    rmp, nb_dend_seg, skcc2, snkcc1, cl_i, cl_o, na_i,
                    na_o, k_i, k_o, tau, time, dgab, conc, syn_pos, pos,
                    volt_clamp=False, clamp_v=False):
    """Function to save the results of a simulation done with the 'syna' function in a h5py dataset.
        It also saves the majority of the parameters of the simulation.

    Args:
        time_mp_soma (tuple): Time array and membrane potential in soma array (syna[0], lenght = 2).
        mp_dend (array): Membrane potential at each synapse positions in dendrite (syna[1], lenght = 1).
        dend_concentration (tuple): Arrays for the ionic concentrations at each recording point in dendrite (Cl-, K+, Na+, GABA - syna[2], lenght = 4)
        soma_concentration (tuple): Arrays for the ionic concentrations in soma (Cl-, K+, Na+, GABA - syna[3], lenght = 4)
        reversal_pot_soma (tuple): Arrays for the reversal potentials in soma (Ecl, Ek, Ena - syna[4], lenght = 3)
        reversal_pot_dend (tuple): Arrays for the reversal potentials at each synapse positions in dendrite (Ecl, Ek, Ena - syna[5], lenght = 3)
        dend_currents_cl (tuple): Arrays for the chloride currents at each synapse positions in dendrite (icl, icl_kcc2, icl_nkcc1, icl_leak, icl_gaba - syna[6], lenght = 5)
        dend_currents_k (tuple): Arrays for the potassium currents at each synapse positions in dendrite (ik, ik_kcc2, ik_nkcc1, ik_leak, ik_nak, ik_hh - syna[7], lenght = 6)
        dend_currents_na (tuple): Arrays for the sodium currents at each synapse positions in dendrite (ina, ina_nkcc1, ina_leak, ina_nak, ina_hh - syna[8], lenght = 5)
        soma_currents_cl (tuple): Arrays for the chloride currents in soma (icl, icl_kcc2, icl_nkcc1, icl_leak - syna[9], lenght = 4)
        soma_currents_k (tuple): Arrays for the potassium currents in soma (ik, ik_kcc2, ik_nkcc1, ik_leak, ik_nak, ik_hh - syna[10], lenght = 6)
        soma_currents_na (tuple): Arrays for the sodium currents in soma (ina, ina_nkcc1, ina_leak, ina_nak, ina_hh - syna[11], lenght = 5)
        dend_currents_other (tuple): Arrays for the HCO3- current and for the total current at each synapse position in dendrite (ihco3, isynapse - syna[12], lenght = 2)
        dend_g_and_o (tuple): Arrays for the conductances at each synapse position in dendrite and for the number of channels in the 3 different open states
                                     (gcl, ghco3, gtotal, open state 1, open state 2, open state 3 - syna[13], lenght = 6) 
        filepath (string): Saving destination.
        dt (tuple): Integration steps (dt1, dt2, dt3)
        sim_lenght (float): Lenght of the whole simulation
        rnum (float): Number of GABAA receptors
        dend_lenght (float): The lenghts of the dendrite (part 1, part 2)
        record_location (array): Position of the recording vectors for the concentration recordings.
        rmp (float): Initial resting membrane potential.
        nb_dend_seg (int): Number of segments in the dendrite (part 1, part 2)
        skcc2 (float): U_kcc2 parameter
        snkcc1 (float): U_nkcc1 parameter
        cl_i (float): Initial intracellular chloride concentration
        cl_o (float): Extracellular chloride concentration
        na_i (float): Initial intracellular sodium concentration
        na_o (float): Extracellular sodium concentration
        k_i (float): Initial intracellular potassium concentration
        k_o (float): Extracellular potassium concentration
        tau (float): Exchange with bath constant for the GABA
        time (float): Time at wich the GABA puff event occurs
        dgab (float): Diffusion coefficient of chloride
        conc (float): Concentration of the GABA puff
        syn_pos (array): Array containing the positions of the synapses
        pos (array): Position of the GABA puff event
        volt_clamp (bool, optional): True if there is the soma is voltage clamped. False if not. Defaults to False.
        clamp_v (bool, optional): Amplitude of the voltage clamp if there is one. Defaults to False.
    """

    h5_file = h5py.File(filepath, "w")

    dataset1 = h5_file.create_dataset("mp_soma", data=time_mp_soma)
    dataset1.attrs["Shape"] = f"Time array mp_soma[0]\nMP at 0.5 in soma mp_soma[1]"

    dataset2 = h5_file.create_dataset("mp_dend", data=mp_dend)
    dataset2.attrs["Shape"] = f"MP at synapses in dendrite mp_dend[0-{int(len(syn_pos)-1)}]"

    dataset3 = h5_file.create_dataset("conc_dend", data=dend_concentration)
    dataset3.attrs["temporal_resolution"] = dt
    dataset3.attrs["simulation_length"] = sim_lenght
    dataset3.attrs["number_of_receptors"] = rnum
    dataset3.attrs["dendrite_full_length_um"] = dend_lenght
    dataset3.attrs["recording_locations"] = record_location
    dataset3.attrs["voltage_clamped"] = volt_clamp
    dataset3.attrs["voltage_clamp_mV"] = clamp_v
    dataset3.attrs["initial_membrane_potential_mV"] = rmp
    dataset3.attrs["number_dendrite_segments"] = nb_dend_seg
    dataset3.attrs["kcc2_strength"] = skcc2
    dataset3.attrs["nkcc1_strength"] = snkcc1
    dataset3.attrs["cli_0"] = cl_i
    dataset3.attrs["clo_0"] = cl_o
    dataset3.attrs["nai_0"] = na_i
    dataset3.attrs["nao_0"] = na_o
    dataset3.attrs["ki_0"] = k_i
    dataset3.attrs["ko_0"] = k_o
    dataset3.attrs["puff_time"] = time
    dataset3.attrs["puff_position"] = pos
    dataset3.attrs["puff_taubath"] = tau
    dataset3.attrs["puff_longitudinal_diffusion_coefficient"] = dgab
    dataset3.attrs["puff_concentration"] = conc
    dataset3.attrs["Synapses positions"] = syn_pos
    dataset3.attrs["Shape"] = f"""Cloride concentration in dendrite arrays conc_dend[0][0-{int(len(record_location)-1)}]
    Potassium concentration in dendrite arrays conc_dend[1][0-{int(len(record_location)-1)}]
    Sodium concentration in dendrite arrays conc_dend[2][0-{int(len(record_location)-1)}]
    GABA concentration out of dendrite arrays conc_dend[3][0-{int(len(record_location)-1)}]"""

    dataset4 = h5_file.create_dataset("conc_soma", data=soma_concentration)
    dataset4.attrs["Shape"] = f"""Chloride concentration at 0.5 in soma conc_soma[0]
    Potassium concentration at 0.5 in soma conc_soma[1]
    Sodium concentration at 0.5 in soma conc_soma[2]
    External GABA concentration at 0.5 in soma conc_soma[3]"""

    dataset5 = h5_file.create_dataset("e_soma", data=reversal_pot_soma)
    dataset5.attrs["Shape"] = f"""Ecl at 0.5 in soma e_soma[0]
    Ek at 0.5 in soma e_soma[1]
    Ena at 0.5 in soma e_soma[2]"""

    dataset6 = h5_file.create_dataset("e_dend", data=reversal_pot_dend)
    dataset6.attrs["Shape"] = f"""Ecl at synapses in dendrite e_dend[0][0-{int(len(syn_pos)-1)}]
    Ek at synapses in dendrite e_dend[1][0-{int(len(syn_pos)-1)}]
    Ena at synapses in dendrite revpot[2][0-{int(len(syn_pos)-1)}]"""   

    dataset7 = h5_file.create_dataset("syn_current_cl", data=dend_currents_cl)
    dataset7.attrs["Shape"] = f"""Icl at synapses syn_current_cl[0][0-{int(len(syn_pos)-1)}]
    Icl_kcc2 at synapses syn_current_cl[1][0-{int(len(syn_pos)-1)}]
    Icl_nkcc1 at synapses syn_current_cl[2][0-{int(len(syn_pos)-1)}]
    Icl_leak at synapses syn_current_cl[3][0-{int(len(syn_pos)-1)}]
    Icl_syn at synapses syn_current_cl[4][0-{int(len(syn_pos)-1)}]
    Icl_clc2 at synapses syn_current_cl[5][0-{int(len(syn_pos)-1)}]"""   

    dataset8 = h5_file.create_dataset("syn_current_k", data=dend_currents_k)
    dataset8.attrs["Shape"] = f"""Ik at synapses syn_current_k[0][0-{int(len(syn_pos)-1)}]
    Ik_kcc2 at synapses syn_current_k[1][0-{int(len(syn_pos)-1)}]
    Ik_nkcc1 at synapses syn_current_k[2][0-{int(len(syn_pos)-1)}]
    Ik_leak at synapses syn_current_k[3][0-{int(len(syn_pos)-1)}]
    Ik_nak at synapses syn_current_k[4][0-{int(len(syn_pos)-1)}]
    Ik_hh at synapses syn_current_k[5][0-{int(len(syn_pos)-1)}]"""   

    dataset9 = h5_file.create_dataset("syn_current_na", data=dend_currents_na)
    dataset9.attrs["Shape"] = f"""Ina at synapses syn_current_na[0][0-{int(len(syn_pos)-1)}]
    Ina_nkcc1 at synapses syn_current_na[1][0-{int(len(syn_pos)-1)}]
    Ina_leak at synapses syn_current_na[2][0-{int(len(syn_pos)-1)}]
    Ina_nak at synapses syn_current_na[3][0-{int(len(syn_pos)-1)}]
    Ina_hh at synapses syn_current_na[4][0-{int(len(syn_pos)-1)}]"""   

    dataset10 = h5_file.create_dataset("soma_current_cl", data=soma_currents_cl)
    dataset10.attrs["Shape"] = f"""Icl at 0.5 in soma soma_current_cl[0]
    Icl_kcc2 at 0.5 in soma soma_current_cl[1]
    Icl_nkcc1 at 0.5 in soma soma_current_cl[2]
    Icl_leak at 0.5 in soma soma_current_cl[3]
    Icl_clc2 at 0.5 in soma soma_current_cl[4]"""   

    dataset10 = h5_file.create_dataset("soma_current_k", data=soma_currents_k)
    dataset10.attrs["Shape"] = f"""Ik at 0.5 in soma soma_current_k[0]
    Ik_kcc2 at 0.5 in soma soma_current_k[1]
    Ik_nkcc1 at 0.5 in soma soma_current_k[2]
    Ik_leak at 0.5 in soma soma_current_k[3]
    Ik_hh at 0.5 in soma soma_current_k[4]"""   

    dataset11 = h5_file.create_dataset("soma_current_na", data=soma_currents_na)
    dataset11.attrs["Shape"] = f"""Ina at 0.5 in soma soma_current_na[0]
    Ina_nkcc1 at 0.5 in soma soma_current_na[1]
    Ina_leak at 0.5 in soma soma_current_na[2]
    Ina_hh at 0.5 in soma soma_current_na[3]"""   

    dataset12 = h5_file.create_dataset("syn_current_other", data=dend_currents_other)
    dataset12.attrs["Shape"] = f"""Ihco3 at synapses syn_current_other[0][0-{int(len(syn_pos)-1)}]
    Igaba at synapses syn_current_other[1][0-{int(len(syn_pos)-1)}]"""   

    dataset13 = h5_file.create_dataset("syn_g_and_o", data=dend_g_and_o)
    dataset13.attrs["Shape"] = f"""gcl at synapses syn_g_and_o[0][0-{int(len(syn_pos)-1)}]
    ghco3 at synapses syn_g_and_o[1][0-{int(len(syn_pos)-1)}]
    grel at synapses syn_g_and_o[2][0-{int(len(syn_pos)-1)}]
    open state 1 (O1) at synapses syn_g_and_o[3][0-{int(len(syn_pos)-1)}]
    open state 2 (O2) at synapses syn_g_and_o[4][0-{int(len(syn_pos)-1)}]
    open state 3 (O3) at synapses syn_g_and_o[5][0-{int(len(syn_pos)-1)}]"""   


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
def syna_kcc2_nkcc1(first, cell, puff_pos, puff_time, puff_conc, pos, tau, dgab, rnum, clamp=False,
        clamp_amp=-70, rmp_initial=-72.38, sim_time=10000, record_pos=np.linspace(0.01,0.99,100),
        skip=0, pipette=(1500*ms, 8*mM, 140*mM, 12*mM)):
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

    # Creation des arrays nécessaires --------------------------------
    if len(pos) == 1:
        if pos[0] < 1: position = pos
        else: position = np.linspace(0.01, 0.99, pos[0])
    else:
        position = pos # entre 0 et 1

    if len(rnum) == 1: rnum_ = [rnum[0] for _ in range(len(position))]
    else: rnum_ = rnum

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


    # GABA diffusion parameters -----------------------------------------------
    h.DGab_iondifus = dgab       # Diffusion coefficient of GABA [µm2/ms]
    h.taugaba_iondifus = tau*ms  # Exchange with bath constant ([GABA]=0) [ms]
    h.fhspace_iondifus = 0.3* µm # width of the anulus in wich GABA diffuse [µm]


    # GABA puff event -----------------------------------------------------------------
    stim = h.NetStim()                    # Creation of the structure of the GABA puff
    stim.number = 1                       # Number of events
    stim.start = puff_time                # Time at wich the event occurs
    puff = h.gabpuff(cell.dend(puff_pos)) # Creation of the puff 

    # Makes the connection between the event and the simulation (stim and puff mecanism) 
    netcon = h.NetCon(stim, puff, 0, 0, 0, sec=cell.dend)
    netcon.weight[0] = puff_conc # Define the GABA puff concentration


    # Synapses --------------------------------------------------------------------
    gaba_R = [0]*len(position) # To stock the synapses
    for i in range(len(position)):
        gaba_R[i] = h.gaghk(cell.dend(position[i])) # Creation of the point process
        gaba_R[i].Rnumber = rnum_[i]                # Number of GABA receptors


    # Initialization at 0 of the 'messenger' concentration ----------
    # This 'messenger' makes the link between the puff (puff.mod) and
    # the extracellular GABA concentration (iondiffus.mod)
    for seg in cell.all:
        seg.messi = 0


    # Recording vectors for chloride concentration in dendrite --
    dend_cli = []
    for j in record_pos:
        dend_cli.append(h.Vector().record(cell.dend(j)._ref_cli))


    # Recording vectors for ionic concentrations in soma -
    soma_cli = h.Vector().record(cell.soma(0.5)._ref_cli)
    soma_nai = h.Vector().record(cell.soma(0.5)._ref_nai)
    soma_ki = h.Vector().record(cell.soma(0.5)._ref_ki)


    # Recording vectors for Ecl reversal potential and membrane potential, soma
    soma_ecl = h.Vector().record(cell.soma(0.5)._ref_ecl)
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
    h.continuerun((puff_time-100)*ms)
    h.dt = p.dt3*ms
    h.continuerun(sim_time)


    # Index of GABA puff event in the arrays ------------------
    if first == True:
        puff_ind = 0
        for i, val in enumerate(t):
            if abs(puff_time-val) < abs(puff_time-t[puff_ind]):
                puff_ind = i
    else:
        puff_ind = first
    

    # Index of the time just before the GABA puff event
    stable_conc = int(puff_ind - 2)


    # Stable values of chloride concentrations ---------------------------------------------------------
    stable_chloride = [dend_cli[i][stable_conc] for i in range(len(dend_cli))] + [soma_cli[stable_conc]]

    # Maximum delta chloride -----------------------------------------------------------------------------------------------------
    dend_cli_numpy = [np.array(dend_cli[i]) for i in range(len(dend_cli))]
    soma_cli_numpy = np.array(soma_cli)

    max_chloride_list = [max(dend_cli_numpy[i][puff_ind:]) for i in range(len(dend_cli_numpy))] + [max(soma_cli_numpy[puff_ind:])]
    min_chloride_list = [min(dend_cli_numpy[i][puff_ind:]) for i in range(len(dend_cli_numpy))] + [min(soma_cli_numpy[puff_ind:])]
    max_chloride_list = np.asarray(max_chloride_list)
    min_chloride_list = np.asarray(min_chloride_list)
    stable_chloride = np.asarray(stable_chloride)

    delta_chloride_max = np.max(max_chloride_list - stable_chloride)
    delta_chloride_min = np.min(min_chloride_list - stable_chloride)
    delta_chloride = 0
    if abs(delta_chloride_max) > abs(delta_chloride_min):
        delta_chloride = delta_chloride_max
    else:
        delta_chloride = delta_chloride_min

    return delta_chloride, soma_cli[stable_conc], soma_ki[stable_conc], soma_nai[stable_conc], soma_v[stable_conc], soma_ecl[stable_conc], soma_nai, soma_ki, soma_cli, soma_v, t


# Function 5 (_end_) ------------------------------------------------------------------------------------------------------------------------------------------------