import math
from neuron.units import mV, µm, mM, ms


# Temporal parameters -------------------------------------------------------------------------------------------------------
time_for_stabilization = 600000  # Time for initial stabilization of the cell during simulation [ms]
simulation_lenght = 650000       # Time for the whole simulation. Must be > time_for_stabilization [ms]
dt1 = 5                          # dt1 for (0 to skip) : first stabilisation [ms]
dt2 = 0.144                      # dt2 for (skip to puff time) : just before the puff [ms]
dt3 = 0.144                      # dt3 for (puff time to simulation lenght) : real integration during the GABA puff event [ms]


# GABA puff parameters ----------------------------------------------------------------------------------------------------------------------------
time_of_puff = time_for_stabilization + 100 # Time of the GABA puff event. Must be time_for_stabilization < time_of_puff < simulation lenght  [ms]
position_of_puff = 40                       # Position of the GABA puff event [µm]
concentration_of_puff = 1                   # Concentration of the GABA puff [mM]
Dgaba = 0.765                               # GABA diffusion coefficient [um2/ms]
                                                # Source :
                                                # The structure and diffusion behaviour of the neurotransmitter  
                                                # c-aminobutyric acid (GABA) in neutral aqueous solutions
tau_GABA = 2100                             # GABA exchange with bath constant [ms]
number_of_rec = 40                          # Number of recording vectors for concentrations


# Synapses parameters ------------------------------------------------------------------------------------
dend_lenght_with_synapses = 80 # Dendrite (first part) lenght where there is synapses (< dend_lenght) [µm]
syn_per_micron = 1/4           # Synapse per µm ratio
rnum = 200                     # Number of GABA_A neuroreceptors per synapse


# Voltage clamp parameters ----------------------------------------------------------------------------------------------
clamp = True                            # True : the soma is voltage clamped and False : the soma is not voltage clamped.
clamp_amp = -50                         # Amplitude of the clamp [mV]
pipett = (1500*ms, 8*mM, 140*mM, 12*mM) # Pipette parameters (if the soma is voltage clamped)
                                            # pipett[0] : Exchange constant with the pipette
                                            # pipett[1] : Chloride concentration in the pipette
                                            # pipett[2] : Potassium concentration in the pipette
                                            # pipett[3] : Sodium concentration in the pipette


# Cell geometry -> [soma]>[dend]>[dend2] -----------------------------
soma_lenght = 20   # Soma lenght [µm]
soma_diam = 20     # Soma diameter [µm]
soma_nseg = 1      # Number of segments in soma [-]

dend_lenght = 150  # Dendrite (first part) lenght [µm]
dend_diam = 1      # Dendrite (first part) diameter [µm]
dend_nseg = 251    # Number of segments in dendrite (first part) [-]

dend2_lenght = 330 # Dendrite (second part) lenght [µm]
dend2_diam = 1     # Dendrite (secpnd part) diameter [µm]
dend2_nseg = 21    # Number of segments in dendrite (second part) [-]

fork = False                    # True if you want the cell to have a fork [-]
fork_position = 100/dend_lenght # Position of the fork on the dendrite (first part, between 0 and 1) [-]
fork_lenght = 100               # Fork length [µm]
fork_diam = 1                   # Fork diameter [µm]
fork_nseg = 6                   # Number of segments in fork [-]


# kcc2.mod and nkcc1.mod parameters -------------------------------------------------------------------
# These mecanisms are dependant of the volume and the area of the section, but the value
# entered in each of these sections is ajusted so that the maximum current density created by
# the mecanisms is uniform on the cell. The parameters here are the value for the soma.

U_kcc2 = 1e-4     # Maximum KCC2 pump strength [mM/ms]
U_nkcc1 = 1e-5    # Maximum NKCC1 pump strength [mM/ms]

V = soma_lenght*math.pi*(soma_diam/2)**2                             # Volume of soma [um3]
S = 2*math.pi*(soma_diam/2)*soma_lenght + 2*math.pi*(soma_diam/2)**2 # Surface of soma [um2]
F = 96485.309                                                        # Faraday constant [faraday]
print("Current density (KCC2) [mA/cm2] : ", (U_kcc2*F*V)/(S*1e4))    # Current density caused by KCC2
print("Current density (NKCC1) [mA/cm2] : ", (U_nkcc1*F*V)/(S*1e4))  # Current density caused by NKCC1


# iondif.mod parameters ---------------------------------------------------------
soma_DCl = 2  # Chloride diffusion coefficient in soma [um2/ms]
dend_DCl = 2  # Chloride diffusion coefficient in dendrite (first part) [um2/ms]
dend2_DCl = 2 # Chloride diffusion coefficient in dendrite (second part) [um2/ms]
fork_DCl = 2  # Chloride diffusion coefficient in fork [um2/ms]


# nakpump. mod parameters ---------------------------------------------------------------------
soma_imax = 0.010  # Maximum current caused by the Na-K pump in soma [mA/cm2]
dend_imax = 0.010  # Maximum current caused by the Na-K pump in dendrite (fisrt part) [mA/cm2]
dend2_imax = 0.010 # Maximum current caused by the Na-K pump in dendrite (second part) [mA/cm2]
fork_imax = 0.010  # Maximum current caused by the Na-K pump in fork [mA/cm2]

soma_kmk = 2  # ? [mM]
dend_kmk = 2  # ? [mM]
dend2_kmk = 2 # ? [mM]
fork_kmk = 2  # ? [mM]

soma_kmna = 10  # ? [mM]
dend_kmna = 10  # ? [mM]
dend2_kmna = 10 # ? [mM]
fork_kmna = 10  # ? [mM]


# leak.mod parameters -----------------------------------------------------------------------
soma_gk = 5e-5  # Potassium leak channels conductance in soma [S/cm2]
dend_gk = 5e-5  # Potassium leak channels conductance in dendrite (first part) [S/cm2]
dend2_gk = 5e-5 # Potassium leak channels conductance in dendrite (second part) [S/cm2]
fork_gk = 5e-5  # Potassium leak channels conductance in fork [S/cm2]

soma_gna = 1e-5  # Sodium leak channels conductance in soma [S/cm2]
dend_gna = 1e-5  # Sodium leak channels conductance in dendrite (first part) [S/cm2]
dend2_gna = 1e-5 # Sodium leak channels conductance in dendrite (second part) [S/cm2]
fork_gna = 1e-5  # Sodium leak channels conductance in fork [S/cm2]

soma_gnaother = 0.5e-5  # Other Sodium channels conductance in soma [S/cm2]
dend_gnaother = 0.5e-5  # Other Sodium channels conductance in dendrite (first part) [S/cm2]
dend2_gnaother = 0.5e-5 # Other Sodium channels conductance in dendrite (second part) [S/cm2]
fork_gnaother = 0.5e-5  # Other Sodium channels conductance in fork [S/cm2]

soma_gcl = 0.2e-5  # Chloride leak channels conductance in soma [S/cm2]
dend_gcl = 0.2e-5  # Chloride leak channels conductance in dendrite (first part) [S/cm2]
dend2_gcl = 0.2e-5 # Chloride leak channels conductance in dendrite (second part) [S/cm2]
fork_gcl = 0.2e-5  # Chloride leak channels conductance in fork [S/cm2]


# hh_rat.mod -----------------------------------------------------------------------------------------------
factor = 1/100 # Reduction factor of the HH channels strength

soma_gnabar = 0.12*factor  # Maximum HH channels conductance for sodium in soma [S/cm2]
dend_gnabar = 0.12*factor  # Maximum HH channels conductance for sodium in dendrite (first part) [S/cm2]
dend2_gnabar = 0.12*factor # Maximum HH channels conductance for sodium in dendrite (second part) [S/cm2]
fork_gnabar = 0.12*factor  # Maximum HH channels conductance for sodium in fork [S/cm2]

soma_gkbar = 0.036*factor  # Maximum HH channels conductance for potassium in soma [S/cm2]
dend_gkbar = 0.036*factor  # Maximum HH channels conductance for potassium in dendrite (first part) [S/cm2]
dend2_gkbar = 0.036*factor # Maximum HH channels conductance for potassium in dendrite (second part) [S/cm2]
fork_gkbar = 0.036*factor  # Maximum HH channels conductance for potassium in fork [S/cm2]


# clc2.mod ----------------------------------------------------------------------------------
soma_gclc2 = 5e-5  # Maximum conductance of CLC-2 channels in soma [S/cm2] 
dend_gclc2 = 5e-5  # Maximum conductance of CLC-2 channels in dendrite (first part) [S/cm2] 
dend2_gclc2 = 5e-5 # Maximum conductance of CLC-2 channels in dendrite (second part) [S/cm2]
fork_gclc2 = 5e-5  # Maximum conductance of CLC-2 channels in fork [S/cm2]

ptau = 300         # Parameter in the equation of the opened probabilities [ms] 
vhalf = 15         # Parameter in the equation of the opened probabilities [mV]
vslope = -14       # Parameter in the equation of the opened probabilities [mV]


# Other parameters ----------------------------------------------------------------
axial_resistance = 100   # Axial resistance in all the cell [Ohm * cm]
membrane_capacitance = 1 # Membrane capacitance [micro Farads/cm2]
temperature = 23         # Temperature during the simulation [degC]

intial_cli = 6.5844  # Initial intracellular chloride concentration in cell [mM]
initial_ki = 115.04  # Initial intracellular potassium concentration in cell [mM]
initial_nai = 13.530 # Initial intracellular sodium concentration in cell [mM]
clo = 130.5          # Extracellular chloride concentration [mM]
ko = 3.5             # Extracellular potassium concentration [mM]
nao = 147.25         # Extracellular sodium concentration [mM]


# Fix points --------------------------------------------
# U_kcc2 = 1e-6     # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1.5e-4  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Unclamped
cli, ki, nai = 27.731, 179.972, 12.939

# U_kcc2 = 1e-6     # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1.5e-4  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Clamped at -90 mV
cli, ki, nai = 21.31797, 173.96972, 15.3054648

# U_kcc2 = 1e-4   # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1e-5  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Unclamped
cli, ki, nai = 4.90666, 140.033, 14.8446

# U_kcc2 = 1e-4   # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1e-5  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Clamped at -80 mV
cli, ki, nai = 5.058, 142.89, 14.4018

# U_kcc2 = 1e-4   # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1e-5  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Clamped at -40 mV
cli, ki, nai = 7.135, 107.62, 13.468

# U_kcc2 = 1e-4   # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1e-5  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 1e-5
# Clamped at -50 mV
cli, ki, nai = 6.5844, 115.04, 13.530

# U_kcc2 = 1e-4   # Maximum KCC2 pump strength [mM/ms]
# U_nkcc1 = 1e-5  # Maximum NKCC1 pump strength [mM/ms]
# Leak cl = 0.2e-5
# Clamped at -50 mV
#cli, ki, nai = 6.5844, 115.04, 13.530