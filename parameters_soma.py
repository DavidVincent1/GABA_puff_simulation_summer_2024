import math
from neuron.units import ms, mM


# Temporal parameters -------------------------------------------------------------------------------------------------------
time_for_stabilization = 600000  # Time for initial stabilization of the cell during simulation [ms]
simulation_lenght = 650000       # Time for the whole simulation. Must be > time_for_stabilization [ms]
dt1 = 5                          # dt1 for (0 to skip) : first stabilisation [ms]
dt2 = 0.025                      # dt2 for (skip to puff time) : just before the puff [ms]
dt3 = 0.025                      # dt3 for (puff time to simulation lenght) : real integration during the GABA puff event [ms]


# Potassic choc ---------------------------------
tchoc = 600100  # Time of the potassic event [ms]
kchoc = 13.5    # Concentration [mM]
tauchoc = 1000 # Time constant [ms]


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


# nakpump. mod parameters ---------------------------------------------------------------------
soma_imax = 0.018 # Maximum current caused by the Na-K pump in soma [mA/cm2]
soma_kmk = 2      # ? [mM]
soma_kmna = 10    # ? [mM]


# leak.mod parameters -----------------------------------------------------------------------
soma_gk = 9e-5           # Potassium leak channels conductance in soma [S/cm2]
soma_gna = 1.8e-5        # Sodium leak channels conductance in soma [S/cm2]
soma_gnaother = 0.9e-5   # Other Sodium channels conductance in soma [S/cm2]
soma_gcl = 0.36e-5       # Chloride leak channels conductance in soma [S/cm2]


# hh_rat.mod -----------------------------------------------------------------------------------------------
factor = 1/100 # Reduction factor of the HH channels strength

soma_gnabar = 0.2162*factor  # Maximum HH channels conductance for sodium in soma [S/cm2]
soma_gkbar = 0.06486*factor  # Maximum HH channels conductance for potassium in soma [S/cm2]


# clc2.mod ----------------------------------------------------------------------------------
soma_gclc2 = 1.8e-5  # Maximum conductance of CLC-2 channels in soma [S/cm2] 

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