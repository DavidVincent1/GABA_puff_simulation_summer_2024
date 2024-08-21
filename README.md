# GABA_puff_simulation_summer_2024



## How to start
1) Download NEURON (https://neuron.yale.edu/neuron/download)

2) Download the mod files which are in the 'mod_files' folder. To do so, go to the NEURON folder and open the 'mknrndll' file. Then, choose the destination of the folder containing the mod files and click on 'Make nrnmech.dll'.

3) Read the rest of this README, the 'Equation_description.pdf' file and the 'fichiers_mod.pdf'.

4) Use 'GABA_puff_simulation.py' or 'kcc2_nkcc1_grid_analysis.py' to run simulation.


## Experimental_data folder
This folder contains the experimental data (from Justin Hamel) analyzed during the summer of 2024. It also contains scripts used to do the analysis. 


### lifetime_to_conc_approx.py
Input 1: Excel sheet containing lifetime measurements in different segments of a dendrite.
Input 2: Simulation results (h5py dataset).

This script calculates an approximation of the chloride concentration. Both simulation results and experimental data are normalized by the maximum and the minimum of the combined arrays and are compared.

Output: Multiple graphs.


### Multiple_tif_to_stack.py
Input: Path of a folder containing multiple tiff files of a time series.

The script combines the tiff files into a stack. The stack is saved in the input folder.

Output: The saved stack.


### tau_different_pixel_size.py
Input: Path of a tiff file. The tiff file corresponds to a stack showing the photon counts associated with the GABA during an experiment.

Monoexponential and biexponential fits are first done on the sum of all photons counts in the data. The script then creates two new arrays corresponding to the original data with a modified pixel size. One pixel of the first new array corresponds to a 5x5 pixels region of the original data. One pixel of the second new array corresponds to a 10x10 pixels region of the original data. Then, a threshold is applied to the data so that only the pixels containing a certain number of photon counts are conserved for further analysis. By doing so, the pixels containing too much noise and the pixels which contain too much photon counts due to cross talk are excluded. Then, a curve fit is done on the remaining pixels and the calculated lifetimes are extracted. Multiple fiting parameters and results are printed to the terminal.

Output 1: Graph showing the total number of counts in each pixel of the original data. (figure : Image)
Output 2: Graph showing the data with new pixel sizes. (figure with two panels : Pixel_reduction)
Output 3: Graph showing which pixels were used to fit. The yellow pixels have a value of one and were used and the other pixels have a value of 0 and were not used. (figure : Fited_pixels)
Output 4: Graph showing the calculated lifetime per pixel for the data with new pixel sizes. (Fited_pixels_tau)
Output 5: Graph showing the monoexponential and biexponential fits on all of the photon counts. (Curve_fit_all_counts)
Output 6: Graph showing all the curve fits done on pixels of the data with new pixel sizes. (Fited_pixels_all_curves)
Output 7: Graph showing photon counts and curve fits done on pixels near the GABA puff and far from the GABA puff. (Fit_near_and_far_from_puff)  


## mod_files folder
The folder contains mod files essential to the model. See 'fichiers_mod.pdf' for more details on the way these files work. See 'Equation_description.pdf' for details on the mathematical equations. 


## The other scripts


### parameters.py
Every parameter of the simulation are here. This file is used as an input for the simulation files.


### function.py
The functions and classes used by other scripts are here.


### GABA_puff_simulation.py
Input : See 'parameters.py'.

Run a GABA puff simulation. The recorded values during all simulation length are : the intracellular concentrations (chloride, potassium and sodium), the extracellular concentration of GABA, the different chloride currents, the different potassium currents, the different sodium currents, the HCO3- and synapse currents, the reversal potentials, the memrane potential, the conductances at synapses and the number of GABA open channels at synapses (in different states). See 'syn' function in the 'function.py'.

Output: HDF5 dataset saved in the dataset folder.


### data_reader_animation_mega_sim.py
Input: Simulation results (HDF5 dataset).

The script generates animations. At the beginning of the file, the desired animation must be chosen. It must also be specified if the animation is saved. The possible animations are : chloride intracellular concentration and GABA extracelular function as a function of the position on the dendrite, chloride currents at each synapses as a function of the position on the dendrite, the number of open GABA channels at synapses as a function of the position on the dendrite and a combination of the three previous animations.

Output: The chosen animations are shown and saved (if you choose to).


### data_reader_mega_sim
Input: Simulation results (HDF5 dataset).

Generates figures. At the beginning of the file, the desired figures must be chosen. Be careful not to choose too many figures at the same time. This can cause memory issues and crash.

Output: The chosen figures.


### kcc2_nkcc1_grid_analysis.py
Input 1: The path of the folder where you want the pickle files to be saved.
Input 2: The dimension of the grid (two integers -> 'nb_kcc2' and 'nb_nkcc1')
Input 3: The range of values for 'U_kcc2' and 'U_nkcc1' parameters (two tuples -> 'range_kcc2' and 'range_nkcc1')

Run multiple simulations with different kcc2 and nkcc1 strengths. It only keeps certain recorded values including the intracellular concentrations just before the GABA puff, the membrane potential just before the GABA puff and the reversal potential of chloride just before GABA puff. It also calculates the maximum chloride intracellular concentration variation due to the GABA puff. The script uses multiprocessing to run parallel simulations. This can slow down your computer. Also, all the graphs are saved as pickle files in the input folder. This makes it possible to edit them through the 'data_reader_pickle.py' file.

Output 1: Grid of the chloride intracellular concentrations before the GABA puff for all the simulations (figure : Stable chloride)
Output 2: Grid of the potassium intracellular concentrations before the GABA puff for all the simulations (figure : Stable potassium)
Output 3: Grid of the sodium intracellular concentrations before the GABA puff for all the simulations (figure : Stable sodium)
Output 4: Grid of the membrane potential before the GABA puff for all the simulations (figure : Stable mp)
Output 5: Grid of the chloride reversal potential before the GABA puff for all the simulations (figure : Stable Ecl)
Output 6: Grid of the maximum chloride intracellular concentration variation due to the GABA puff for all the simulations (figure : Delta chloride)
Output 7: Chloride concentration in soma as a function of time for all simulations (figure : Chloride steady state)
Output 8: Potassium concentration in soma as a function of time for all simulations (figure : Potassium steady state)
Output 9: Sodium concentration in soma as a function of time for all simulations (figure : Sodium steady state)
Output 10: Membrane potential in soma as a function of time for all simulations (figure : Resting mp)
Output 11: All graphs are saved as pickle file in the input folder.

### data_reader_pickle.py
Input: Path of the folder containing pickle files from a run of 'kcc2_nkcc1_grid_analysis.py'.

Open the graphs in the folder. This makes it possible to edit them if you want to change something.

Outputs: See Output 1 to Output 10 in 'kcc2_nkcc1_grid_analysis.py'.