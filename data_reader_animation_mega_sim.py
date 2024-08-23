import h5py
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, FFMpegWriter
from function import show_info_sim


# Graph parameters -----------------------------
rcParams.update({'font.size': 22})
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'


# Loading the h5py dataset ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# INPUT
#f = h5py.file(r"Path of the dataset file", 'r')
path = r"Path of the dataset file"
f = h5py.File(path, 'r')


# To see the keys of the datatset
#print(list(f.keys()))


# Offset to skip the initial stabilization of the data and language used on the graph. ----------
decal = 500000    # Th animation will be from decal to simulation lenght - end skip
end_skip = 70000  # in ms
skip = 10         # The animation will only take 1/skip data for the animation
skip2 = 150       # The animation will only take 1/skip2 data for the animation (toward the end)
time_dt1 = 500000 # Duration of the simulation at dt1 in ms
time_dt2 = 100    # Duration of the simulation at dt2 in ms
time_dt3 = 89900  # Duration of the simulation at dt3 in ms


# Choose the animation you want -------------------------------------------------------------------------------------
show_info = 0               # Show information on the simulation
anim_chlore = 0             # Animation for chloride and GABA concentration as a function of the position on dendrite
save_anim_chlore = 0        # If you want to save the animation (can take time)
anim_icl = 0                # Animation for chloride currents as a function of the position on dendrite (2 panels)
save_anim_icl = 0           # If you want to save the animation (can take time)
anim_open_channels = 1      # Animation for the currents created by the synapses
save_anim_open_channels = 1 # If you want to save the animation (can take time)
anim_all_four = 0           # All the above anumation in one figure with four panels
save_anim_all_four = 0      # If you want to save the anumation (can take time)


# Separation of the dataset in arrays ----------------------------------------------------
t_and_soma_v = f["mp_soma"] # Time and MP in soma
#print(len(t_and_soma_v[0]))

dend_v = f["mp_dend"]                 # MP on dendrite (at synapses)
dend_conc = f["conc_dend"]            # Concentrations in dendrite (at recording pos)
dend_e = f["e_dend"]                  # Reversal potentials in dendrite (at synapses)
dend_current_cl = f["syn_current_cl"] # Icl in dendrite (at synapses)
dend_current_k = f["syn_current_k"]   # Ik in dendrite (at synapses)
dend_current_na = f["syn_current_na"] # Ina in dendrite (at synapses)

soma_conc = f["conc_soma"]             # Concentrations in soma
soma_e = f["e_soma"]                   # Reversal potentials in soma
soma_current_cl = f["soma_current_cl"] # Icl in soma
soma_current_k = f["soma_current_k"]   # Ik in soma
soma_current_na = f["soma_current_na"] # Ina in soma

synapse_current = f["syn_current_other"] # Ihco3 and Igaba (Itot) in dendrite (at synapses)
gcl_and_other = f["syn_g_and_o"]         # gcl and ghco3 in dendrite (at synapses)


# Show info on the simulation and on the arrays -----------------------
if show_info == 1:
    show_info_sim(time_mp_soma=t_and_soma_v,
                                    mp_dend=dend_v,
                                    dend_concentration=dend_conc,
                                    soma_concentration=soma_conc,
                                    reversal_pot_soma=soma_e,
                                    reversal_pot_dend=dend_e,
                                    dend_currents_cl=dend_current_cl,
                                    dend_currents_k=dend_current_k,
                                    dend_currents_na=dend_current_na,
                                    soma_currents_cl=soma_current_cl,
                                    soma_currents_k=soma_current_k,
                                    soma_currents_na=soma_current_na,
                                    dend_currents_other=synapse_current,
                                    dend_g=gcl_and_other)


# Other values used later -----------------------------------------------------------------------------
record_pos = dend_conc.attrs["recording_locations"]    # Recording positions for the concentrations [-]
syn_pos = dend_conc.attrs["Synapses positions"]        # Synapses positions [-]
center_ind = int(len(syn_pos)/2-1)                     # Central synapse index [-]
lenght = dend_conc.attrs["dendrite_full_length_um"][0] # Dendrite (first part) lenght [µm]
puff_time = dend_conc.attrs["puff_time"]               # Time at wich the GABA puff occurs [ms]
puff_pos = dend_conc.attrs["puff_position"]            # Position of the GABA puff [-]
sim_time = dend_conc.attrs["simulation_length"]        # Simulation lenght [ms]
rnum = dend_conc.attrs["number_of_receptors"]          # Number of GABA receptors per synapse [-]
kcc2 = dend_conc.attrs["kcc2_strength"]                # U_KCC2 parameter [mM/ms]
nkcc1 = dend_conc.attrs["nkcc1_strength"]              # U_NKCC1 parameter [mM/ms]
dt = dend_conc.attrs["temporal_resolution"]            # Temporal resolution [ms]


# Everything separated -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
t_ = t_and_soma_v[0]  # Time array [ms]
t = (t_ - decal)/1000 # Time array with shifted zero [s]

# Concentration arrays in dendrite and in soma
dend_cli, dend_gabo = dend_conc[0], dend_conc[3]

# Chloride currents arrays in dendrite and in soma
dend_icl, dend_icl_kcc2, dend_icl_nkcc1, dend_icl_leak, dend_icl_synapses = dend_current_cl[0], dend_current_cl[1], dend_current_cl[2], dend_current_cl[3], dend_current_cl[4]

# HCO3 and total currents arrays in dendrite (at synapses)
dend_ihco3, dend_igaba = synapse_current[0], synapse_current[1]

# Chloride and HCO3 conductance at synapses. Open states of the GABA channels.
open_states1, open_states2, open_states3 = gcl_and_other[3], gcl_and_other[4], gcl_and_other[5] 
open_states = open_states1 + open_states2 + open_states3

del t_and_soma_v
del dend_v
del dend_conc
del dend_e
del dend_current_cl
del dend_current_k
del dend_current_na
del soma_conc
del soma_e
del soma_current_cl
del soma_current_k
del soma_current_na
del synapse_current
del gcl_and_other



if anim_chlore == 1:
    # Beginning and end index of the range of the simulation -------------------
    begin = int(decal/dt[0])
    end = int(time_dt1/dt[0] + time_dt2/dt[1] + time_dt3/dt[2] - end_skip/dt[2])
    puff_ind = int(time_dt1/dt[0] + time_dt2/dt[1])
    range_high_precision = [5,100]


    # Creation of the arrays for the animation ------------------------------------------------------------
    # The arrays are separated into three region
    # The first one has a low temporal resolution but the third one has an ever lower resolution.
    # The second region has a higher temporal resolution and is to capture the GABA puff correctly  
    dend_cli1 = dend_cli[:, begin:int(puff_ind-range_high_precision[0])]
    dend_cli2 = dend_cli[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_cli3 = dend_cli[:, int(puff_ind+range_high_precision[1]):end]
    dend_anime1 = dend_cli1[:, ::skip]
    dend_anime2 = dend_cli2
    dend_anime3 = dend_cli3[:, ::skip2]

    dend_gabo1 = dend_gabo[:, begin:int(puff_ind-range_high_precision[0])]
    dend_gabo2 = dend_gabo[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_gabo3 = dend_gabo[:, int(puff_ind+range_high_precision[1]):end]
    gabo_anime1 = dend_gabo1[:, ::skip]
    gabo_anime2 = dend_gabo2
    gabo_anime3 = dend_gabo3[:, ::skip2]

    t1 = t[begin:int(puff_ind-range_high_precision[0])]
    t2 = t[int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    t3 = t[int(puff_ind+range_high_precision[1]):end]
    tanime1 = t1[::skip]
    tanime2 = t2
    tanime3 = t3[::skip2]

    dend_anime = np.concatenate((dend_anime1, dend_anime2, dend_anime3), axis=1)
    gabo_anime = np.concatenate((gabo_anime1, gabo_anime2, gabo_anime3), axis=1)
    new_t =np.concatenate((tanime1, tanime2, tanime3))


    # Array for the label -------------------------------------
    for_label_1 = np.full(dend_anime[0].shape, fill_value=-100)


    # Figure object creation ---------------------------------------------------------------
    fig, ax = plt.subplots()
    ax_gaba = ax.twinx()                              # Second ax for the GABA concentration
    line, = ax.plot([], [], lw=2, color='orange')     # Chloride concentration
    line2, = ax_gaba.plot([], [], lw=2, color='blue') # GABA concentration


    # Curve drawned for the legend only ----------------------------------------
    ax.plot(new_t, for_label_1, color='orange', label=r'$[Cl^-]_i$')
    ax.plot(new_t, for_label_1, color='blue', label=r'$[GABA]_o$')
    ax.plot(new_t, for_label_1, color='black', linestyle='--', label='Synapses')
    ax.legend(loc="upper right")


    # Axis limits and box for the time -----------------------------------------------
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.95)
    text = ax.text(0.02, 0.93, f'', color='black', transform=ax.transAxes, bbox=props)
    ax.set_xlim(0, record_pos[-1]*lenght)
    ax.set_ylim(np.min(dend_anime), np.max(dend_anime))
    ax_gaba.set_ylim(np.min(gabo_anime), np.max(gabo_anime))


    # Axis label ------------------------------------------------
    ax.set_xlabel("Position on dendrite [µm]")
    ax.set_ylabel("Chloride concentration [mM]")
    ax_gaba.set_ylabel("GABA extracellular concentration [mM]")
    for j in syn_pos:
        ax.axvline(j*lenght, color='black', lw=1, linestyle='--')


    # Initialization function for FuncAnimation
    def init():
        line.set_data([], [])
        line2.set_data([], [])
        text.set_text("")
        return line, line2, text


    # Function that is called at each iteration by FuncAnimation
    def animate(i):
        text.set_text(f'Time : {new_t[i]:.2f} s')
        x = record_pos*lenght
        y = dend_anime[:,i]
        line.set_data(x, y)

        x2 = record_pos*lenght
        y2 = gabo_anime[:,i]
        line2.set_data(x2, y2)
        return line, line2, text


    # Animation creation --------------------------------------------------------------------------------
    anim = FuncAnimation(fig, animate, init_func=init, frames=len(dend_anime[0]), interval=50, blit=True)


    # Show the animation before saving ------------------------------------------------------------
    plt.show()
    if save_anim_chlore == 1:
        writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('animation' + f'\Cl_{len(syn_pos)}synapses_kcc2={kcc2}_nkcc1={nkcc1}_-90mV.mp4', writer='ffmpeg')



if anim_icl == 1:
    # Beginning and end index of the range of the simulation -------------------
    begin = int(decal/dt[0])
    end = int(time_dt1/dt[0] + time_dt2/dt[1] + time_dt3/dt[2] - end_skip/dt[2])
    puff_ind = int(time_dt1/dt[0] + time_dt2/dt[1])
    range_high_precision = [5,100]


    # Creation of the arrays for the animation ------------------------------------------------------------------------
    # The arrays are separated into three region
    # The first one has a low temporal resolution but the third one has an ever lower resolution.
    # The second region has a higher temporal resolution and is to capture the GABA puff changes correctly 
    dend_icl1 = dend_icl[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl2 = dend_icl[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl3 = dend_icl[:, int(puff_ind+range_high_precision[1]):end]
    icl_anime1 = dend_icl1[:, ::skip]
    icl_anime2 = dend_icl2
    icl_anime3 = dend_icl3[:, ::skip2]

    dend_icl_syn1 = dend_icl_synapses[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_syn2 = dend_icl_synapses[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_syn3 = dend_icl_synapses[:, int(puff_ind+range_high_precision[1]):end]
    icl_syn_anime1 = dend_icl_syn1[:, ::skip]
    icl_syn_anime2 = dend_icl_syn2
    icl_syn_anime3 = dend_icl_syn3[:, ::skip2]

    dend_icl_kcc21 = dend_icl_kcc2[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_kcc22 = dend_icl_kcc2[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_kcc23 = dend_icl_kcc2[:, int(puff_ind+range_high_precision[1]):end]
    icl_kcc2_anime1 = dend_icl_kcc21[:, ::skip]
    icl_kcc2_anime2 = dend_icl_kcc22
    icl_kcc2_anime3 = dend_icl_kcc23[:, ::skip2]

    dend_icl_nkcc11 = dend_icl_nkcc1[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_nkcc12 = dend_icl_nkcc1[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_nkcc13 = dend_icl_nkcc1[:, int(puff_ind+range_high_precision[1]):end]
    icl_nkcc1_anime1 = dend_icl_nkcc11[:, ::skip]
    icl_nkcc1_anime2 = dend_icl_nkcc12
    icl_nkcc1_anime3 = dend_icl_nkcc13[:, ::skip2]

    dend_icl_leak1 = dend_icl_leak[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_leak2 = dend_icl_leak[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_leak3 = dend_icl_leak[:, int(puff_ind+range_high_precision[1]):end]
    icl_leak_anime1 = dend_icl_leak1[:, ::skip]
    icl_leak_anime2 = dend_icl_leak2
    icl_leak_anime3 = dend_icl_leak3[:, ::skip2]

    t1 = t[begin:int(puff_ind-range_high_precision[0])]
    t2 = t[int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    t3 = t[int(puff_ind+range_high_precision[1]):end]
    tanime1 = t1[::skip]
    tanime2 = t2
    tanime3 = t3[::skip2]

    icl_anim = np.concatenate((icl_anime1, icl_anime2, icl_anime3), axis=1)
    icl_syn_anim = np.concatenate((icl_syn_anime1, icl_syn_anime2, icl_syn_anime3), axis=1)
    icl_kcc2_anim = np.concatenate((icl_kcc2_anime1, icl_kcc2_anime2, icl_kcc2_anime3), axis=1)
    icl_nkcc1_anim = np.concatenate((icl_nkcc1_anime1, icl_nkcc1_anime2, icl_nkcc1_anime3), axis=1)
    icl_leak_anim = np.concatenate((icl_leak_anime1, icl_leak_anime2, icl_leak_anime3), axis=1)
    new_t =np.concatenate((tanime1, tanime2, tanime3))


    # Array for the label --------------------------------------
    for_label_1 = np.full(icl_anim[0].shape, fill_value=-100000)


    # Figure object creation ---------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1)
    scat1 = ax1.scatter([], [], lw=2, color='orange', s=60)
    scat2 = ax1.scatter([], [], lw=2, facecolor='none', edgecolor='darkgreen', s=50)
    scat3 = ax2.scatter([], [], lw=2, facecolor='none', edgecolor='red', s=50)
    scat4 = ax2.scatter([], [], lw=2, facecolor='none', edgecolor='deepskyblue', s=50)
    scat5 = ax2.scatter([], [], lw=2, facecolor='none', edgecolor='orchid', s=50)



    # Curve drawned for the legend only ----------------------------------------------------------------
    ax1.plot(new_t, for_label_1, color='black', label=r'$Synapses$', linestyle='--')
    ax1.scatter(new_t, for_label_1, color='orange', label=r'$I_{cl}$')
    ax1.scatter(new_t, for_label_1, facecolor='none', edgecolor='darkgreen', label=r'$I_{cl,synapse}$')
    ax2.scatter(new_t, for_label_1, facecolor='none', edgecolor='red', label=r'$I_{cl,kcc2}$')
    ax2.scatter(new_t, for_label_1, facecolor='none', edgecolor='deepskyblue', label=r'$I_{cl,nkcc1}$')
    ax2.scatter(new_t, for_label_1, facecolor='none', edgecolor='orchid', label=r'$I_{cl,leak}$')


    # Axis limits and box for the time ----------------------------------------------------------------------------------------------------------------------------------------
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.95)
    text = ax1.text(0.02, 0.85, f'', color='black', transform=ax1.transAxes, bbox=props)
    ax1_min, ax1_max = min([np.min(icl_anim), np.min(icl_syn_anim)]), max([np.max(icl_syn_anim), np.max(icl_anim)])
    ax2_min, ax2_max = min([np.min(icl_kcc2_anim), np.min(icl_nkcc1_anim), np.min(icl_leak_anim)]), max([np.max(icl_kcc2_anim), np.max(icl_nkcc1_anim), np.max(icl_leak_anim)])
    ax1.set_ylim(ax1_min, ax1_max)
    ax2.set_ylim(ax2_min, ax2_max)
    ax1.set_xlabel("Position on dendrite [µm]")
    ax1.set_ylabel(r"Currents [pA]")
    ax2.set_xlabel("Position on dendrite [µm]")
    ax2.set_ylabel(r"Currents [pA]")
    if abs(ax1_min) > abs(ax1_max):
        ax1.legend(loc="lower right")
    else:
        ax1.legend(loc="upper right")
    if abs(ax2_min) > abs(ax2_max):
        ax2.legend(loc="lower right")
    else:
        ax2.legend(loc="upper right")

    for j in syn_pos:
        ax1.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax2.axvline(j * lenght, color='black', lw=1, linestyle='--')


    # Initialization function for FuncAnimation -------
    def init():
        scat1.set_offsets(np.empty((0, 2)))
        scat2.set_offsets(np.empty((0, 2)))
        scat3.set_offsets(np.empty((0, 2)))
        scat4.set_offsets(np.empty((0, 2)))
        scat5.set_offsets(np.empty((0, 2)))
        text.set_text("")
        return scat1, scat2, scat3, scat4, scat5, text


    # Function that is called at each iteration by FuncAnimation
    def animate(i):
        text.set_text(f'Time : {new_t[i]:.2f} s')
        x = syn_pos*lenght

        y = icl_anim[:,i]
        y2 = icl_syn_anim[:,i]
        y3 = icl_kcc2_anim[:,i]
        y4 = icl_nkcc1_anim[:,i]
        y5 = icl_leak_anim[:,i]

        scat1.set_offsets(np.c_[x, y])
        scat2.set_offsets(np.c_[x, y2])
        scat3.set_offsets(np.c_[x, y3])
        scat4.set_offsets(np.c_[x, y4])
        scat5.set_offsets(np.c_[x, y5])
        return scat1, scat2, scat3, scat4, scat5, text


    # Animation creation ------------------------------------------------------------------------------
    anim = FuncAnimation(fig, animate, init_func=init, frames=len(icl_anim[0]), interval=50, blit=True)


    # Show the animation before saving ---------------------------------------------------------------
    plt.show()
    if save_anim_icl == 1:
        writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('animation' + f'\iCl_{len(syn_pos)}synapses_kcc2={kcc2}_nkcc1={nkcc1}_-20mV.mp4', writer='ffmpeg')



if anim_open_channels:
    # Beginning and end index of the range of the simulation -------------------
    begin = int(decal/dt[0])
    end = int(time_dt1/dt[0] + time_dt2/dt[1] + time_dt3/dt[2] - end_skip/dt[2])
    puff_ind = int(time_dt1/dt[0] + time_dt2/dt[1])
    range_high_precision = [5,100]


    # Creation of the arrays for the animation ------------------------------------------------------------------------
    # The arrays are separated into three region
    # The first one has a low temporal resolution but the third one has an ever lower resolution.
    # The second region has a higher temporal resolution and is to capture the GABA puff changes correctly 
    open_states11 = open_states1[:, begin:int(puff_ind-range_high_precision[0])]
    open_states12 = open_states1[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    open_states13 = open_states1[:, int(puff_ind+range_high_precision[1]):end]
    o1_anime1 = open_states11[:, ::skip]
    o1_anime2 = open_states12
    o1_anime3 = open_states13[:, ::skip2]

    open_states21 = open_states2[:, begin:int(puff_ind-range_high_precision[0])]
    open_states22 = open_states2[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    open_states23 = open_states2[:, int(puff_ind+range_high_precision[1]):end]
    o2_anime1 = open_states21[:, ::skip]
    o2_anime2 = open_states22
    o2_anime3 = open_states23[:, ::skip2]

    open_states31 = open_states3[:, begin:int(puff_ind-range_high_precision[0])]
    open_states32 = open_states3[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    open_states33 = open_states3[:, int(puff_ind+range_high_precision[1]):end]
    o3_anime1 = open_states31[:, ::skip]
    o3_anime2 = open_states32
    o3_anime3 = open_states33[:, ::skip2]

    open_states1 = open_states[:, begin:int(puff_ind-range_high_precision[0])]
    open_states2 = open_states[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    open_states3 = open_states[:, int(puff_ind+range_high_precision[1]):end]
    o_anime1 = open_states1[:, ::skip]
    o_anime2 = open_states2
    o_anime3 = open_states3[:, ::skip2]

    t1 = t[begin:int(puff_ind-range_high_precision[0])]
    t2 = t[int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    t3 = t[int(puff_ind+range_high_precision[1]):end]
    tanime1 = t1[::skip]
    tanime2 = t2
    tanime3 = t3[::skip2]

    o1_anime = np.concatenate((o1_anime1, o1_anime2, o1_anime3), axis=1)*rnum
    o2_anime = np.concatenate((o2_anime1, o2_anime2, o2_anime3), axis=1)*rnum
    o3_anime = np.concatenate((o3_anime1, o3_anime2, o3_anime3), axis=1)*rnum
    o_anime = np.concatenate((o_anime1, o_anime2, o_anime3), axis=1)*rnum
    new_t =np.concatenate((tanime1, tanime2, tanime3))


    # Array for the label ---------------------------------------
    for_label_1 = np.full(o1_anime[0].shape, fill_value=-100000)


    # Figure object creation ----------------------------------------------------
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    scat1 = ax1.scatter([], [], lw=2, facecolor='none', color='black', s=50)
    scat2 = ax2.scatter([], [], lw=2, facecolor='none', edgecolor='black', s=50)
    scat3 = ax3.scatter([], [], lw=2, facecolor='none', edgecolor='black', s=50)
    scat4 = ax4.scatter([], [], lw=2, facecolor='black', edgecolor='black', s=50)


    # Curve drawned for the legend only --------------------------------------------
    ax1.plot(new_t, for_label_1, color='black', label=r'$Synapses$', linestyle='--')
    ax1.legend(loc="upper right")


    # Axis limits and box for the time -------------------------------------------------
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.95)
    text = ax1.text(0.02, 0.85, f'', color='black', transform=ax1.transAxes, bbox=props)
    ax1.set_ylim(np.min(o1_anime), np.max(o1_anime))
    ax2.set_ylim(np.min(o2_anime), np.max(o2_anime))
    ax3.set_ylim(np.min(o3_anime), np.max(o3_anime))
    ax4.set_ylim(np.min(o_anime), np.max(o_anime))
    ax1.set_xlabel("Position on dendrite [µm]")
    ax2.set_xlabel("Position on dendrite [µm]")
    ax3.set_xlabel("Position on dendrite [µm]")
    ax4.set_xlabel("Position on dendrite [µm]")
    ax1.set_ylabel("# of opened\nchannels (state 1)")
    ax2.set_ylabel("# of opened\nchannels (state 2)")
    ax3.set_ylabel("# of opened\nchannels (state 3)")
    ax4.set_ylabel("# of opened\nchannels")

    for j in syn_pos:
        ax1.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax2.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax3.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax4.axvline(j * lenght, color='black', lw=1, linestyle='--')


    # Initialization function for FuncAnimation
    def init():
        scat1.set_offsets(np.empty((0, 2)))
        scat2.set_offsets(np.empty((0, 2)))
        scat3.set_offsets(np.empty((0, 2)))
        scat4.set_offsets(np.empty((0, 2)))
        text.set_text("")
        return scat1, scat2, scat3, scat4, text


    # Function that is called at each iteration by FuncAnimation
    def animate(i):
        text.set_text(f'Time : {new_t[i]:.2f} s')
        x = syn_pos*lenght

        y = o1_anime[:,i]
        y2 = o2_anime[:,i]
        y3 = o3_anime[:,i]
        y4 = o_anime[:,i]

        scat1.set_offsets(np.c_[x, y])
        scat2.set_offsets(np.c_[x, y2])
        scat3.set_offsets(np.c_[x, y3])
        scat4.set_offsets(np.c_[x, y4])
        return scat1, scat2, scat3, scat4, text


    # Animation creation ------------------------------------------------------------------------------
    anim = FuncAnimation(fig, animate, init_func=init, frames=len(o1_anime[0]), interval=50, blit=True)


    # Show the animation before saving ----------------------------------------------------------------------
    plt.show()
    if save_anim_open_channels == 1:
        writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('animation' + f'\Open_states_{len(syn_pos)}synapses_kcc2={kcc2}_nkcc1={nkcc1}_-20mV.mp4', writer='ffmpeg')



if anim_all_four:
    # Beginning and end index of the range of the simulation -------------------
    begin = int(decal/dt[0])
    end = int(time_dt1/dt[0] + time_dt2/dt[1] + time_dt3/dt[2] - end_skip/dt[2])
    puff_ind = int(time_dt1/dt[0] + time_dt2/dt[1])
    range_high_precision = [5,100]


    # Creation of the arrays for the animation --------------------------------------------------------------------------------
    # The arrays are separated into three region
    # The first one has a low temporal resolution but the third one has an ever lower resolution.
    # The second region has a higher temporal resolution and is to capture the GABA puff changes correctly 
    dend_cli1 = dend_cli[:, begin:int(puff_ind-range_high_precision[0])]
    dend_cli2 = dend_cli[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_cli3 = dend_cli[:, int(puff_ind+range_high_precision[1]):end]
    dend_anime1 = dend_cli1[:, ::skip]
    dend_anime2 = dend_cli2
    dend_anime3 = dend_cli3[:, ::skip2]

    dend_gabo1 = dend_gabo[:, begin:int(puff_ind-range_high_precision[0])]
    dend_gabo2 = dend_gabo[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_gabo3 = dend_gabo[:, int(puff_ind+range_high_precision[1]):end]
    gabo_anime1 = dend_gabo1[:, ::skip]
    gabo_anime2 = dend_gabo2
    gabo_anime3 = dend_gabo3[:, ::skip2]

    dend_icl1 = dend_icl[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl2 = dend_icl[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl3 = dend_icl[:, int(puff_ind+range_high_precision[1]):end]
    icl_anime1 = dend_icl1[:, ::skip]
    icl_anime2 = dend_icl2
    icl_anime3 = dend_icl3[:, ::skip2]

    dend_icl_syn1 = dend_icl_synapses[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_syn2 = dend_icl_synapses[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_syn3 = dend_icl_synapses[:, int(puff_ind+range_high_precision[1]):end]
    icl_syn_anime1 = dend_icl_syn1[:, ::skip]
    icl_syn_anime2 = dend_icl_syn2
    icl_syn_anime3 = dend_icl_syn3[:, ::skip2]

    dend_icl_kcc21 = dend_icl_kcc2[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_kcc22 = dend_icl_kcc2[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_kcc23 = dend_icl_kcc2[:, int(puff_ind+range_high_precision[1]):end]
    icl_kcc2_anime1 = dend_icl_kcc21[:, ::skip]
    icl_kcc2_anime2 = dend_icl_kcc22
    icl_kcc2_anime3 = dend_icl_kcc23[:, ::skip2]

    dend_icl_nkcc11 = dend_icl_nkcc1[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_nkcc12 = dend_icl_nkcc1[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_nkcc13 = dend_icl_nkcc1[:, int(puff_ind+range_high_precision[1]):end]
    icl_nkcc1_anime1 = dend_icl_nkcc11[:, ::skip]
    icl_nkcc1_anime2 = dend_icl_nkcc12
    icl_nkcc1_anime3 = dend_icl_nkcc13[:, ::skip2]

    dend_icl_leak1 = dend_icl_leak[:, begin:int(puff_ind-range_high_precision[0])]
    dend_icl_leak2 = dend_icl_leak[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    dend_icl_leak3 = dend_icl_leak[:, int(puff_ind+range_high_precision[1]):end]
    icl_leak_anime1 = dend_icl_leak1[:, ::skip]
    icl_leak_anime2 = dend_icl_leak2
    icl_leak_anime3 = dend_icl_leak3[:, ::skip2]

    open_states1 = open_states[:, begin:int(puff_ind-range_high_precision[0])]
    open_states2 = open_states[:, int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    open_states3 = open_states[:, int(puff_ind+range_high_precision[1]):end]
    o_anime1 = open_states1[:, ::skip]
    o_anime2 = open_states2
    o_anime3 = open_states3[:, ::skip2]

    t1 = t[begin:int(puff_ind-range_high_precision[0])]
    t2 = t[int(puff_ind-range_high_precision[0]):int(puff_ind+range_high_precision[1])]
    t3 = t[int(puff_ind+range_high_precision[1]):end]
    tanime1 = t1[::skip]
    tanime2 = t2
    tanime3 = t3[::skip2]

    dend_anime = np.concatenate((dend_anime1, dend_anime2, dend_anime3), axis=1)
    gabo_anime = np.concatenate((gabo_anime1, gabo_anime2, gabo_anime3), axis=1)
    icl_anim = np.concatenate((icl_anime1, icl_anime2, icl_anime3), axis=1)
    icl_syn_anim = np.concatenate((icl_syn_anime1, icl_syn_anime2, icl_syn_anime3), axis=1)
    icl_kcc2_anim = np.concatenate((icl_kcc2_anime1, icl_kcc2_anime2, icl_kcc2_anime3), axis=1)
    icl_nkcc1_anim = np.concatenate((icl_nkcc1_anime1, icl_nkcc1_anime2, icl_nkcc1_anime3), axis=1)
    icl_leak_anim = np.concatenate((icl_leak_anime1, icl_leak_anime2, icl_leak_anime3), axis=1)
    o_anime = np.concatenate((o_anime1, o_anime2, o_anime3), axis=1)*rnum
    new_t =np.concatenate((tanime1, tanime2, tanime3))


    # Array for the label -------------------------------------
    for_label_1 = np.full(o_anime[0].shape, fill_value=-100000)


    # Figure object creation ---------------------------------------------------------------
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Top left panel
    ax_gaba = ax1.twinx()                              # Second ax for the GABA concentration
    line, = ax1.plot([], [], lw=2, color='orange')     # Chloride concentration
    line2, = ax_gaba.plot([], [], lw=2, color='blue') # GABA concentration

    # Top and bottom right panels
    scat1 = ax2.scatter([], [], lw=2, color='orange', s=60)
    scat2 = ax2.scatter([], [], lw=2, facecolor='none', edgecolor='darkgreen', s=50)
    scat3 = ax4.scatter([], [], lw=2, facecolor='none', edgecolor='red', s=50)
    scat4 = ax4.scatter([], [], lw=2, facecolor='none', edgecolor='deepskyblue', s=50)
    scat5 = ax4.scatter([], [], lw=2, facecolor='none', edgecolor='orchid', s=50)

    # Bottom left panel
    scat6 = ax3.scatter([], [], lw=2, facecolor='black', edgecolor='black', s=50)


    # Curve drawned for the legend only --------------------------------------------------------------
    ax1.plot(new_t, for_label_1, color='black', label=r'$Synapses$', linestyle='--')
    ax1.plot(new_t, for_label_1, color='orange', label=r'$[Cl^-]_i$')
    ax1.plot(new_t, for_label_1, color='blue', label=r'$[GABA]_o$')
    ax1.legend(loc="lower right", fontsize=14)

    ax2.scatter(new_t, for_label_1, color='orange', label=r'$I_{cl}$')
    ax2.scatter(new_t, for_label_1, facecolor='none', edgecolor='darkgreen', label=r'$I_{cl,synapse}$')
    ax4.scatter(new_t, for_label_1, facecolor='none', edgecolor='red', label=r'$I_{cl,kcc2}$')
    ax4.scatter(new_t, for_label_1, facecolor='none', edgecolor='deepskyblue', label=r'$I_{cl,nkcc1}$')
    ax4.scatter(new_t, for_label_1, facecolor='none', edgecolor='orchid', label=r'$I_{cl,leak}$')


    # Axis limits and box for the time ----------------------------------------------------------------------------------------------------------------------------------------
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.95)
    text = ax1.text(0.02, 0.85, f'', color='black', transform=ax1.transAxes, bbox=props)

    ax1.set_ylim(np.min(dend_anime), np.max(dend_anime))
    ax_gaba.set_ylim(np.min(gabo_anime), np.max(gabo_anime))

    ax2_min, ax2_max = min([np.min(icl_anim), np.min(icl_syn_anim)]), max([np.max(icl_syn_anim), np.max(icl_anim)])
    ax4_min, ax4_max = min([np.min(icl_kcc2_anim), np.min(icl_nkcc1_anim), np.min(icl_leak_anim)]), max([np.max(icl_kcc2_anim), np.max(icl_nkcc1_anim), np.max(icl_leak_anim)])
    ax2.set_ylim(ax2_min, ax2_max)
    ax4.set_ylim(ax4_min, ax4_max)
    if abs(ax2_min) > abs(ax2_max):
        ax2.legend(loc="lower right", fontsize=14)
    else:
        ax2.legend(loc="upper right", fontsize=14)
    if abs(ax4_min) > abs(ax4_max):
        ax4.legend(loc="lower right", fontsize=14)
    else:
        ax4.legend(loc="upper right", fontsize=14)

    ax3.set_ylim(np.min(o_anime), np.max(o_anime))

    ax1.set_xlim(0, record_pos[-1]*lenght)
    ax_gaba.set_xlim(0, record_pos[-1]*lenght)
    ax2.set_xlim(0, syn_pos[-1]*lenght+5)
    ax3.set_xlim(0, syn_pos[-1]*lenght+5)
    ax4.set_xlim(0, syn_pos[-1]*lenght+5)

    ax1.set_xlabel("Position on dendrite [µm]", fontsize=18)
    ax2.set_xlabel("Position on dendrite [µm]", fontsize=18)
    ax3.set_xlabel("Position on dendrite [µm]", fontsize=18)
    ax4.set_xlabel("Position on dendrite [µm]", fontsize=18)
    ax1.set_ylabel(r"$[Cl^-]_i$ [mM]", fontsize=18)
    ax_gaba.set_ylabel(r"$[GABA]_o$ [mM]", fontsize=18)
    ax2.set_ylabel("Current [pA]", fontsize=18)
    ax3.set_ylabel("# of opened\nchannels", fontsize=18)
    ax4.set_ylabel("Current [pA]", fontsize=18)

    for j in syn_pos:
        ax1.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax2.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax3.axvline(j * lenght, color='black', lw=1, linestyle='--')
        ax4.axvline(j * lenght, color='black', lw=1, linestyle='--')


    # Initialization function for FuncAnimation ---------------------------
    def init():
        line.set_data([], [])
        line2.set_data([], [])
        scat1.set_offsets(np.empty((0, 2)))
        scat2.set_offsets(np.empty((0, 2)))
        scat3.set_offsets(np.empty((0, 2)))
        scat4.set_offsets(np.empty((0, 2)))
        scat5.set_offsets(np.empty((0, 2)))
        scat6.set_offsets(np.empty((0, 2)))
        text.set_text("")
        return line, line2, scat1, scat2, scat3, scat4, scat5, scat6, text


    # Function that is called at each iteration by FuncAnimation -------
    def animate(i):
        text.set_text(f'Time : {new_t[i]:.2f} s')
        x = syn_pos*lenght
        x2 = record_pos*lenght

        y = dend_anime[:,i]
        ygaba = gabo_anime[:,i]
        y1 = icl_anim[:,i]
        y2 = icl_syn_anim[:,i]
        y3 = icl_kcc2_anim[:,i]
        y4 = icl_nkcc1_anim[:,i]
        y5 = icl_leak_anim[:,i]
        y6 = o_anime[:,i]

        line.set_data(x2, y)
        line2.set_data(x2, ygaba)

        scat1.set_offsets(np.c_[x, y1])
        scat2.set_offsets(np.c_[x, y2])
        scat3.set_offsets(np.c_[x, y3])
        scat4.set_offsets(np.c_[x, y4])
        scat5.set_offsets(np.c_[x, y5])
        scat6.set_offsets(np.c_[x, y6])
        return line, line2, scat1, scat2, scat3, scat4, scat5, scat6, text


    # Animation creation -----------------------------------------------------------------------------
    anim = FuncAnimation(fig, animate, init_func=init, frames=len(o_anime[0]), interval=50, blit=True)


    # Show the animation before saving ------------------------------------------------------------------
    plt.show()
    if save_anim_all_four == 1:
        writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('animation' + f'\All_four_{len(syn_pos)}synapses_kcc2={kcc2}_nkcc1={nkcc1}_-20mV.mp4', writer='ffmpeg')