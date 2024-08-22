import numpy as np
import matplotlib.pyplot as plt
import tifffile
from scipy.optimize import curve_fit
from matplotlib import rcParams, colors


# Graph parameters ---------------
rcParams.update({'font.size': 22})


# GABA puff intensity measurement file -----------------------------------------------------------------------------
# INPUT
# Path of a tiff file. Corresponds to a stack showing the photon counts associated to the GABA during an experiment.
path = f'Experimental_data\Tiff Intensité\C1.tif'
data = tifffile.imread(path)
data_5pix = data[:, :int(data.shape[1]//5 * 5), :int(data.shape[2]//5 * 5), 0]
data_10pix = data[:, :int(data.shape[1]//10 * 10), :int(data.shape[2]//10 * 10), 0]


# Data reshaped for the different pixel sizes ---------------------------------------------------------------------------------
data_reshaped_5pix = data_5pix.reshape(data_5pix.shape[0], int(data_5pix.shape[1]/5), 5, int(data_5pix.shape[2]/5), 5)
data_reshaped_10pix = data_10pix.reshape(data_10pix.shape[0], int(data_10pix.shape[1]/10), 10, int(data_10pix.shape[2]/10), 10)

data_summed_5pix = data_reshaped_5pix.sum(axis=(2, 4))[:, :, data_5pix.shape[1]//5:data_5pix.shape[2]//5]
data_summed_10pix = data_reshaped_10pix.sum(axis=(2, 4))[:, :, data_10pix.shape[1]//10:data_10pix.shape[2]//10]

print('Shape of the data with 5x5 pixels   : ', data_summed_5pix.shape)
print('Shape of the data with 10x10 pixels : ', data_summed_10pix.shape)
print('')


# Time array of the data -----------------------------------
time_array = np.arange(0, 119.536, 0.496)
time_array_fit = np.arange(0, 88.784, 0.496)


# Image intensity ------------------------------
im_no_time = data.sum(axis=0)
im_no_time_5pix = data_summed_5pix.sum(axis=0)
im_no_time_10pix = data_summed_10pix.sum(axis=0)


# Arrays used later to stuck information --------
# Arrays to stuck wich pixel is used for the fits
im_check = np.zeros_like(im_no_time_5pix)
im_check_10pix = np.zeros_like(im_no_time_10pix)

# Arrays to stuck lifetime
tau_5pix = np.zeros_like(im_no_time_5pix, dtype=np.float64)
tau_10pix = np.zeros_like(im_no_time_10pix, dtype=np.float64)


# Fit on all counts -------------------------------------------------------------------------------------------------------------
time_no_im = data_5pix.sum(axis=1).sum(axis=1)
time_no_im_norm = (time_no_im-time_no_im.min())/(time_no_im.max()-time_no_im.min())
max_ind = np.argmax(time_no_im_norm)

# Monoexponential function with noise
def expo(t, tau, a):
    y = 1 + a*(np.exp(-(t)/tau) -1)
    return y

# Biexponential function with noise
def biexpo(t, tau, a, tau2, a2):
    y = 1 + a*(np.exp(-(t)/tau) -1) + a2*(np.exp(-(t)/tau2) -1)
    return y

# Monoexponential fit on all counts
popt_exp, pcov_exp = curve_fit(expo, time_array_fit, time_no_im_norm[max_ind:], bounds=([0, 0], [10000, 2]))
fit = expo(time_array_fit, *popt_exp)

ssr = np.sum((time_no_im_norm[max_ind:] - fit) ** 2) # SSR
sst = np.sum((time_no_im_norm[max_ind:] - np.mean(time_no_im_norm[max_ind:])) ** 2) # SST
r_squared = 1 - (ssr / sst) # R^2

print('Monoexponential Fit on all counts')
print('tau : ', popt_exp[0])
print('a   : ', popt_exp[1])
print('R^2 : ', r_squared)
print('')

# Biexponential fit on all counts
popt_exp2, pcov_exp2 = curve_fit(biexpo, time_array_fit, time_no_im_norm[max_ind:], bounds=([0, 0, 0, 0], [10000, 2, 10000, 2]))
fit2 = biexpo(time_array_fit, *popt_exp2)

ssr = np.sum((time_no_im_norm[max_ind:] - fit2) ** 2) # SSR
sst = np.sum((time_no_im_norm[max_ind:] - np.mean(time_no_im_norm[max_ind:])) ** 2) # SST
r_squared = 1 - (ssr / sst) # R^2

print('Biexponential fit on all counts')
print(r'tau  : ', popt_exp2[0])
print(r'a    : ', popt_exp2[1])
print(r'tau2 : ', popt_exp2[2])
print(r'a2   : ', popt_exp2[3])
print(r'R^2  : ', r_squared)
print('')


# Fit per pixel 5x5 ----------------------------------------------------------------------------------------------------------------------------------------
# Arrays used to stuck results
curve_near_far_from_puff = []
curve_near_far_from_puff_time = []
data_near_far_from_puff = []
data_near_far_from_puff_time = []
data_pix = []
fit_pix = []
param_pix = []
time_pix = []
time_for_graph = []
R_2 = []
data_norm_pix = 0
for i in range(data_summed_5pix.shape[1]):
    for j in range(data_summed_5pix.shape[2]):
        # Index of the maximum in the pixel
        max_ind_ = data_summed_5pix[:, i, j].argmax()

        if 58 < max_ind_ < 68:                                # Index range of the GABA puff                              
            if 1500 < data_summed_5pix[:, i, j].sum() < 4000: # Exclude too low and too high number of photons

                val = data_summed_5pix[max_ind_:, i, j] # Counts before maximum ignored
                im_check[i, j] = 1                      # Confirms that a fit was done at that pixel

                # Normalization of the counts of the pixel
                data_norm_pix = (data_summed_5pix[:, i, j]-data_summed_5pix[:, i, j].min())/(data_summed_5pix[:, i, j].max()-data_summed_5pix[:, i, j].min())

                # Time array with zero shifted to maximum value
                time = np.linspace(0, 119.536-((max_ind_+1)*0.496), len(data_norm_pix[max_ind_:]))

                # Monoexponential fit
                popt_exp, pcov_exp = curve_fit(expo, time, data_norm_pix[max_ind_:], bounds=([0, 0], [10000, 2]))
                tau_5pix[i, j] = popt_exp[0]

                # Stucking the information used later
                data_pix.append(data_norm_pix)
                time_for_graph.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_5pix[:, i, j])))
                fitt = expo(time, *popt_exp)
                fit_pix.append(fitt)
                param_pix.append((popt_exp[0], popt_exp[1]))
                time_pix.append(time)
                R_2.append(1 - (np.sum((data_norm_pix[max_ind_:] - fitt) ** 2) / np.sum((data_norm_pix[max_ind_:] - np.mean(data_norm_pix[max_ind_:])) ** 2)))

                # Curves near puff and far from puff
                if i == 22 and j == 38:
                    curve_near_far_from_puff.append(fitt)
                    curve_near_far_from_puff_time.append(time)
                    data_near_far_from_puff.append(data_norm_pix)
                    data_near_far_from_puff_time.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_5pix[:, i, j])))
                elif i == 10 and j == 51:
                    curve_near_far_from_puff.append(fitt)
                    curve_near_far_from_puff_time.append(time)
                    data_near_far_from_puff.append(data_norm_pix)
                    data_near_far_from_puff_time.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_5pix[:, i, j])))


# Fit per pixel 10x10 ----------------------------------------------------------------------------------------------------------------------------------------
# Arrays used to stuck results
data_pix10 = []
fit_pix10 = []
param_pix10 = []
time_pix10 = []
time_for_graph10 = []
R_2_10 = [] # 30, 24
data_norm_pix10 = 0
for i in range(data_summed_10pix.shape[1]):
    for j in range(data_summed_10pix.shape[2]):
        # Index of the maximum in the pixel
        max_ind_ = data_summed_10pix[:, i, j].argmax()

        if 58 < max_ind_ < 68:                                     # Index range of the GABA puff
            if 1500*4 < data_summed_10pix[:, i, j].sum() < 4000*4: # Exclude too low and too high number of photons (4 times the threshold of 5x5 pixels)

                val = data_summed_10pix[max_ind_:, i, j] # Counts before maximum ignored
                im_check_10pix[i, j] = 1                 # Confirms that a fit was done at that pixel

                # Normalization of the counts of the pixel
                data_norm_pix10 = (data_summed_10pix[:, i, j]-data_summed_10pix[:, i, j].min())/(data_summed_10pix[:, i, j].max()-data_summed_10pix[:, i, j].min())

                # Time array with zero shifted to maximum value
                time = np.linspace(0, 119.536-((max_ind_+1)*0.496), len(data_norm_pix10[max_ind_:]))

                # Monoexponential fit
                popt_exp, pcov_exp = curve_fit(expo, time, data_norm_pix10[max_ind_:], bounds=([0, 0], [10000, 2]))
                tau_10pix[i, j] = popt_exp[0]

                # Stucking the information used later
                data_pix10.append(data_norm_pix10)
                time_for_graph10.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_10pix[:, i, j])))
                fitt = expo(time, *popt_exp)
                fit_pix10.append(fitt)
                param_pix10.append((popt_exp[0], popt_exp[1]))
                time_pix10.append(time)
                R_2_10.append(1 - (np.sum((data_norm_pix10[max_ind_:] - fitt) ** 2) / np.sum((data_norm_pix10[max_ind_:] - np.mean(data_norm_pix10[max_ind_:])) ** 2)))

                # Curves near puff and far from puff
                if i == 10 and j == 20:
                    curve_near_far_from_puff.append(fitt)
                    curve_near_far_from_puff_time.append(time)
                    data_near_far_from_puff.append(data_norm_pix10)
                    data_near_far_from_puff_time.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_10pix[:, i, j])))
                elif i == 3 and j == 33:
                    curve_near_far_from_puff.append(fitt)
                    curve_near_far_from_puff_time.append(time)
                    data_near_far_from_puff.append(data_norm_pix10)
                    data_near_far_from_puff_time.append(np.linspace(-((max_ind_+1)*0.496), 119.536-((max_ind_+1)*0.496), len(data_summed_10pix[:, i, j])))


# Mean lifetime for the 5x5 pixels ----------------------------------------------------------------------------------
# List of tau
for_mean_5pix = np.asarray([param_pix[i][0] for i in range(len(param_pix))])

# List of R^2
R_2_for_mean_5pix = np.asarray(R_2) 
R_2_mean_5pix = R_2_for_mean_5pix.mean() # Mean R^2

# Mean tau and R^2 keeping the best fits
mean_with_high_R_2_5pix = []
mean_tau_with_high_R_2_5pix = []
for r, val in enumerate(R_2):
    if val > R_2_mean_5pix:
        mean_with_high_R_2_5pix.append(val)
        mean_tau_with_high_R_2_5pix.append(for_mean_5pix[r])

mean_with_high_R_2_5pix = np.asarray(mean_with_high_R_2_5pix)
mean_tau_with_high_R_2_5pix = np.asarray(mean_tau_with_high_R_2_5pix)

print("Mean tau/R^2 for 5 x 5 pixels (mean ± std)")
print('Mean tau                    : ', for_mean_5pix.mean(), "±", for_mean_5pix.std())
print('Mean R2                     : ', R_2_mean_5pix, "±", R_2_for_mean_5pix.std())
print('Mean R2 keeping best fits   : ', mean_with_high_R_2_5pix.mean(), "±", mean_with_high_R_2_5pix.std())
print('Mean tau keeping best fits  : ', mean_tau_with_high_R_2_5pix.mean(), "±", mean_tau_with_high_R_2_5pix.std())
print('')


# Mean lifetime for the 10x10 pixels ---------------------------------------------------------------------------
# List of tau
for_mean_10pix = np.asarray([param_pix10[i][0] for i in range(len(param_pix10))]) 

# List of R^2
R_2_for_mean_10pix = np.asarray(R_2_10)
R_2_mean_10pix = R_2_for_mean_10pix.mean() # Mean R^2

# Mean tau and R^2 keeping the best fits
mean_with_high_R_2_10pix = []
mean_tau_with_high_R_2_10pix = []
for r, val in enumerate(R_2_10):
    if val > R_2_mean_10pix:
        mean_with_high_R_2_10pix.append(val)
        mean_tau_with_high_R_2_10pix.append(for_mean_10pix[r])

mean_with_high_R_2_10pix = np.asarray(mean_with_high_R_2_10pix)
mean_tau_with_high_R_2_10pix = np.asarray(mean_tau_with_high_R_2_10pix)

print("Mean tau/R^2 for 10 x 10 pixels (mean ± std)")
print('Mean tau                   : ', for_mean_10pix.mean(), for_mean_10pix.std())
print('Mean R2                    : ', R_2_mean_10pix, R_2_for_mean_10pix.std())
print('Mean R2 keeping best fits  : ', mean_with_high_R_2_10pix.mean(), mean_with_high_R_2_10pix.std())
print('Mean tau keeping best fits : ', mean_tau_with_high_R_2_10pix.mean(), mean_tau_with_high_R_2_10pix.std())
print('')


# Graphs ---------------------------------------------------------------------------------------------------------
# Intensity graph (initial data)
plt.figure("Image")
plt.imshow(im_no_time)
cax = plt.axes((0.1, 0.25, 0.8, 0.035))
plt.colorbar(cax=cax, orientation='horizontal', label='Photon counts [-]')


# Intensity graph (with 5x5 and 10x10 pixels)
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title('Pixel_reduction')
ax[0].set_title("5 x 5 pixels")
im1 = ax[0].imshow(im_no_time_5pix)
ax[1].set_title("10 x 10 pixels")
im2 = ax[1].imshow(im_no_time_10pix)
fig.colorbar(im1, ax=ax[0], orientation='vertical', fraction=0.046, pad=0.04, label='Photon counts')
fig.colorbar(im2, ax=ax[1], orientation='vertical', fraction=0.046, pad=0.04, label='Photon counts')


# Show wich pixel was used to fit
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title('Fited_pixels')
ax[0].set_title("5 x 5 pixels")
ax[0].imshow(im_check)
ax[1].set_title("10 x 10 pixels")
ax[1].imshow(im_check_10pix)


# Lifetime per pixel
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title('Fited_pixels_tau')
ax[0].set_title("5 x 5 pixels")
ax[1].set_title("10 x 10 pixels")
images = []
data_ = [tau_5pix, tau_10pix]

# This code is not from me
for i in range(2):
    # Generate data with a range that varies from one plot to the next.
    images.append(ax[i].imshow(data_[i]))

# Find the min and max of all colors for use in setting the color scale.
vmin = min(image.get_array().min() for image in images)
vmax = max(image.get_array().max() for image in images)
norm = colors.Normalize(vmin=vmin, vmax=vmax)
for im in images:
    im.set_norm(norm)
fig.colorbar(images[0], ax=ax, orientation='vertical', fraction=.1, label=r"$\tau [s]$")

def update(changed_image):
    for im in images:
        if (changed_image.get_cmap() != im.get_cmap()
                or changed_image.get_clim() != im.get_clim()):
            im.set_cmap(changed_image.get_cmap())
            im.set_clim(changed_image.get_clim())
for im in images:
    im.callbacks.connect('changed', update)


# Curve fits on all photon counts
plt.figure('Curve_fit_all_counts')
plt.plot(time_array, time_no_im_norm, color='blue')
plt.plot(time_array[max_ind:], fit, color='red', label='Monoexponential')
plt.plot(time_array[max_ind:], fit2, color='orange', label='Biexponential')
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Normalized intensity [-]")


# Show all the curve fits
fig, ax = plt.subplots(2, 1)
fig.canvas.manager.set_window_title('Fited_pixels_all_curves')
ax[0].set_title("Curve fit 5x5 pixels")
for i,j in zip(time_pix, fit_pix):
    ax[0].plot(i, j)
ax[1].set_title("Curve fit 10x10 pixels")
for i,j in zip(time_pix10, fit_pix10):
    ax[1].plot(i, j)
ax[0].set_ylabel("Normalized\nintensity [-]")
ax[1].set_ylabel("Normalized\nintensity [-]")
ax[1].set_xlabel("Time [s]")


# Pixel 5x5 : (24,30) and (10,58) -> ax[0,0] and ax[0,1]
# Pixel 10x10 : (10,20) and (3,33) -> ax[1,0] and ax[1,1]
# This graph showes that pixel near the puff corresponds to better fit
fig, ax = plt.subplots(2,2)
fig.canvas.manager.set_window_title('Fit_near_and_far_from_puff')
ax[0,0].set_title("5x5 pixel, near puff")
ax[0,1].set_title("5x5 pixel, far from puff")
ax[1,0].set_title("10x10 pixel, near puff")
ax[1,1].set_title("10x10 pixel, far from puff")
ax[0,0].plot(curve_near_far_from_puff_time[1], curve_near_far_from_puff[1])
ax[0,0].plot(data_near_far_from_puff_time[1], data_near_far_from_puff[1])
ax[0,1].plot(curve_near_far_from_puff_time[0], curve_near_far_from_puff[0])
ax[0,1].plot(data_near_far_from_puff_time[0], data_near_far_from_puff[0])
ax[1,0].plot(curve_near_far_from_puff_time[3], curve_near_far_from_puff[3])
ax[1,0].plot(data_near_far_from_puff_time[3], data_near_far_from_puff[3])
ax[1,1].plot(curve_near_far_from_puff_time[2], curve_near_far_from_puff[2])
ax[1,1].plot(data_near_far_from_puff_time[2], data_near_far_from_puff[2])

ax[1,0].set_xlabel("Time [s]")
ax[1,1].set_xlabel("Time [s]")

ax[0,0].set_ylabel("Normalized\ncounts")
ax[1,0].set_ylabel("Normalized\ncounts")


plt.show()