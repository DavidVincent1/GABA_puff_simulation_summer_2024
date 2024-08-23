import numpy as np
import tifffile
import os


# INPUT
# Path of the folder where the  tiff files are located
path = "Experimental_data\Tiff Intensité\C1"

# All tiff files are put into a list ---------------------------------------------------
data = []
for file in os.listdir(path):
        data.append(tifffile.imread(path + f'\{file}'))


# Conversion of the list of tiff files into a stack
stacked_data = np.stack(data)


# Saving the stack ----------------------------
tifffile.imwrite(path + '\C1.tif', stacked_data)