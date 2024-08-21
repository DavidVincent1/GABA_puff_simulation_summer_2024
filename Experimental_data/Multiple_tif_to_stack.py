import numpy as np
import tifffile
import os


# Path of the folder where the  tiff files are located
path = "Données_expérimentales\Tiff Intensité\C1"

# All tiff files are put into a list ---------------------------------------------------
data = []
for file in os.listdir(path):
        data.append(tifffile.imread(f'Données_expérimentales\Tiff Intensité\C1/{file}'))


# Conversion of the list of tiff files into a stack
stacked_data = np.stack(data)