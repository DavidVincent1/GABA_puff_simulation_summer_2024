import pickle
import matplotlib.pyplot as plt
import os


# INPUT
# This script opens the pickle files and display the graphs
#filepath = r"Path of the folder where the pickle files are"
filepath = r"dataset\pickle_kcc2=1e-06-0.0005_nkcc1=1e-06-0.0005_gclc2=0.0_unclamp"
graphs = []
for file in os.listdir(filepath):
    print(file)
    graphs.append(pickle.load(open(filepath + r'''\\''' + file, 'rb')))

graphs[0].show()

plt.show()