import pickle
import matplotlib.pyplot as plt
import os


# This script opens the pickle files and display the graphs
#filepath = r"Path of the folder where the pickle files are"
filepath = r"Path of the folder where the pickle files are"
graphs = []
for file in os.listdir(filepath):
    print(file)
    graphs.append(pickle.load(open(filepath + r'''\\''' + file, 'rb')))

graphs[0].show()

plt.show()