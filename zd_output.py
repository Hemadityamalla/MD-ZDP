import numpy as np
import h5py
import matplotlib.pyplot as plt

#Read the h5 file

fi = h5py.File("streamer_2d_output.h5", "r")

#Obtain the list of databases
timeGroup_name = list(fi.keys())[0]
time_dbs_names = list(fi.get(timeGroup_name))
timeGroup = fi.get(timeGroup_name)
time_dbs = list(timeGroup.items())
#Extract the values from the database
large_array = [x[1][...] for x in time_dbs]


#Obtain the electron density for each timestep
t_steps = np.arange(len(time_dbs_names))
e_idx = 4 #Not 5 because the indexing starts from 0
e_dens = np.array([x[e_idx] for x in large_array])
plt.plot(t_steps,e_dens)
plt.yscale("log")
plt.show()

fi.close()
