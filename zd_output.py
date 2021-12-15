import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt

pr = argparse.ArgumentParser(
    description='''Plot the output of the ZeroD plasma sim of given file''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog='''Usage: python3 zd_output.py filename'''
)

pr.add_argument("data_file", type=str, help="Filename to plot")
args = pr.parse_args()
fileName = args.data_file
#Read the h5 file
fi = h5py.File(fileName+".h5", "r")
#Read the variable names file
fi_vars = open(fileName+"_var_names.txt").read().rstrip().split("\n")

#Obtain the list of databases
timeGroup_name = list(fi.keys())[0]
time_dbs_names = list(fi.get(timeGroup_name))
timeGroup = fi.get(timeGroup_name)
time_dbs = list(timeGroup.items())
#Extract the values from the database
large_array = [x[1][...] for x in time_dbs]


# Obtaining the timesteps
tstart = 0
var_end = 0
for i,val in enumerate(fi_vars):
    if val.startswith("----"):
        var_end = i
        tstart = i+1
        break
t_step = fi_vars[tstart:]
t_steps = [float(x.lstrip()) for x in t_step]




fig, axs = plt.subplots(2,2)
#First plot has the gas properties
axs[0,0].set_title("Gas Properties")
for i in range(3):
    var_vals = np.array([x[i] for x in large_array])
    axs[0,0].plot(t_steps,var_vals, label=fi_vars[i])
axs[0,0].legend()
#Second plot has the electric field
axs[0,1].set_title("Electric field")
axs[0,1].plot(t_steps, [x[3] for x in large_array], label=fi_vars[3])
axs[0,1].legend()
#Third plot has all the species stuff
scaling = 1e-6
cmap = plt.cm.get_cmap("gist_rainbow")
axs[1,1].set_title("Specie densities")
for i in range(4,var_end):
    c = cmap(float(i)/var_end)
    var_vals = np.array([x[i] for x in large_array])
    axs[1,1].plot(t_steps, var_vals*scaling, color=c, label=fi_vars[i])
axs[1,1].set_yscale("log")
axs[1,1].set_ylim([1e0, 1e20])
axs[1,1].set_ylabel("density*"+str(scaling))
axs[1,1].legend()
plt.show()

fi.close()
