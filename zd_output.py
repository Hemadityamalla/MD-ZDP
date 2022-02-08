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
#Read the file containing variable names, reaction names, and time steps
fi_vars = open(fileName+"_var_names.txt").read().rstrip().split(23*"-")

#Obtain the list of databases
timeGroup_name = list(fi.keys())[0]
rateGroup_name, specieGroup_name = list(fi.get(timeGroup_name))
time_dbs_names = list(fi.get(timeGroup_name+"/"+rateGroup_name))
#Extract the specie densities data
specieGroup = fi.get(timeGroup_name+"/"+specieGroup_name)
time_dbs = list(specieGroup.items())
specie_array = [x[1][...] for x in time_dbs]
#Extracting the rates data
rateGroup = fi.get(timeGroup_name+"/"+rateGroup_name)
time_dbs = list(rateGroup.items())
rate_array = [x[1][...] for x in time_dbs]


#Obtaining the variable names
var_names = fi_vars[0].split("\n")
var_end = len(var_names)-1
#Obtaining and cleaning the reaction names
reac_names = [x.strip() for x in fi_vars[1].split("\n")]
reac_names.remove("")
reac_names.remove("")
#Obtaining and cleaning the timesteps
t_steps = [x.strip() for x in fi_vars[2].split("\n")]
t_steps.remove("")
t_steps = list(map(float, t_steps))





fig, axs = plt.subplots(2,1)
#First plot has the gas properties
axs[0].set_title("Gas Properties")
for i in range(3):
    var_vals = np.array([x[i] for x in specie_array])
    axs[0].plot(t_steps,var_vals, label=var_names[i])
axs[0].legend()
#Second plot has the electric field
axs[1].set_title("Electric field")
axs[1].plot(t_steps, [x[3] for x in specie_array], label=var_names[3])
axs[1].legend()
#Third plots have all the species stuff
scaling = 1e-6
cmap = plt.cm.get_cmap("tab20")
#Obtain the number of chemical species to generate plots
n_spec = var_end - 4
spec_per_plot = 5
n_spec_plots = int(n_spec/spec_per_plot) + 1
#figure,axes=plt.subplots(n_spec_plots,1)
#axes[0].set_title("Specie densities")
#i_plot = 0
for i_specie in range(4,var_end):
    #c = cmap(float(i_specie)/var_end)
    var_vals = np.array([x[i_specie] for x in specie_array])
    i_plot = int((i_specie-4)/spec_per_plot)
    plt.figure(i_plot+2)
    #plt.plot(t_steps, var_vals*scaling, color=c,  
    #    label=fi_vars[i_specie])
    plt.plot(t_steps, var_vals*scaling,label=var_names[i_specie])

    plt.legend()
    



#Plotting the reaction rates
reac_per_plot = 6
n_reac_plots = int(len(reac_names)/reac_per_plot) + 1
for i_reac in range(0, len(reac_names)):
    var_vals = np.array([x[i_reac] for x in rate_array])
    i_reac_plot = int(i_reac/reac_per_plot)
    plt.figure(i_reac_plot+n_spec_plots)
    plt.plot(t_steps, var_vals*scaling, label=reac_names[i_reac])
    plt.legend()

plt.legend()
plt.show()
fi.close()
