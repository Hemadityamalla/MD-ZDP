# MD-ZDP

- A simple code to perform 0D chemistry simulations occuring in electrical discharges.
- Some of the modules used in this code were taken from afivo-streamer (give the link)
- TODO: Improve the initial condition inputs. Right now, we can just specify the electron density and the first positive ion density.

#### How to run the code

- `make`
- `./zdPlasmaChem streamer_2d.cfg`

#### Immediate TODOs
- Solve the M and X/Ar neutral gas variable duplications.
- Functionality to add initial specie density-- e, pos, negative, (anything else?)
- Check and remove the transport data module
- Make the constant/variable gas density consistent throughtout the code.
- Add example with a repetitive pulse 

#### Not-so immediate TODOs
- Add options to write out source/sink terms for the species
- Output rate constants, etc
- ADD UNIT AND REGRESSION TESTS
- Think of just generating txt files because the hd5 files being generated are too large!


### Future work
- Write a CMakeLists.txt for this project.
- Add the functionality to use dvode modules for the time integration.
- Option to output the rate constants as a matrix
- Option to output the rates/rate constants as a function of some parameters.
- FULLY DOCUMENT HOW THE HDF5 STUFF WAS INSTALLED AND THE PROBLEMS
### Installing HDF5 for writing the output
- I used HDF5 to write the data. For that, I used the following tutorial: https://www.youtube.com/watch?v=ysqvDunaXNA . It also has instructions on how to install the HDF5 libraries in your system
