# MD-ZDP

- A simple code to perform 0D chemistry simulations occuring in electrical discharges.
- Some of the modules used in this code were taken from afivo-streamer (give the link)
- TODO: Improve the initial condition inputs. Right now, we can just specify the electron density and the first positive ion density.

#### How to run the code

- `make`
- `./zdPlasmaChem streamer_2d.cfg`


### Future work
- Write a CMakeLists.txt for this project.
- Add the functionality to use dvode modules for the time integration.
- Option to output the rate constants as a matrix
- Option to output the rates/rate constants as a function of some parameters.
