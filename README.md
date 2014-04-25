## onezone

A simple one-zone model of Type I X-ray bursts in python. The accretion rate in Eddington units and the time to run in seconds should be given on the command line. 

Note that the code uses simplified microphysics, e.g. it assumes ideal gas for the heat capacity, sets the free-free Gaunt factor to 1, and includes only the triple alpha reaction. Nonetheless it nicely illustrates the expected range of bursting behavior: flashes at low accretion rates, stable burning at high accretion rates, and oscillatory burning near marginal stability. This code is similar to the one-zone models shown in [Heger et al. (2007)](http://arxiv.org/abs/astro-ph/0511292).

e.g. `python onezone.py 1.0 1e4` shows a series of flashes

![lightcurve](https://github.com/andrewcumming/onezone/raw/master/mdot=1.png)
![column depth-temperature](https://github.com/andrewcumming/onezone/raw/master/mdot=1_yT.png)

e.g. `python onezone.py 6.0 1e4` shows the evolution to stable burning

![column depth-temperature](https://github.com/andrewcumming/onezone/raw/master/mdot=6_yT.png)

