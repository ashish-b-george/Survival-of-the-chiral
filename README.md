# Survival-of-the-chiral
This repository contains code used for simulations in the paper  "Chirality provides a direct fitness advantage and facilitates intermixing in cellular aggregates"

Description of the files are below:

"model_simulations.c" is the code used to implement the model on the lattice used in all simulations except for the last SI figure with off-lattice simulations. Small appropriate modifications were made to code for the simulation required. This could be changing the initial conditions from circular to rectangular, well-mixed to de-mixed, etc.  Changes were also made if the parameters had to be varied in a specific way like if the chiralities were always opposite, equal or summed to a constant. Some of these changes are documented as comments at the end of the code. 

"loopScript_oppch_N.sh" is a sample bash script that loops over a set of parameters and submits jobs. It does this by modifying a child bash script in the run folder and passes the parameters to the child bash script which runs the code.

"qsubScript_oppch_N.sh" is the child bash script mentioned above that runs the code using the correct compiler and provides the appropriate linkers


