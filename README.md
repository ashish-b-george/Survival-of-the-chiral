# Survival-of-the-chiral
This repository contains code used for simulations in the paper  "Chirality provides a direct fitness advantage and facilitates intermixing in cellular aggregates"

Description of the files are below:

"Figures_python_notebook.ipynb" is a python notebook used to plot all the figures in the paper.

-----------------Code for simulations on lattice----------------
"model_simulations.c" is the code used to implement the model on the lattice used in all simulations except for the last SI figure with off-lattice simulations. Small appropriate modifications were made to code for the simulation required. This could be changing the initial conditions from circular to rectangular, well-mixed to de-mixed, etc.  Changes were also made if the parameters had to be varied in a specific way like if the chiralities were always opposite, equal or summed to a constant. Some of these changes are documented as comments at the end of the code. 

-----------------Bash scripts to submit the on lattice code----------------
"loopScript_oppch_N.sh" is a sample bash script that loops over a set of parameters and submits jobs. It does this by modifying a child bash script in the run folder and passes the parameters to the child bash script which runs the code.

"qsubScript_oppch_N.sh" is the child bash script mentioned above that runs the code using the correct compiler and provides the appropriate linkers


-----------------Off Lattice Simulations----------------
"off_lattice_population.hpp", "off_lattice_main.cpp", "off_lattice_oppch63_f0.sh" are used for off lattice simulations.
"off_lattice_population_class.py", "submit_python_analysis.sh", "off_lattice_avg_velocity_cluster", "off_lattice_cluster_plot.py", "off_lattice_parameter_analysis"  are python scripts to analyze and plot the off lattice simulations. 

"off_lattice_population.hpp" is the header file for the C++ class that runs the simulations. It has a few different ways of implementing chiral migration in range expansions. We only used one for the paper. The results in the paper are independent of the different ways.
"off_lattice_main.cpp" is the main function that runs the code as required.
"off_lattice_oppch63_f0.sh" is the bash script for submitting parameters
"off_lattice_population_class.py" is a python class that is used by the other anlysis scripts to analyze the data from off-lattice simulations.


-----------------various script for analyzing on lattice simulation results ----------------
"new_extract_all_walls_line_cluster_nojumps_destpath.py", "new_extract_all_walls_line_cluster_nojumps.py"	and "new_extract_wall_circle.py" are used to extract the domain boundaries in rectangular and circular expansions respectively

"new_fstar_fbar_fig.py" is used for the f* vs fbar plot

"new_neutralHdecay_htime_avg_cluster.py", "plot_htime.py", "new_neutralHdecay_htime_plot.py" plots Heterozygosity in time from simulation in Fig 8 adn SI figs. The first file averages heterozygosity over runs with same parameters but different realizations of the noise, the second and third make plots.

"new_plot_1concentration.py" plots the image of the colony from simulation output.

"wall_extraction_bulge_analysis.py", "bulge_slope_analysis_plots.py", "wall_extraction_bulge_analysisPlots.py", "bulge_slope_analysis.py" were directly of indirectly involved in  Fig S2.
