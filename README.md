# JAMOVET

These files are for the manuscript "A novel approach for estimating vaccine efficacy for infections with multiple disease outcomes: application to a COVID-19 vaccine trial".

The R scripts for the analysis of the University of Oxford/AstraZeneca sponsored trial, COV002 are:
 - "1. COV002 VE analysis GitHub", which runs the models and generates the results for the analysis of COV002, with and without bias and covariate adjustment, and
 - "2. "COV002 paper figures GitHub", which generates the tables and figures presented in the manuscript.

The R scripts for the simulation analysis are "Simulation_setup_scenario1/2/3/4" and each file beginning with "modscript". The simulation setup files create the 4 datasets for each of the four scenarios described in the paper. All separate analyses are run in "modscript_separate", while each of the joint analyses are run in the respective joint modscript files.
