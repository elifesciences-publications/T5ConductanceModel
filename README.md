# T5ConductanceModel

This repository contains the Matlab code to generate the model data used in "The computation of directional selectivity in the Drosophila OFF motion pathway", by Gruntman, Romani, and Reiser (2019).

In order to execute this code, the Matlab Optimization Toolbox (https://www.mathworks.com/products/optimization.html) needs to be installed, and the Matlab working directory needs to be set to the folder containing the files in this repository.

Running the model optimization:
Within Matlab, from the folder containing all the code and .mat files, run the command:
optimize_model(rng_seed), where rng_seed is a string containing an integer used as a seed for the random number generator (e.g. optimize_model("42")). This seed is used to randomly select an initial starting point for the optimization in the space of model parameters.

The optimization code generates a subfolder named rng_seed (e.g. 42) with seventeen .mat files containing the optimization results for the seventeen experimentally recorded T5 neurons. For instance, the results for neuron #1 can be found in result_cell_1.mat

Each file contains the simulated membrane potential traces in response to the simulated presentation of (i) single bars (ii) minimal motion stimuli (iii) moving bars (iv) static gratings (v) drifting gratings. 
