# agSimulations: replication code for Azevedo and Gottlieb's paper "Perfect Competition in Markets with Adverse Selection"

Here you find the replication code for Azevedo and Gottlieb's paper "Perfect Competition in Markets with Adverse Selection". All code was written in MATLAB, using MATLAB's OOP implementation. While only the source code is included here, the populations used to generate the figures in the paper can be made available on request.

In the present version of the repository, the scripts are divided in five folders:

## ./classes

In this folder you can find the main classes used to run the examples. More documentation can be found inside the files.

- model.m: abstract class that defines all properties and methods that should be implemented by a model.
	- healthcaralognormalmodel.m: subclass that defines the properties and methods of the linear case.
	- healthcaralognormalmodel_nl: subclass that defines the properties and methods of the nonlinear case.
- population.m: defines an object that holds a population (generated according to a model class) and implements the algorithms to find optimal and efficient prices and allocations.

## ./figuresManuscriptCompetition

Codes used to generate the linear example figures included in the manuscript.

- run_directory.sh: bash script used to run the whole folder.
- calculations_parallel.m: creates the populations and run the main calculations in parallel
- welfare_table.m: prints the tex tables used in the paper
- plot_figures.m: plots the figures. Uses plotFunctions/plotEquilibriumAndOptimum.m as an auxiliary function for some figures.

## ./figuresNonLinearModel

- run_directory.sh: bash script used to run the whole folder.
- create_population_lognormal.m: creates the populations in parallel. Do not forget to set the desired number of workers, this program is very CPU intensive and may take a few days if running with few cores.
- calculate_equilibrium_parallel.m: runs the main calculations.
- prepare_figures.m: saves the data used for the plots.
- plot_figures.m: plots the figures.
- print* - print tex tables and values, and some .txt tables with additional information.

## ./tests

Several test scripts and a folder containing the codes used for a exploratory analysis of the nonlinear problem. Very useful if you plan on implementing different models.

## ./plotFunctions

Auxiliary functions for plotting:

- num2bank.m : ???
- ./export_fig : ???

## Authors 

Eduardo Azevedo and Rafael Mourão

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.


