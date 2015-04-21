# agSimulations: replication code for "Perfect Competition in Markets with Adverse Selection" by Eduardo Azevedo and Daniel Gottlieb.

The scripts replicate all the figures and tables in the paper. The most important methods are implemented as general MATLAB classes. These can be adapted to other models of firm competition with adverse selection. Developed by Eduardo Azevedo and Rafael Mourao.

The scripts are divided in five folders:

## `./classes`

In this folder you can find the classes with the most important functions.

- `model.m`: abstract class that defines all properties and methods that should be implemented by a model of perfect competition with adverse selection.
	- `healthcaralognormalmodel.m`: subclass for the health insurance model in the paper with linear contracts and normally distributed losses.
	- `healthcaralognormalmodel_nl`: subclass for the health insurance model in the paper with nonlinear contracts and lognormally distributed losses.
- `population.m`: A population object describes the preferences of a finite number of consumers. This class has methods for finding competitive equilibria and optima.

## `./figuresManuscriptCompetition`

Generates figures and tables for the linear contracts health insurance model in the paper.

- `run_directory.sh`: bash script used to run the whole folder.
- `calculations_parallel.m`: creates the populations and runs the main calculations in parallel
- `welfare_table.m`: prints the `.tex` tables used in the paper
- `plot_figures.m`: plots the figures. Uses `plotFunctions/plotEquilibriumAndOptimum.m` as an auxiliary function for some figures.

## `./figuresNonLinearModel`

Generates figures and tables for the nonlinear contracts health insurance model in the paper.

- `run_directory.sh`: bash script used to run the whole folder.
- `create_population_lognormal.m`: creates the populations in parallel. Do not forget to set the desired number of workers, this program is very CPU intensive and may take a few days if running with few cores.
- `calculate_equilibrium_parallel.m`: runs the main calculations.
- `prepare_figures.m`: saves the data used for the plots.
- `plot_figures.m`: plots the figures.
- `print*` - print tex tables and values, and some `.txt` tables with additional information.

## `./tests`

Several test scripts and a folder containing the codes used for a exploratory analysis of the nonlinear problem. Very useful if you plan on implementing different models.

## `./plotFunctions`

Auxiliary functions for plotting:

- `num2bank.m` : ???
- `./export_fig` : ???

## Authors 

Eduardo Azevedo and Rafael Mourão

## Running

All code was written in MATLAB, using MATLAB's OOP implementation. Tested in an unix server, MATLAB 2014b with a parallel cluster set up. Some of the code uses KNITRO 9.1.0, but it will probably run even if knitro is not installed.

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.


