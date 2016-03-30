#!/bin/bash
#$ -cwd
#$ -m ae
#$ -M eduaze@gmail.com
#$ -q bigram
#$ -j y

#This script runs the simulations for the non-linear model, plotting all figures and tables
matlab -nodisplay < create_population_lognormal.m
matlab -nodisplay < calculate_equilibrium_parallel.m
matlab -nodisplay < print_values.m
matlab -nodisplay < print_tex_tables.m
matlab -nodisplay < prepare_figures.m
matlab < plot_figures.m
