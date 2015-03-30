#!/bin/bash
qsub -N population -m ae -M rafael@mourao.org -M eazevedo@wharton.upenn.edu -q bigram -b y -j y "matlab -nodisplay < create_population_lognormal.m"
qsub -N calculation -m ae -M rafael@mourao.org -M eazevedo@wharton.upenn.edu -hold_jid population -q bigram -b y -j y "matlab < calculate_equilibrium_parallel.m"