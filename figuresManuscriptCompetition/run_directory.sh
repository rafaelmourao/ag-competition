#!/bin/bash
#$ -cwd
#$ -m ae
#$ -M eduaze@gmail.com
#$ -q bigram
#$ -j y

export MATLAB_PARALLEL_RAM=31G
matlab -nodisplay < calculations_parallel.m;
matlab -nodisplay < welfare_table.m;
matlab -nodisplay < comparative_statics.m
matlab < plot_figures.m;

ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start; sleep 120;'
ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start'
