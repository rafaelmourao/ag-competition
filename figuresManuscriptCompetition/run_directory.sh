#!/bin/bash

matlab -nodisplay < calculations_parallel.m;
matlab -nodisplay < welfare_table.m;
matlab < plot_figures.m;

ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start; sleep 120;'
ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start'
