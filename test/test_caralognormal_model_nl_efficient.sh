#!/bin/bash

matlab < test_caralognormal_model_nl_efficient.m

ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start; sleep 120;'
ssh unix.wharton.upenn.edu 'dropbox stop; sleep 20; dropbox start'
