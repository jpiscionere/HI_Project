#!/bin/bash

#sudo port install clang-3.8 +openmp
#port select --list clang
#sudo port select --set clang mp-clang-3.8

clang -std=c99 group_finder.c gridlink1D.c utils.c progressbar.c -I../utils/ -o group_finder  `gsl-config --cflags` `gsl-config --libs` -Wall -Wextra -g  -O3 -march=native -Ofast -DDARWIN_C_SOURCE -DISO_C_SOURCE -DPOSIX_C_SOURCE=2 -fopenmp 
