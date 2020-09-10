# Baryogenesis through Asymmetric Hawking Radiation from Primordial Black Holes as Dark Matter

Code repository for the PBH and baryon density fractions and yield from asymmetric Hawking radiation using the broad mass spectrum of 2001.04371. 

## Overview and running 

This repository contains the c++ codes used to reproduce the results of Section III A and B of 1404.0113 (bh_hook.cpp) and the results of <> (bh_clean.cpp). A python notebook is also included which was used to generate the plots of <>.  

The codes bh_clean.cpp and bh_hook.cpp can be run from the command line using a similar command to: 

>g++ -lgsl  -I/Users/bbose/Desktop/Black_hole/git  bh_clean.cpp 

replacing the include directory to the header files directory. Then just run the executable produced with

>./a.out 

The code is fairly well commented but if you run into any issues, contact Ben Bose at nebblub@gmail.com. 

## Requirements
### C++ Compiler and automake
The codes are written in C++, so you'll need a relatively modern C++ compiler.

### GSL
The code also makes use of [gsl](http://www.gnu.org/software/gsl/). You will need to have a version of gsl installed (gsl 2. or later) needed for integrations of the mass spectrum. 


