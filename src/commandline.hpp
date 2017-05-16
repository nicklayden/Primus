#ifndef COMMANDLINE_HPP
#define COMMANDLINE_HPP
/*
  Command line options for Barnes Hut simulation program.
*/

/* Standard Libraries */
#include <iostream>
#include <vector>

/* Boost Libraries */
#include "boost/program_options.hpp"

int CommandLineSettings(int argc, char** argv,double* timestep,\
                        uint* Nparticles ,uint* Niterations ,float* theta,\
                        float* softener, bool* stats, uint* freq, uint* fpsmax, \
                        bool* drawNodes);

#endif /* end of include guard: COMMANDLINE_HPP */
