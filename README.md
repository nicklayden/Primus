# Primus

Working version of the Barnes-Hut NBody simulator in the Cpp repository. 

Primus now uses SFML for real-time rendering of the simulation. It is built to be run from the command line, along with a number of optional flags.

A CUDA accelerated version of the program will be added eventually. Currently, cuda is a bit of a nightmare to work with across systems.

## Example gif
Just a small gif showing the underlying Barnes hut tree, along with the Particles in the simulation. This one includes:
Mars, Earth, Sun, and Jupiter, as well as 10 "asteroid" sized particles.
![alt text](https://github.com/nicklayden/Primus/blob/master/Peek%202017-05-15%2023-55.gif "Nbody Simulation")


## Installation
Clone this repository
```
 git clone http:://github.com/nicklayden/Primus
 cd Primus
```
Invoke the makefile
```
make
```
## Running
To run the program with default settings
```
./primus
```

