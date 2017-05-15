/*

  First Integrator class for use in the Barnes Hut simulator program.
  Stick with Verlet Integration for now.

*/
#include "Particle.hpp"


struct state
{
    double x, y, vx, vy, ax, ay;
};

/*
  VerletIntegrator class, takes a particle object on construction, uses the
  properties of the particle to set up an integration state for the particle.

  In this case, particle is loaded, the state is updated (force calculated again)
  then the verlet method is implemented.
*/
class VerletIntegrator
{
  public:
    //constructors
    VerletIntegrator(Particle* particle);

    //container that holds the current state (current (x,y), (vx,vy), (ax,ay))
    state currentState;
    state nextState;

    //functions
    void stepOne();
    void stepTwo();


};
