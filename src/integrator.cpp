/*

  verlet integrator cpp file.

*/
#include "integrator.hpp"

VerletIntegrator::VerletIntegrator(Particle* particle)
{
    std::cout << "Creating Integrator" << std::endl;
    currentState.x = particle->rx;
    currentState.y = particle->ry;
    currentState.vx = particle->vx;
    currentState.vy = particle->vy;
    currentState.ax = particle->ax;
    currentState.ay = particle->ay;

}
