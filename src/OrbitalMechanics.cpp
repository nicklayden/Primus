/*
  Orbital Mechanics equations.
    -Calculating aphelion and perihelion velocities,
    -Energy calculations?
    -Other shit.

*/


#include "OrbitalMechanics.hpp"
#include <math.h>
#include <vector>

double perihelionVelocity(double G, double msun, double rperi, double e)
{
    return sqrt(G*msun*(1.0-e)/rperi);
}

double aphelionVelocity(double G, double msun, double rap, double e)
{
    return sqrt(G*msun*(1.0-e)/rap);
}

std::vector<double> tiltedAphelion(double vap, double angle)
{
  // Calculates the vx and vy of a particle at a given d,e and angle.
  std::vector<double> velocity;
  double vx;
  double vy;
  vx = vap*sin(angle);
  vy = vap*cos(angle);
  velocity.push_back(vx);
  velocity.push_back(vy);
  return velocity;
}

std::vector<double> tiltedPerihelion(double vperi, double angle)
{
  // Calculates the vx and vy of a particle at a given d,e and angle.
  std::vector<double> velocity;
  double vx;
  double vy;
  vx = vperi*sin(angle);
  vy = vperi*cos(angle);
  velocity.push_back(vx);
  velocity.push_back(vy);
  return velocity;
}

double CircularVelocity(double G, double msun, double r)
{
  return sqrt(G*msun/r);
}
