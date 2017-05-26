#ifndef ORBITALMECHANICS_HPP
#define ORBITALMECHANICS_HPP
#include <vector>
#include "Particle.hpp"

typedef std::vector<Particle*> Body_ctr;

// Orbital parameter functions
double perihelionVelocity(double G, double msun, double rperi, double e);
double aphelionVelocity(double G, double msun, double rap, double e);
std::vector<double> tiltedAphelion(double vap, double angle);
std::vector<double> tiltedPerihelion(double vperi, double angle);
double CircularVelocity(double G, double msun, double r);

// Solar system initialization
void SpawnSolarSystemPlanets(Body_ctr& bodies);

// Particle Disc Creators
void initParticleDisk(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles);
void makeCircularDisc(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles);

// Conserving Solar Momentum before simulation starts
void FixSunMomentum(Body_ctr bodies);

// Random number generators
double unirandomval(double min, double max);
double normrandomval(double mean, double variance);


#endif /* end of include guard: ORBITALMECHANICS_HPP */
