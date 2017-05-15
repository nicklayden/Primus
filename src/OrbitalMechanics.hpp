#ifndef ORBITALMECHANICS_HPP
#define ORBITALMECHANICS_HPP
#include <vector>

double perihelionVelocity(double G, double msun, double rperi, double e);
double aphelionVelocity(double G, double msun, double rap, double e);
std::vector<double> tiltedAphelion(double vap, double angle);
std::vector<double> tiltedPerihelion(double vperi, double angle);

double CircularVelocity(double G, double msun, double r);

#endif /* end of include guard: ORBITALMECHANICS_HPP */
