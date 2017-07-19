#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "AstroConstants.hpp"

/**
  Particle class container for simulation particles. Currently in 2D
  and each particle has basic properties of a celestial body.


*/

///template<class T>
class Particle
{
    public:
        Particle(){};
        Particle(double rx, double ry, double mass)
        :rx(rx), ry(ry), mass(mass)
        {}
        Particle(double rx, double ry, double mass, double density)
        :rx(rx), ry(ry), mass(mass), density(density)
        {}
        Particle(double rx, double ry, double vx, double vy, double mass, double density)
        :rx(rx), ry(ry),vx(vx),vy(vy),mass(mass),density(density)
        {}

        ~Particle() {}

        void resetAcceleration();
        double rx, ry;
        double vx, vy;
        double KE,PE;
        double ax=0;
        double ay=0;
        double mass;
        double density;
        bool nulled = false;
        short int type=2; //Particle type: 0 null (dead) 1 dust 2 planet 3 star
        // Only calculate forces between particles if the sum of their type is >2.
        // -> no dust-dust calculations.
        long int numForceCalcs;
        double radius=cbrt(3*mass/(4*constants::pi*density));

        bool fixedPosition=false;

        std::string name= "basic";

        friend std::ostream& operator<<(std::ostream& out, const Particle& p)
        {
            /// Overloading << so that we can send print particle properties
            /// in some stream, either to std::cout or a file I guess.
            return out << p.rx << ' ' << p.ry << ' ' << \
                          p.mass << ' ' << p.mass*p.ax << ' ' << p.mass*p.ay <<\
                          " " << p.vx << " " << p.vy << '\n';
        }
};


#endif /* end of include guard: PARTICLE_HPP */
