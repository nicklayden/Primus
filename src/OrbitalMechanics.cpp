/*
  Orbital Mechanics equations.
    -Calculating aphelion and perihelion velocities,
    -Energy calculations?
    -Other shit.

*/


#include "OrbitalMechanics.hpp"
#include "planets.hpp"
#include <random>
#include <math.h>
#include <vector>

using namespace constants;


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

void SpawnSolarSystemPlanets(Body_ctr& bodies)
{
  // Particles here NEED to be dynamically allocated!!!!
  // Stack memory allocated here is freed upon exiting the function.
  // That obviously results in a SEGFAULT when the main program tries to use
  // the addresses inside this function.

  // Particle* Moon = new Particle(moon_init_pos,0,0,moon_tangent_vel,moon_mass, moon_density);
  // Moon->name = "Moon";
  using namespace Planets;

  Particle* Earth = new Particle(au,0,0,earth_tangent_vel,earth_mass,earth_density);
  Earth->name = "Earth";

  Particle* Mars = new Particle(mars::aphelion*au,0,0,mars::v_aphelion,mars::mass,2);
  Mars->name = "Mars";

  Particle* Sun = new Particle(0,0,0,0,solar_mass,solar_density);
  // Sun->fixedPosition = true;
  Sun->name = "Sun";

  Particle* Jupiter = new Particle(jupiter_aphelion,0,0,jupiter_v_ap,jupiter_mass,jupiter_density);
  Jupiter->name = "Jupiter";

  bodies.push_back(Sun);
  bodies.push_back(Earth);
  // bodies.push_back(Jupiter);
  // bodies.push_back(Mars);

}

void initParticleDisk(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles)
{
  // initializes particles in a uniform disk between r and r+dr
  // with masses m+dm, and velocities +/- 5000 m/s
  double rx, ry, radius, mass, angle, randomVel, randomVx, randomVy, angle2;
  for (size_t i = 0; i < Nrngparticles; i++) {
      // cout <<"\rSpawning particle " << i+1 << " out of " << Nrngparticles << "   ";
      // cout.flush();
      radius = unirandomval(r, r+dr);
      angle = unirandomval(0,2*pi);
      angle2 = unirandomval(0,2*pi);
      mass = unirandomval(mass_min, mass_max);
      randomVel = 30000*unirandomval(-1,1);
      // randomVx = randomVel*cos(angle2);
      // randomVy = randomVel*sin(angle2);
      randomVx = 0;
      randomVy = 0;
      rx = radius*cos(angle);
      ry = radius*sin(angle);
      Particle* obj = new Particle(rx,ry,randomVx,randomVy,mass,2000);
      bodies.push_back(obj);
   }
}

void makeCircularDisc(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles)
{
  // initializes particles in a uniform disk between r and r+dr
  // with masses m+dm, and velocities +/- 5000 m/s
  double rx, ry, radius, mass, angle, orbvel, vx, vy;
  for (size_t i = 0; i < Nrngparticles; i++) {
      // cout <<"\rSpawning particle " << i+1 << " out of " << Nrngparticles << "   ";
      // cout.flush();
      radius = unirandomval(r, r+dr);
      angle = unirandomval(0,2*pi);
      rx = radius*cos(angle);
      ry = radius*sin(angle);
      orbvel = CircularVelocity(G,solar_mass,radius);
      vx = -orbvel*sin(angle);
      vy = orbvel*cos(angle);
      mass = unirandomval(mass_min, mass_max);
      Particle* obj = new Particle(rx,ry,vx,vy,mass,2000);
      bodies.push_back(obj);
   }
}

void FixSunMomentum(Body_ctr bodies)
{
  // Function that makes sure all initial momentum created in the system is constant
  // So that the sun doesn't drift from the center of the simulation.
  // the sun's momentum should balance all of the bodies in the solar system,
  // so the total momentum is zero.
  // p_sun + sum(p_i) = 0,  i = (1,N) for all N particles minus the sun.
  // thus: p_sun = -sum(p_i)
  // Vector equations:
  // components: m_sun[(v_xs)i + (v_ys)j] = -sum(M_i[(v_xi)i + (v_yi)j])
  double vx_sun = 0;
  double vy_sun = 0;
  std::cout << "Fixing Sun's momentum." << '\n';
  for (size_t i = 0; i < bodies.size()-1; i++) {
    vx_sun -= (bodies[i]->mass * bodies[i]->vx);
    vy_sun -= (bodies[i]->mass * bodies[i]->vy);
  }
  vx_sun /= solar_mass;
  vy_sun /= solar_mass;
  if (bodies[0]->name == "Sun") {
    bodies[0]->vx = vx_sun;
    bodies[0]->vy = vy_sun;
    std::cout << "Sun's new velocity: ";
    std::cout << "v = (" << vx_sun << ")i + (" << vy_sun << ")j" << " (m/s)\n";
  }
}

double unirandomval(double min, double max)
{
  /// This function uses a Mersenne Twister to generate
  /// a uniform real distribution between min and max
  std::random_device rand_device;
  std::mt19937 gen(rand_device());
  std::uniform_real_distribution<double> distr(min, max);
  return distr(gen);
}

double normrandomval(double mean, double variance)
{
  std::random_device rand_device;
  std::mt19937 gen(rand_device());
  std::normal_distribution<double> distr(mean,variance);
  return distr(gen);
}
