/*
  MOON COLLISION AND ASTEROID INTERACTION SIMULATOR
  Uses the barnes hut program to simulate asteroid impacts on the moon,
  and the effect that earth has on asteroid trajectories on close flybys.




  IDEA INTRO:
    Have the moon and earth fixed at moon(0,au), and earth(moondist,au)
    Launch asteroids past the earth and moon system, track trajectories.

*/
/*Standard libraries*/
#include <vector>
#include <random>
#include <chrono>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

/*Graphics library SFML*/
#include <SFML/Graphics.hpp>



/*Custom libraries and modules*/
#include "AstroConstants.hpp"
#include "planets.hpp"
#include "NodeTree.hpp"
#include "Particle.hpp"
#include "simulator.hpp"
#include "commandline.hpp"
#include "gui.hpp"
#include "OrbitalMechanics.hpp"


/*Prototypes -- move to header file*/
double unirandomval(double min, double max);
double normrandomval(double mean, double variance);
void VerletStepOne(Body_ctr bodies, double timestep, uint Nparticles);
void VerletStepTwo(Body_ctr bodies, double timestep, uint Nparticles);
void initParticleDisk(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles);
void InitDoubleBinarySystem(Body_ctr& bodies);
void SpawnSolarSystemPlanets(Body_ctr& bodies);
void singleLeap(Body_ctr bodies, double timestep);
void makeCircularDisc(Body_ctr& bodies, double r, double dr, double mass_min, double mass_max, uint Nrngparticles);
void makeCircularDisc2(Body_ctr& bodies, double r, double dr, double discmass, uint Nrngparticles);
void infoPointers(std::vector<sf::Text*>& textptrs, sf::Font* font);
void WindowDrawText(sf::RenderWindow* window, std::vector<sf::Text*> textpointers, std::vector<double> simvals);
void EnterNodeTree(std::vector<std::vector<double> >& nodes, NodeTree* NODE);
void GetBoxDimensions(std::vector<std::vector<double> >& nodes, NodeTree* NODE);
void FixSunMomentum(Body_ctr bodies);
/*Global parameters */


/*Namespace contamination*/
using namespace constants;

/*Typedefs*/
typedef std::chrono::high_resolution_clock Clock;
namespace sc = std::chrono;

int main(int argc, char** argv)
{
  // Setup simulation parameters. DEFAULT VALUES!
  NodeTree GlobalNode;
  float theta = 0.5;
  double timestep = 86400;
  double G = mks_universal_G;
  float softener = 100000;
  double TotalEnergy,PE,KE;
  // uint Nparticles = 4;
  uint Nrngparticles = 400;
  uint Niterations = 1000;
  uint k = 0;
  double flop;
  double simtimer=0;
  // Initial Node Box bounds. Could possibly just be simulation hard limits
  double w = 30*au;
  double h = 30*au;
  double x = -15*au;
  double y = -15*au;

  bool MonitorStats = false;
  bool drawNodes = false;
  uint dumpfrequency = 10;
  uint fpsmax = 0;

  std::vector<std::vector<double> > Nodes;


  // Capture command line arguments
  CommandLineSettings(argc, argv,&timestep,&Nrngparticles,&Niterations,&theta,\
                      &softener, &MonitorStats, &dumpfrequency, &fpsmax, &drawNodes);
  /*--------------------------------------------------------------------------*/
  // Gui window parameters & setup
  sf::RenderWindow window(sf::VideoMode(1080,1080), "N-Body Simulator");
  sf::CircleShape temp(1e10);
  sf::View SimView(sf::Vector2f(0,0), sf::Vector2f(w/au,h/au));
  window.setPosition(sf::Vector2i(0,0));
  window.setView(SimView);
  if (fpsmax != 0) {
    window.setFramerateLimit(fpsmax);
  }
  // Window 2 for text:
  sf::RenderWindow window2(sf::VideoMode(500,150), "Real-time Parameters");
  window2.setPosition(sf::Vector2i(800,0));
  // sf::View SimView2(sf::Vector2f(0,0), sf::Vector2f(100,100));
  // window2.setView(SimView2);
  // window2.setFramerateLimit(60);


  /*--------------------------------------------------------------------------*/
  // Gui text overlay:
  std::vector<sf::Text*> textpointers;
  std::vector<double> simvalues;
  sf::Text textStep;
  sf::Font font;
  if (!font.loadFromFile("OpenSans-Regular.ttf")) {
    return 1;
  }

  infoPointers(textpointers,&font);

  textStep.setFont(font);
  textStep.setString("Test");
  textStep.setCharacterSize(12);
  textStep.setColor(sf::Color::Red);

  sf::Text timer;
  timer.setFont(font);
  timer.setString("0");
  timer.setCharacterSize(12);
  timer.setColor(sf::Color::Red);
  timer.setPosition(0,12);

  sf::Text simspeed;
  simspeed.setFont(font);
  simspeed.setString("0");
  simspeed.setCharacterSize(12);
  simspeed.setColor(sf::Color::Red);
  simspeed.setPosition(0,24);

  sf::Text simtimetext;
  simtimetext.setFont(font);
  simtimetext.setString("0");
  simtimetext.setCharacterSize(12);
  simtimetext.setColor(sf::Color::Red);
  simtimetext.setPosition(0,36);


  sf::Text forcesperstep;
  forcesperstep.setFont(font);
  forcesperstep.setString("0");
  forcesperstep.setCharacterSize(12);
  forcesperstep.setColor(sf::Color::Red);
  forcesperstep.setPosition(0,48);

  sf::Text flopsTimer;
  flopsTimer.setFont(font);
  flopsTimer.setString("0");
  flopsTimer.setCharacterSize(12);
  flopsTimer.setColor(sf::Color::Red);
  flopsTimer.setPosition(0,60);

  sf::Text currentTimestep;
  currentTimestep.setFont(font);
  currentTimestep.setString("0");
  currentTimestep.setCharacterSize(12);
  currentTimestep.setColor(sf::Color::Red);
  currentTimestep.setPosition(0,72);

  sf::Text timeChanges;
  timeChanges.setFont(font);
  timeChanges.setString("0");
  timeChanges.setCharacterSize(12);
  timeChanges.setColor(sf::Color::Red);
  timeChanges.setPosition(0,84);

  sf::Text maxchanges;
  maxchanges.setFont(font);
  maxchanges.setString("0");
  maxchanges.setCharacterSize(12);
  maxchanges.setColor(sf::Color::Red);
  maxchanges.setPosition(0,96);

  /*--------------------------------------------------------------------------*/
  // File streams
  std::ofstream positionOutput("positions.txt");
  std::ofstream EnergyTracking("Energy.txt");

  /*--------------------------------------------------------------------------*/
  // Initialize particles in the container.
  Body_ctr bodies;
  SpawnSolarSystemPlanets(bodies);
  // InitDoubleBinarySystem(bodies);
  makeCircularDisc(bodies,5.1*au,0.2*au,1e15,1e18,Nrngparticles);
  makeCircularDisc(bodies,1.5*au,3.6*au,1e18,1.1e20,Nrngparticles);
  // initParticleDisk(bodies,7*au,1*au,1e18,1e20,Nrngparticles);
  // Fix solar momentum to conserve.
  FixSunMomentum(bodies);

  // // Initialize Simulation object
  Simulator simulation(&GlobalNode,bodies,theta,timestep,G, softener,w,h,x,y);
  simulation.Parameters();
  /*--------------------------------------------------------------------------*/
  /*
    IMPORTANT: Currently, the for loop in the event loop doesn't do anything.
    The simulation runs continuously until the window is closed.
    --Shitty fix with a break statement when k increments past Niterations.
  */
  auto begintimeclock = Clock::now();
  while (window.isOpen() && window2.isOpen()) {
    //Check window events that have occurred since the last iteration.
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      } else if (event.type == sf::Event::KeyPressed) {
        if (event.key.code == sf::Keyboard::A && drawNodes == true) {
          drawNodes = false;
        } else if (event.key.code == sf::Keyboard::A && drawNodes == false) {
          drawNodes = true;
        }
      }



    } // end event check
    while (window2.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window2.close();
      }
    } // end event check
    /*------------------------------------------------------------------------*/

    // Main simulation loop.
      auto beginsimloop = Clock::now();
      window.clear(sf::Color::Black);
      window2.clear(sf::Color::Black);
      /*-Integration steps-*/
      simulation.SimulateForces(true);
      VerletStepOne(simulation.bodies, simulation.itimestep, simulation.bodies.size());
      simulation.SimulateForces(false);
      forcesperstep.setString(NumberToString("Forces per step:",simulation.forceCounter,""));
      flop = simulation.forceCounter;
      VerletStepTwo(simulation.bodies, simulation.itimestep, simulation.bodies.size());
      // Energy Monitoring
      if ( MonitorStats && k%dumpfrequency == 0) {
        KE = simulation.TotalKineticEnergy();
        PE = simulation.TotalPotentialEnergy();
        TotalEnergy = KE + PE;
        EnergyTracking << TotalEnergy << " " << k  << " " << KE << " " << PE << std::endl;
      }
      simulation.TotalForcesCalculated();

      /*-Drawing steps-*/
      // std::cout << "Number of Nodes: " << Nodes.size() << '\n';
      if (drawNodes){
        // Get node box dimensions
        EnterNodeTree(Nodes, &GlobalNode);
        DrawBoxes(&window, Nodes);
      }
      DrawParticles(&window,simulation.bodies);
      // window.draw(textStep);
      window.display();
      textStep.setString(NumberToString("Timestep:",k," "));

      auto endtime = Clock::now();
      auto simtime = sc::duration_cast<std::chrono::milliseconds>(endtime-begintimeclock).count()/1e3;
      auto looptime = sc::duration_cast<std::chrono::milliseconds>(endtime-beginsimloop).count()/1e3;
      simtimer += simulation.itimestep;
      timer.setString(NumberToString("Elapsed time:",simtime,"s"));
      simspeed.setString(NumberToString("Simulation speed:",1./looptime,"steps/s"));
      flopsTimer.setString(NumberToString("Forces/s:",flop/looptime,""));
      simtimetext.setString(NumberToString("sim time:",simtimer,"s"));
      currentTimestep.setString(NumberToString("Current timestep:",simulation.itimestep,"s"));
      timeChanges.setString(NumberToString("Time changeS:",simulation.timechanges,""));
      maxchanges.setString(NumberToString("Max changes (sqrt(steps)):",sqrt(k),""));
      window2.draw(maxchanges);
      window2.draw(currentTimestep);
      window2.draw(timeChanges);
      window2.draw(timer);
      window2.draw(simtimetext);
      window2.draw(simspeed);
      window2.draw(textStep);
      window2.draw(forcesperstep);
      window2.draw(flopsTimer);
      // WindowDrawText(&window2,textpointers,simvalues);
      window2.display();
      if (drawNodes) {
        Nodes.clear();
      }
      k++;
      simulation.forceCounter = 0;
  }
  // auto endtime = Clock::now();
  // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-begintimeclock).count()/1e3;
  std::cout << "Collisions: " << simulation.collisioncounter << std::endl;

  std::cout << "Destroying allocated Particles." << '\n';
  for (size_t i = 0; i < bodies.size(); i++) {
    delete bodies[i];
  }
}



void VerletStepOne(Body_ctr bodies, double timestep, uint Nparticles)
{
    // Function for updating the position in the Verlet Integration Scheme
    // MUST BE EXECUTED IN A PAIR WITH VerletStepTwo.
    // #pragma omp parallel for
    for (size_t k = 0; k < bodies.size(); k++) {
      if (bodies[k]->fixedPosition == false){
        bodies[k]->vx += 0.5*bodies[k]->ax*timestep; // v(t+0.5dt) += 0.5*a(t)*dt
        bodies[k]->vy += 0.5*bodies[k]->ay*timestep; // same
        bodies[k]->rx += bodies[k]->vx*timestep; // x(t+dt) += v(t+0.5dt)*dt
        bodies[k]->ry += bodies[k]->vy*timestep; // same
      }
    } // after this loop, all positions are now x(t + dt), velocities unchanged
}

void VerletStepTwo(Body_ctr bodies, double timestep, uint Nparticles)
{
    // Function for updating the velocity in the Verlet Integration Scheme
    // MUST BE EXECUTED IN A PAIR WITH VerletStepOne.
    // #pragma omp parallel for
    for (size_t k = 0; k < bodies.size(); k++) {
      if (bodies[k]->fixedPosition == false){
        bodies[k]->vx += 0.5*bodies[k]->ax*timestep;
        bodies[k]->vy += 0.5*bodies[k]->ay*timestep;
      }
    } // now velocities have been updated. this is the end of the integration step.
}

void singleLeap(Body_ctr bodies, double timestep)
{
  for (size_t i = 0; i < bodies.size(); i++) {
    bodies[i]->vx += bodies[i]->ax*timestep;
    bodies[i]->vy += bodies[i]->ay*timestep;
    bodies[i]->rx += bodies[i]->vx*timestep;
    bodies[i]->ry += bodies[i]->vy*timestep;
  }
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
      randomVx = randomVel*cos(angle2);
      randomVy = randomVel*sin(angle2);
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

// void makeCircularDisc2(Body_ctr& bodies, double r, double dr, double discmass, uint Nrngparticles)
// {
//   double massAvg = 1e15, massdev;
//   double rx, ry, radius, mass, angle, orbvel, vx, vy;
//   for (size_t i = 0; i < Nrngparticles; i++) {
//       radius = unirandomval(r, r+dr);
//       angle = unirandomval(0,2*pi);
//       rx = radius*cos(angle);
//       ry = radius*sin(angle);
//       orbvel = CircularVelocity(G,solar_mass,radius);
//       vx = -orbvel*sin(angle);
//       vy = orbvel*cos(angle);
//       mass = unirandomval(massAvg, massdev);
//       Particle* obj = new Particle(rx,ry,vx,vy,mass,20000);
//       bodies.push_back(obj);
//    }
// }

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
  bodies.push_back(Jupiter);
  bodies.push_back(Mars);

}

void InitDoubleBinarySystem(Body_ctr& bodies)
{
  // Creates a system where a binary star system is bound to another binary
  // star system. Thus making this, a quadrinary(?) star system.
  // Got this idea from Chris MacMackin's program A'tuin.
  Particle* s1 = new Particle(au,-5e-2*au,21061,66601,solar_mass,2);
  s1->name = "s1";
  std::cout << "radius: " << s1->radius << std::endl;
  Particle* s2 = new Particle(au,5e-2*au,-21061,66601,solar_mass,2);
  s2->name = "s2";
  Particle* s3 = new Particle(-au,-5e-2*au,21061,-66601,solar_mass,2);
  s3->name = "s3";
  Particle* s4 = new Particle(-au,5e-2*au,-21061,-66601,solar_mass,2);
  s4->name = "s4";
  std::cout << "particles: " << bodies.size() << std::endl;
  // bodies.push_back(Sun);
  bodies.push_back(s1);
  bodies.push_back(s2);
  bodies.push_back(s3);
  bodies.push_back(s4);
}

void infoPointers(std::vector<sf::Text*>& textptrs, sf::Font* font)
{
  // Create pointers to text for the simulation stats window.
  // append to a vector full of sf::Text pointers.

  // sf::Font* font;
  // if (!font.loadFromFile("OpenSans-Regular.ttf")) {
  //   // return 1;
  // }
  sf::Text* textpointer = new sf::Text;
  textpointer->setFont(*font);
  textpointer->setString("Test value");
  textpointer->setCharacterSize(12);
  textpointer->setColor(sf::Color::Red);
  textpointer->setPosition(0,90);

  // sf::Text* textStep = new sf::Text;
  // textStep.setFont(*font);
  // textStep.setString("Test");
  // textStep.setCharacterSize(12);
  // textStep.setColor(sf::Color::Red);
  //
  // sf::Text* timer = new sf::Text;
  // timer->setFont(*font);
  // timer->setString("0");
  // timer->setCharacterSize(12);
  // timer->setColor(sf::Color::Red);
  // timer->setPosition(0,12);
  //
  // sf::Text* simspeed = new sf::Text;
  // simspeed->setFont(*font);
  // simspeed->setString("0");
  // simspeed->setCharacterSize(12);
  // simspeed->setColor(sf::Color::Red);
  // simspeed->setPosition(0,24);
  //
  // sf::Text* simtimetext = new sf::Text;
  // simtimetext->setFont(*font);
  // simtimetext->setString("0");
  // simtimetext->setCharacterSize(12);
  // simtimetext->setColor(sf::Color::Red);
  // simtimetext->setPosition(0,36);
  //
  //
  // sf::Text* forcesperstep = new sf::Text;
  // forcesperstep->setFont(*font);
  // forcesperstep->setString("0");
  // forcesperstep->setCharacterSize(12);
  // forcesperstep->setColor(sf::Color::Red);
  // forcesperstep->setPosition(0,48);
  //
  // sf::Text* flopsTimer = new sf::Text;
  // flopsTimer->setFont(*font);
  // flopsTimer->setString("0");
  // flopsTimer->setCharacterSize(12);
  // flopsTimer->setColor(sf::Color::Red);
  // flopsTimer->setPosition(0,60);
  //
  // sf::Text* currentTimestep = new sf::Text;
  // currentTimestep->setFont(*font);
  // currentTimestep->setString("0");
  // currentTimestep->setCharacterSize(12);
  // currentTimestep->setColor(sf::Color::Red);
  // currentTimestep->setPosition(0,72);
  //
  // textptrs.push_back(textStep);
  // textptrs.push_back(timer);
  // textptrs.push_back(simspeed);
  // textptrs.push_back(simtimetext);
  // textptrs.push_back(forcesperstep);
  // textptrs.push_back(flopsTimer);
  // textptrs.push_back(currentTimestep);
  textptrs.push_back(textpointer);
}

void WindowDrawText(sf::RenderWindow* window, std::vector<sf::Text*> textpointers, std::vector<double> simvals)
{
  // take in simulation values and text drawn on second window and draw it.
  for (size_t i = 0; i < textpointers.size(); i++) {
    // *textpointers[i] = simvalues
    window->draw(*textpointers[i]);
  }
}
void EnterNodeTree(std::vector<std::vector<double> >& nodes, NodeTree* NODE)
{
  //  Re-implement as member to NodeTree? Probably.
  // Enter from the toplevel node, and recursively add each Box's dimensions
  // to the nodes vector.
  for (size_t i = 0; i < NODE->ChildNode.size(); i++) {
    GetBoxDimensions(nodes,NODE->ChildNode[i]);
  }
}

void GetBoxDimensions(std::vector<std::vector<double> >& nodes, NodeTree* NODE)
{
  // Gets the dimensions of every Node box in the Barnes Hut Tree.
  std::vector<double> temp;
  temp.push_back(NODE->nodexLocation);
  temp.push_back(NODE->nodeyLocation);
  temp.push_back(NODE->nodeWidth);
  temp.push_back(NODE->nodeHeight);
  // Only pick out external nodes for now. (Boxes around particles only)
  // if (NODE->nodeBodies.size() >= 1) {
    nodes.push_back(temp);
  // }
  if ( NODE->treeDepth < 60) {
    EnterNodeTree(nodes, NODE);
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
