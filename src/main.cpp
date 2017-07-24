/*
  MOON COLLISION AND ASTEROID INTERACTION SIMULATOR
  Uses the barnes hut program to simulate asteroid impacts on the moon,
  and the effect that earth has on asteroid trajectories on close flybys.




  IDEA INTRO:
    Have the moon and earth fixed at moon(0,au), and earth(moondist,au)
    Launch asteroids past the earth and moon system, track trajectories.




  Custom install script for packaging the physics/numerics library into
  a single import??

  something like: create a custom main.cpp file, one that can 'swap' in
  or out with the makefile for compilation. ORRRR, precompile it as a header
  ORRRRRRRRR, header library???

  honestly these are all shit ideas, but maybe just for shits...






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
#include "display.hpp"


/*Prototypes -- move to header file*/
void VerletStepOne(Body_ctr bodies, double timestep, uint Nparticles);
void VerletStepTwo(Body_ctr bodies, double timestep, uint Nparticles);
void WindowDrawText(sf::RenderWindow* window, std::vector<sf::Text*> textpointers, std::vector<double> simvals);
void EnterNodeTree(std::vector<std::vector<double> >& nodes, NodeTree* NODE);
void GetBoxDimensions(std::vector<std::vector<double> >& nodes, NodeTree* NODE);
void nullifyMergedParticles(Body_ctr& bodies);
/*Global parameters */

/* Experimental Functions */
// What a stupid fuckin name..
sf::Vector2f CoordTransformSimulationToWindow(Particle* p, sf::RenderWindow* window, double w, double , double x, double y);


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
  double timestep = 86400/3;
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
  bool drawNodes = true;
  bool displayStats = true;
  uint dumpfrequency = 10;
  uint fpsmax = 0;

  std::vector<std::vector<double> > Nodes;


  // Capture command line arguments

  CommandLineSettings(argc, argv,&timestep,&Nrngparticles,&Niterations,&theta,\
                      &softener, &MonitorStats, &dumpfrequency, &fpsmax, &drawNodes);
  /*--------------------------------------------------------------------------*/
  // Gui window parameters & setup
  sf::RenderWindow window(sf::VideoMode(1000,1000), "N-Body Simulator");
  sf::CircleShape temp(1e10);
  sf::View SimView(sf::Vector2f(0,0), sf::Vector2f(w,h));
  sf::View StatsView;
  window.setPosition(sf::Vector2i(0,0));
  window.setView(SimView);
  if (fpsmax != 0) {
    window.setFramerateLimit(fpsmax);
  }


  /*--------------------------------------------------------------------------*/
  // Gui text overlay:
  std::vector<sf::Text*> textpointers;
  std::vector<double> simvalues;
  sf::Text textStep;
  sf::Font font;
  if (!font.loadFromFile("OpenSans-Regular.ttf")) {
    return 1;
  }

  // Display test(font,w,h);


  textStep.setFont(font);
  textStep.setString("Test");
  textStep.setCharacterSize(12);
  textStep.setColor(sf::Color::White);

  sf::Text timer;
  timer.setFont(font);
  timer.setString("0");
  timer.setCharacterSize(12);
  timer.setColor(sf::Color::White);
  timer.setPosition(0,12);

  sf::Text simspeed;
  simspeed.setFont(font);
  simspeed.setString("0");
  simspeed.setCharacterSize(12);
  simspeed.setColor(sf::Color::White);
  simspeed.setPosition(0,24);

  sf::Text simtimetext;
  simtimetext.setFont(font);
  simtimetext.setString("0");
  simtimetext.setCharacterSize(12);
  simtimetext.setColor(sf::Color::White);
  simtimetext.setPosition(0,36);


  sf::Text forcesperstep;
  forcesperstep.setFont(font);
  forcesperstep.setString("0");
  forcesperstep.setCharacterSize(12);
  forcesperstep.setColor(sf::Color::White);
  forcesperstep.setPosition(0,48);

  sf::Text flopsTimer;
  flopsTimer.setFont(font);
  flopsTimer.setString("0");
  flopsTimer.setCharacterSize(12);
  flopsTimer.setColor(sf::Color::White);
  flopsTimer.setPosition(0,60);

  sf::Text currentTimestep;
  currentTimestep.setFont(font);
  currentTimestep.setString("0");
  currentTimestep.setCharacterSize(12);
  currentTimestep.setColor(sf::Color::White);
  currentTimestep.setPosition(0,72);

  sf::Text timeChanges;
  timeChanges.setFont(font);
  timeChanges.setString("0");
  timeChanges.setCharacterSize(12);
  timeChanges.setColor(sf::Color::White);
  timeChanges.setPosition(0,84);

  sf::Text maxchanges;
  maxchanges.setFont(font);
  maxchanges.setString("0");
  maxchanges.setCharacterSize(12);
  maxchanges.setColor(sf::Color::White);
  maxchanges.setPosition(0,96);

  sf::Text earthposition;
  sf::Vector2f earthPosPixels;
  earthposition.setFont(font);
  earthposition.setString("0");
  earthposition.setCharacterSize(12);
  earthposition.setColor(sf::Color::White);
  earthposition.setPosition(0,96);

  /*--------------------------------------------------------------------------*/
  // File streams
  std::ofstream positionOutput("positions.txt");
  std::ofstream EnergyTracking("Energy.txt");

  /*--------------------------------------------------------------------------*/
  // Initialize particles in the container.
  Body_ctr bodies;
  SpawnSolarSystemPlanets(bodies);

  // test particle.
  // Particle* testparticle = new Particle(-3*au,3*au,1e15,2000);
  // double angle = pi/4.0; // angle is 45 degrees at (3,3).
  // double rap = au*sqrt(18.0);
  // std::vector<double> velocity = eccentricityTilt(G,solar_mass,rap,0.9,angle);
  // testparticle->vx = velocity[0];
  // testparticle->vy = velocity[1];
  // std::cout << *testparticle;
  // bodies.push_back(testparticle);
  makeVariableEccentricity(bodies,4*au,1*au,1e15,1e16,0.2,0.3, Nrngparticles);

  // InitDoubleBinarySystem(bodies);
  // makeCircularDisc(bodies,5.0*au,0.2*au,1e15,1e17,Nrngparticles);
  // makeCircularDisc(bodies,1.2*au,3.5*au,1e12,1.1e15,Nrngparticles);
  // initParticleDisk(bodies,2*au,8*au,1e18,1e20,Nrngparticles);
  // Fix solar momentum to conserve.
  FixSunMomentum(bodies);

  // Plots text of a particle's name next to it as it moves.
  // CoordTransformSimulationToWindow(bodies[4],&window,w,h,x,y);


  // // Initialize Simulation object
  Simulator simulation(&GlobalNode,bodies,theta,timestep,G, softener,w,h,x,y);
  simulation.Parameters();
  // SFMLDisplay testwindow(&simulation, &font);
  /*--------------------------------------------------------------------------*/
  auto begintimeclock = Clock::now();
  while (window.isOpen()) {
    //Check window events that have occurred since the last iteration.
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
        // KEYPRESS EVENTS
      } else if (event.type == sf::Event::KeyPressed) {
        // Keyboard handling
          // A key
        if (event.key.code == sf::Keyboard::A && drawNodes == true) {
          drawNodes = false;
        } else if (event.key.code == sf::Keyboard::A && drawNodes == false) {
          drawNodes = true;
          // S key
        } else if (event.key.code == sf::Keyboard::S && displayStats == true) {
          displayStats = false;
        } else if (event.key.code == sf::Keyboard::S && displayStats == false) {
          displayStats = true;
        }
      }
    } // end event check

    /*------------------------------------------------------------------------*/
    // test.MainLoop(simulation.bodies);

    // Main simulation loop.
      auto beginsimloop = Clock::now();
      // Reset view to simulation boundary
      window.setView(SimView);
      window.clear(sf::Color::Black);
      /*-Integration steps-*/
      simulation.SimulateForces(true);
      VerletStepOne(simulation.bodies, timestep, simulation.bodies.size());
      simulation.SimulateForces(false);
      forcesperstep.setString(NumberToString("Forces per step:",simulation.forceCounter,""));
      flop = simulation.forceCounter;
      VerletStepTwo(simulation.bodies, timestep, simulation.bodies.size());
      // Check particles with type=0 and zero all parameters.
      // nullifyMergedParticles(bodies);
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
      textStep.setString(NumberToString("Timestep:",k," "));

      // test.MainLoop(simulation.bodies);

      // Get pixel position of a particle body.
      // earthPosPixels = CoordTransformSimulationToWindow(simulation.bodies[2],&window,w,h,x,y);
      // earthposition.setString(simulation.bodies[2]->name);
      // earthposition.setPosition(earthPosPixels.x,earthPosPixels.y);

      auto endtime = Clock::now();
      auto simtime = sc::duration_cast<std::chrono::milliseconds>(endtime-begintimeclock).count()/1e3;
      auto looptime = sc::duration_cast<std::chrono::milliseconds>(endtime-beginsimloop).count()/1e3;
      simtimer += simulation.itimestep;
      timer.setString(NumberToString("Elapsed time:",simtime,"s"));
      simspeed.setString(NumberToString("Simulation speed:",1./looptime,"steps/s"));
      flopsTimer.setString(NumberToString("Forces/s:",flop/looptime,""));
      simtimetext.setString(NumberToString("sim time:",simtimer,"s"));
      currentTimestep.setString(NumberToString("Current timestep:",simulation.itimestep,"s"));
      timeChanges.setString(NumberToString("Time changes:",simulation.timechanges,""));
      // maxchanges.setString(NumberToString("Suggested max timestep changes (sqrt(steps)):",sqrt(k),""));
      // View the statistics on screen.
      if (displayStats) {
        window.setView(StatsView);
        // window.draw(maxchanges);
        window.draw(currentTimestep);
        window.draw(timeChanges);
        window.draw(timer);
        window.draw(simtimetext);
        window.draw(simspeed);
        window.draw(textStep);
        window.draw(forcesperstep);
        window.draw(flopsTimer);
        // window.draw(earthposition);
      }
      window.display();
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
  // if (NODE->nodeBodies.size() == (1 | 0))  {
    nodes.push_back(temp);
  // }
  if ( NODE->treeDepth < 60) {
    EnterNodeTree(nodes, NODE);
  }
}


void nullifyMergedParticles(Body_ctr& bodies)
{
  for (size_t i = 0; i < bodies.size(); i++) {
    if (bodies[i]->type == 0 && bodies[i]->nulled == false) {
      bodies[i]->vx = 0;
      bodies[i]->vy = 0;
      bodies[i]->rx = 0;
      bodies[i]->ry = 0;
      bodies[i]->nulled = true;
      std::cout << "Nullified Particle: " << bodies[i] << '\n';
    }
  }
}

sf::Vector2f CoordTransformSimulationToWindow(Particle* p, sf::RenderWindow* window, double w, double h, double x, double y)
{
  /* Transform coordinates from simulation coords (in physical units) to render window coordinates
     (pixels)
     Things we need:
        - Position of the particle in question.
        - The width and height of the initial simulation. NOT CURRENT WINDOW PARAMETERS
          so, DON'T use w,h,x,y from "sim".
        - The parameters of the drawing window (width and height)
  */
  // positions of the particle in the SIMULATION
  sf::Vector2f simP(p->rx, p->ry);
  sf::Vector2f winP; // Coordinates of the particle in window coord system.(pixels)

  sf::Vector2u winSize = window->getSize();

  float pixperaux = winSize.x/w; // number of pixels per AU in X dimension
  float pixperauy = winSize.y/h; // number of pixels per AU in Y dimension

  winP.x = 0.5*winSize.x + pixperaux*simP.x;
  winP.y = 0.5*winSize.y + pixperauy*simP.y;


  // std::cout << "Window size: " << winSize.x << " " << winSize.y << std::endl;
  // std::cout << "Particle sim position: " << simP.x << " " << simP.y << std::endl;
  // std::cout << "Particle win position: " << winP.x << " " << winP.y << std::endl;

  return winP;
}
