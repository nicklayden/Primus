/*
  Barnes hut simulator class cpp file

*/

#include <iostream>
#include "simulator.hpp"
#include <cmath>
#include <fstream>
#include <limits>


Simulator::Simulator(NodeTree* node, Body_ctr bodies, float theta, double timestep, double G, float softener, float w, float h, float x, float y)
:bodies(bodies),theta(theta),timestep(timestep), G(G), softener(softener), w(w), h(h), x(x), y(y)
{
    std::cout << "Barnes-Hut Simulator created." << std::endl;
    MainNode = node;
    CreateNodeTree(w,h,x,y);
    collisioncounter = 0;
}

Simulator::~Simulator()
{
    std::cout << "Barnes-Hut Simulator destroyed." << std::endl;
}


void Simulator::ScanNodeTree(NodeTree* node, Particle* particle)
{
  if (node->nodeBodies.size() != 0){
    //Distances 6 flops
    double dx = (node->nodeCOMx - particle->rx);
    double dy = (node->nodeCOMy - particle->ry);
    double distsqr = dx*dx + dy*dy;
    double dist = sqrt(distsqr);

    if (((node->nodeWidth/dist < theta) || (node->hasChild == false)) ) {
      AddForces(particle, node);
    } else {
        // Recursive scan of tree in each child. set loop index to 8
        // for 3d implementation.
      for (size_t k = 0; k < 4; k++) {
        ScanNodeTree(node->ChildNode[k], particle);
      }
    } // end recursive calc
  } // end node size check
}

void Simulator::ForceSum(bool CheckTimestep)
{
  /*Function to call when computing forces for a single timestep.*/
  // std::cout << "Calculating forces..." << std::endl;
  double accel;
  double vel;
  double dt = std::numeric_limits<double>::max();
  double mindt = std::numeric_limits<double>::max();
  #pragma omp parallel for
  for (uint i = 0; i < bodies.size(); i++) {
    ScanNodeTree(MainNode, bodies[i]);

    // calc accel and vel for each particle
    accel = sqrt(bodies[i]->ax*bodies[i]->ax + bodies[i]->ay*bodies[i]->ay);
    vel = sqrt(bodies[i]->vx*bodies[i]->vx + bodies[i]->vy*bodies[i]->vy);

    if (accel == 0) {
      dt = std::numeric_limits<double>::max();
    } else if (vel == 0) {
      dt = 1.5e4/accel;
    } else {
      dt = 1e-1*vel/accel;
    }
    mindt = std::min(mindt,dt);
    mindt = std::max(mindt,hardmindt);
  }
  // update timestep for simulation.
  if (CheckTimestep) {
    if (mindt < itimestep) {
      // std::cout << itimestep << " " << mindt << " " << dt << std::endl;
      itimestep = reducedt*mindt;
      timechanges++;
      // std::cout << "Reduced timestep size to: " << itimestep << std::endl;
    } else if (mindt > changdtif*itimestep) {
      itimestep = reducedt*mindt;
      timechanges++;
      // std::cout << "Increased timestep size to: " << itimestep << std::endl;
    }
  }

}

void Simulator::AddForces(Particle* particle, NodeTree* node)
{
  /*Literally add the forces between particles and nodes on the Quad-Tree*/
  /// calculates the acceleration on the particle from a node COM.
  /// if statement to check if a particle is not itself.
  /// via memory address comparison.

  // Current implementation: calc force only if ptype > 0.
  // "merged" particles have ptype = 0.
  if (particle->type > 0) {
    if (node->nodeBodies.size() == 1 && particle != node->nodeBodies[0]) {

      double dx = particle->rx - node->nodeBodies[0]->rx;
      double dy = particle->ry - node->nodeBodies[0]->ry;
      double distsqr = dx*dx + dy*dy + softener*softener;
      double dist = sqrt(distsqr);

      double acceleration = -node->nodeMass*G/distsqr;
      if (dist <= particle->radius + node->nodeBodies[0]->radius + softener && particle->type > 0 && node->nodeBodies[0]->type > 0) {
        // std::cout << "Collision between " << particle->name << " and " << node->nodeBodies[0]->name << std::endl;
        collisioncounter +=1;
        MergeParticles(particle,node->nodeBodies[0]);
        std::cout << "Merging particles: " << particle->name << " " << node->nodeBodies[0]->name << std::endl;
      }

      particle->ax += acceleration*dx/dist;
      particle->ay += acceleration*dy/dist;
      particle->numForceCalcs++;

    } else if (node->nodeBodies.size() > 1 ) {
      double dx = particle->rx - node->nodeCOMx;
      double dy = particle->ry - node->nodeCOMy;
      double distsqr = dx*dx + dy*dy + softener*softener;
      double dist = sqrt(distsqr);
      double acceleration = -node->nodeMass*G/distsqr;

      particle->ax += acceleration*dx/dist;
      particle->ay += acceleration*dy/dist;
      particle->numForceCalcs++;
    }
  }
}

void Simulator::ResetForces()
{
  /*Reset all attractions on each particle, prepare for next step*/
  for (size_t j = 0; j < bodies.size(); j++) {
    bodies[j]->ax = 0.0;
    bodies[j]->ay = 0.0;
  }
}

void Simulator::Parameters()
{
  std::cout << "Parameters for current simulation:" << std::endl;
  std::cout << "Theta: " << this->theta << std::endl;
  std::cout << "Timestep: " << itimestep << std::endl;
  std::cout << "Softening: " << softener << std::endl;
  std::cout << "G: " << G << std::endl;
  std::cout << "Particle Count: " << bodies.size() << std::endl;
}

void Simulator::CreateNodeTree(double w, double h,double x, double y)
{
  // std::cout << "Test creating node tree" << std::endl;
  MainNode->createSubNode(bodies,w,h,x,y);
  // std::cout << "Node Populated." << std::endl;
}

void Simulator::CreateNodeTree()
{
  MainNode->createSubNode(bodies,this->w,this->h,this->x,this->y);
}

void Simulator::RescaleNodeTree()
{
  double minX = bodies[0]->rx;
  double minY = bodies[0]->ry;
  double maxX = bodies[0]->rx;
  double maxY = bodies[0]->ry;

  for (size_t i = 0; i < bodies.size(); i++) {
    if(bodies[i]->rx < minX) {
      minX = bodies[i]->rx;
    }
    if(bodies[i]->ry < minY) {
      minY = bodies[i]->ry;
    }
    if(bodies[i]->rx > maxX) {
      maxX = bodies[i]->rx;
    }
    if(bodies[i]->ry > maxY) {
      maxY = bodies[i]->ry;
    }
  }
  // std::cout << "Rescaled node values: " << minX << " " << minY << " " << maxX << " " << maxY << std::endl;
  this->w = maxX - minX;
  this->h = maxY - minY;
  this->x = minX;
  this->y = minY;
}

void Simulator::TotalForcesCalculated()
{
  for (size_t i = 0; i < bodies.size(); i++) {
    forceCounter += bodies[i]->numForceCalcs;
    bodies[i]->numForceCalcs = 0;
  }
}

double Simulator::TotalKineticEnergy()
{
  double KE = 0;
  for (size_t i = 0; i < bodies.size(); i++) {
    KE += 0.5*bodies[i]->mass*(bodies[i]->vx*bodies[i]->vx + bodies[i]->vy*bodies[i]->vy);
    // std::cout << bodies[i]->vx*bodies[i]->vx + bodies[i]->vy*bodies[i]->vy << std::endl;
  }
  return KE;
}

double Simulator::TotalPotentialEnergy()
{
  std::ofstream PEsheet("PEsheet.txt");
  double PE = 0;
  double distance,ri,rj;
  for (size_t i = 0; i < bodies.size(); i++) {
    for (size_t j = 0; j < bodies.size(); j++) {
      if (i != j) {
        ri = sqrt(bodies[i]->rx*bodies[i]->rx + bodies[i]->ry*bodies[i]->ry);
        rj = sqrt(bodies[j]->rx*bodies[j]->rx + bodies[j]->ry*bodies[j]->ry);
        distance = sqrt(std::abs(ri - rj));
        PE += -G*bodies[i]->mass*bodies[j]->mass/distance;
      }
    }
    PEsheet << PE << " " << bodies[i]->rx << " " << bodies[i]->ry << "\n";
  }
  return PE;
}

void Simulator::MergeParticles(Particle* p1, Particle* p2)
{
  /*
   Method called when a collision is detected. particles are merged,
   and the smaller particle is destroyed.
   NEEDS:
          Body_ctr resizing (N-1)
          conserve momentum (pure inelastic collision)
          redraw node structure
          move p1 to the barycenter of both particles
          merge into the more massive particle
  */
  double newvx, newvy, newpx, newpy, newmass, newtype;
  newtype = std::max(p1->type, p2->type);

  // combine masses
  newmass = p1->mass + p2->mass;

  // find centre of mass
  newpx = (p1->mass*p1->rx + p2->mass*p2->rx)/newmass;
  newpy = (p1->mass*p1->ry + p2->mass*p2->ry)/newmass;

  // conserve momentum
  newvx = (p1->mass*p1->vx + p2->mass*p2->vx)/newmass;
  newvy = (p1->mass*p1->vy + p2->mass*p2->vy)/newmass;

  // Overwrite p1 properties with the new 'merged' properties
  p1->type = newtype;
  p1->mass = newmass;
  p1->rx = newpx;
  p1->ry = newpy;
  p1->vx = newvx;
  p1->vy = newvy;

  // Delete p2 from the simulation (downgrade to null particle)
  p2->type = 0;
}

void Simulator::SimulateForces(bool updatedt)
{
  RescaleNodeTree();
  ResetForces();
  MainNode->resetNodeTree();
  CreateNodeTree();
  ForceSum(updatedt);
  TotalForcesCalculated();
}
