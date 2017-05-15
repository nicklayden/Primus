#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP
/*

  Barnes hut particle simulator. Calculates all forces on particles in
  simulation

*/
#include <vector>
#include <limits>
#include <algorithm>

#include "Particle.hpp"
#include "NodeTree.hpp"

typedef std::vector<Particle*> Body_ctr;


class Simulator
{
  public:
    //constructors
    Simulator(NodeTree* node, Body_ctr bodies, float theta, double timestep, double G, float softener, float w, float h, float x, float y);
    ~Simulator();

    //functions
    void CreateNodeTree(double w, double h,double x, double y);
    void CreateNodeTree();
    void ForceSum(bool CheckTimestep);
    void ResetForces();
    void Parameters();
    void RescaleNodeTree();
    void TotalForcesCalculated();
    double TotalKineticEnergy();
    double TotalPotentialEnergy();
    void MergeParticles(Particle* p1, Particle* p2);
    void SimulateForces(bool updatedt);

    //containers/class parameters
    Body_ctr bodies;
    NodeTree* MainNode;
    float theta;
    double timestep;
    double itimestep = std::numeric_limits<double>::max();
    float softener;
    double G;
    double w,h,x,y; //NodeTree size parameters
    uint collisioncounter;
    uint forceCounter = 0;
    double changdtif = 4.0;
    double reducedt = 0.5;
    double hardmindt = 86400;
    int timechanges = 0;


  private:
    void ScanNodeTree(NodeTree* node, Particle* particle);
    void AddForces(Particle* body, NodeTree* node);

};


#endif /* end of include guard: SIMULATOR_HPP */
