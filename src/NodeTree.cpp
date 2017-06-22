#include "NodeTree.hpp"
#include <iostream>
#include <fstream>

NodeTree::NodeTree()
{
    /// Default constructor
    treeDepth = 0;
    hasChild = false;
}

NodeTree::NodeTree(uint newDepth)
{
    /// Child nodes call this constructor when creating
    treeDepth = newDepth;
    hasChild = false;
}

void NodeTree::generateChildren()
{
  /// Initialize 4 containers, one for each quad of the Node.
  std::vector<Particle*> quad1, quad2, quad3, quad4;
  /// Sort the input bodies into the 4 quadrants.
  for (uint i = 0; i < this->nodeBodies.size(); i++) {
      if (nodeBodies[i]->rx < (nodexLocation + nodeWidth/2)) {
          if (nodeBodies[i]->ry < (nodeyLocation + nodeHeight/2)) {
              quad1.push_back(nodeBodies[i]); /// this is quad 1
          } else {
              quad3.push_back(nodeBodies[i]); /// this is quad 3
          }
      } else {
          if (nodeBodies[i]->ry < (nodeyLocation + nodeHeight/2)) {
              quad2.push_back(nodeBodies[i]); /// this is quad 2
          } else {
              quad4.push_back(nodeBodies[i]); /// this is quad 4
          }
      } 
  } /// end for

  /// Create new children Nodes out of each of the quadrants
  /// by recursively calling this function until a quad contains 1 or 0 particles
  /// NOTE:
  ///  I HAVE USED HEAP ALLOCATED OBJECTS FOR THE SUB NODES. THEY MUST BE DELETED
  ///  THROUGH A LOOP IN THE RESET FUNCTION. THIS MIGHT BE BETTER TO CHANGE ENTIRELY
  ///  AND GO BACK TO STACK POINTERS.
  NodeTree* q1 = new NodeTree(this->treeDepth + 1);
  NodeTree* q2 = new NodeTree(this->treeDepth + 1);
  NodeTree* q3 = new NodeTree(this->treeDepth + 1);
  NodeTree* q4 = new NodeTree(this->treeDepth + 1);

  /// Initialize the child nodes with the width height parameters they contain:
  q1->createSubNode(quad1,nodeWidth/2, nodeHeight/2,nodexLocation, nodeyLocation);
  q2->createSubNode(quad2,nodeWidth/2, nodeHeight/2,nodexLocation+nodeWidth/2, nodeyLocation);
  q3->createSubNode(quad3,nodeWidth/2, nodeHeight/2,nodexLocation, nodeyLocation+nodeHeight/2);
  q4->createSubNode(quad4,nodeWidth/2, nodeHeight/2,nodexLocation+nodeWidth/2, nodeyLocation+nodeHeight/2);


  /// Storing sub nodes inside the ChildNode container.
  ChildNode.push_back(q1);
  ChildNode.push_back(q2);
  ChildNode.push_back(q3);
  ChildNode.push_back(q4);
  ///std::cout << "Depth: " << treeDepth << " Child Size: " << ChildNode.size() << " Bodies: " << quad1.size() + quad2.size() + quad3.size() + quad4.size() << std::endl;
  /// std::cout << std::endl;
  this->hasChild = true;
  /// flag to tell Node that it has children.
}

void NodeTree::createSubNode(Body_ctr childBodies, double width, double height, double xpos, double ypos)
{
    /// Initializing parameters for the current NodeTree.
    // std::cout << "Creating sub node." << std::endl;
    // std::cout << "Bodies entered: " << childBodies.size() << std::endl;
    // std::ofstream nodefile("nodepositions.txt", std::fstream::app);
    // std::ofstream comfile("compositions.txt", std::fstream::app);
    // nodefile << xpos << ' ' << ypos << ' ' << width << ' ' << height << ' ' << childBodies.size() << ' ' << treeDepth <<std::endl;
    // std::vector<double> nodexywh;
    // nodexywh.push_back(xpos); nodexywh.push_back(ypos); nodexywh.push_back(width);
    // nodexywh.push_back(height); nodexywh.push_back(childBodies.size());
    // nodePositions.push_back(nodexywh);
    this->nodeBodies = childBodies;
    this->nodexLocation = xpos;
    this->nodeyLocation = ypos;
    this->nodeWidth = width;
    this->nodeHeight = height;

    double subMass=0;
    double ycenter=0;
    double xcenter=0;

    if (childBodies.size() > 1 ) {
      for (uint i = 0; i < childBodies.size(); i++) {
        subMass += childBodies[i]->mass;
        xcenter += childBodies[i]->rx*childBodies[i]->mass;
        ycenter += childBodies[i]->ry*childBodies[i]->mass;
      }
      /// Set the center of mass of the system.
      nodeMass = subMass;
      nodeCOMx = xcenter/subMass;
      nodeCOMy = ycenter/subMass;
    } else if (childBodies.size() == 1){
      nodeMass = childBodies[0]->mass;
      nodeCOMx = childBodies[0]->rx;
      nodeCOMy = childBodies[0]->ry;
    }
    /// comfile << nodeCOMx << ' ' << nodeCOMy << std::endl;
    ///std::cout << "quad COM: " << nodeCOMx << " " << nodeCOMy << std::endl;
    /// If we currently have more than 1 body in this node, keep calling generateChildren
    /// on the current node until we have 1 or 0 particles in each node.
    if (this->nodeBodies.size() > 1 && treeDepth < maxDepth) {
        generateChildren();
    }
}

void NodeTree::returnTree()
{
  // Function to return the bounding boxes on all nodes in the simulation.
}


void NodeTree::resetNodeTree()
{
    nodeBodies.clear();
    // clearBodies();
    for (size_t i = 0; i < ChildNode.size(); i++) {
        ChildNode[i]->resetNodeTree();
    }
    for (size_t i = 0; i < ChildNode.size(); i++) {
       delete ChildNode[i];
    }
    ChildNode.clear();
    this->hasChild = false;
    ///std::cout << "NodeTree has been completely reset." << std::endl;
}

void NodeTree::clearBodies()
{
  for (size_t i = 0; i < nodeBodies.size(); i++) {
    delete nodeBodies[i];
  }
  nodeBodies.clear();
}

NodeTree::~NodeTree()
{
  // Handle all memory allocated in the std::vector<> containers
  // delete all allocations inside std::vector<NodeTree*>
  resetNodeTree();
  clearBodies();
  // std::cout << "Destroyed Node stucture" << '\n';

}
