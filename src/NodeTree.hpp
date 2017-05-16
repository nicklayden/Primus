#ifndef NODETREE_HPP
#define NODETREE_HPP
#include "Particle.hpp"
#include <vector>

/*
  NodeTree:
            Takes in a container of particles and creates a node/leaf tree
            based on the Barnes-Hut algorithm. All particles are subdivided
            into boxed regions called 'nodes' until there exists one or
            zero particles in each node throughout the tree. The relative
            positions of the particles and nodes can then be used to
            simplify the calculation of forces for NBody programs.

  Capabilities:
            Currently, this program is set up for a 2D implementation of the
            algorithm.

*/

/// container for a collection of bodies.
typedef std::vector<Particle*> Body_ctr;
typedef unsigned int uint;

class NodeTree
{
    public:

        /// Member functions of NodeTree class
        NodeTree();
        NodeTree(uint newDepth);
        ~NodeTree();
        void generateChildren(); /// creates 4(8) children to store particles
        void createSubNode(Body_ctr childBodies, double width, double height, double xpos=0, double ypos=0);
        void resetNodeTree(); /// Recursive reset of Node and all children inside
        void returnTree();
        void clearBodies();

        /// container for all particles in simulation.
        Body_ctr nodeBodies;
        /// recursively generated child node structure.
        std::vector<NodeTree*> ChildNode;

        /// Member variables of NodeTree class
        bool hasChild; /// flag to check if a current node has children or not
        uint treeDepth; /// current depth into the tree.

        double nodexLocation; /// location of bottom left corner of node
        double nodeyLocation; /// location of bottom left corner of node
        double nodeWidth; /// width of current node.
        double nodeHeight; /// height of current node. equivalent to nodeWidth
        double nodeMass; /// Total mass contained within current node.
        double nodeCOMx; /// center of mass of current node in x
        double nodeCOMy; /// center of mass of current node in y
        std::vector<std::vector<double> > nodePositions;
        friend std::ostream& operator<<(std::ostream& out, const NodeTree& p)
        {
            /// Overloading << so that we can send print particle properties
            /// in some stream, either to std::cout or a file I guess.
            return out << p.nodexLocation << ' ' << p.nodeyLocation << ' ' << \
                          p.nodeWidth << " " << p.nodeHeight << std::endl;
        }
    private:
        uint maxDepth = 60; /// Recursion limit for generating trees.
                            /// Set private so any possible inheriting class
                            /// cannot alter this value.

};


#endif /* end of include guard: NODETREE_HPP */
