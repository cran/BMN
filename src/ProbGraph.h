/****************************************************************
***
*** Write a class that saves the structure of a graph,
*** determines the triangulated graph and generates a join tree
*** the graph will be saved as a adjacency matrix as well as in
*** adjacency list form to save time when processing
*** it also saves a matrix of the corresponding thetas
*** the matrix is being saved in standard R form as a vector
***
***
****************************************************************/

#ifndef _PROBGRAPH_
#define _PROBGRAPH_

#include <list>
#include <vector>
#include <iostream>
#include <utility>
#include "JunctionTree.h"
#include <time.h>

using namespace std;

// node/separator data structure
typedef pair< vector<int> , vector<int> > CliqueSepPair;

class ProbGraph
{
public:
    time_t quittingTime;
    vector<bool> adjMat; // the adjacency matrix as a columnwise vector
    vector<double> thetaMat; // the theta matrix as a columnwise vector
    vector<bool> eliminated; // logical vector that saves if a node has been eliminated 
    vector< list<int> > adjList; // the adjacency list of currently not eliminated nodes
    int size;
    list<JunctionTree> jtl;

    // reset the adjList according to what is in the adjMat
    void setAdjList();
    // initialize eliminated to false
    void setEliminated();
    // checks if a node is a simplical node
    bool isSimplical(int nodeNum);
    // find the size of the family of the node 
    int familySize(int nodeNum); 
    // return the current family as a vector (for potentials later)
    vector<int> family(int nodeNum);
    // return the neighbours of the node
    vector<int> neighbours(int nodeNum);
    // fill in all edges necessary to make nodeNum's family a clique into adjMat and adjList
    void fillIn(int nodeNum); 
    // deletes node nodeNum from the adjList and notes in eliminated that it has been eliminated
    void eliminateNode(int nodeNum); 
    // check if all nodes have been eliminated
    bool allEliminated();
    // checks if the second sorted vector is a subset of the first sorted vector
    bool isSubset(const vector<int> &x, const vector<int> &y);

    // find the next node to eliminate
    int findNextElimNode();
    // see if there is a node in the set with neighbours only within the set itself
    // return -1 if no such node is found
    int findNodeWithCertainNeighbours(vector<int> nodeVec);

    // generate a list of nodes/sep that represent the elimination sequence in the 
    // graph. This list will later be used to generate the junction tree.
    list<CliqueSepPair> generateCliqueSepList();
    
   // generate a list of junction trees (more than one if the graph is separable)
   void makeJoinTrees();

    //  initialize the potentials in the tree; numNodes is the number of nodes in the 
    void initializePotentials();


public:
    // initialize the object and generate the adjacency list
    ProbGraph(const int size, const int *graphAdjMat, const double *graphThetaMat, time_t quittingTime); 
    // copy constructor
    ProbGraph(const ProbGraph& x);

    // destructor
    ~ProbGraph();
    
    // collect the marginal probabilities from all trees in the object
    void marginalProbs(vector<bool> &varInTree, vector<Distribution>& margDists);
    
    // enter evidence into all trees
    void enterEvidence(const int var, const int value);
    
    // writes the expectation and covariance matrix. Given pointers have to have sufficient
    // space allocated for them
    void getExpectSecMom(double *expectation, double* secMomMat);
    
    // the same as getExpectSecMom, however only for pairs of variables with an edge between them
    // for the original graph
    void getExpectSecMomActive(double* secMomMat);

    // writes the expectation and covariance matrix. Given pointers have to have sufficient
    // space allocated for them
    void getExpectSecMomSingleVar(int var, double *expectation, double* secMomVec);
    
    // calculate the normalization constant for the model
    double getNormalizationConstant();
    
    // print the probGraph object
    void print(ostream &out);

};






#endif
