/*************************************************
***
*** Implement the junction tree
***
*************************************************/
#ifndef _JUNCTIONTREE_
#define _JUNCTIONTREE_

#include <vector>
#include <list>
#include <iostream>
#include <time.h>
#include "Potential.h"

using namespace std;

const int noClique = -1;

struct Separator
{
    vector<int> variables;
    Potential pot;
    bool full;
    int from, to;
};

struct TreeNode
{
    vector<int> variables;
    Potential pot;
    list<int> out;
    list<int> in;
};

// class that implements the junction tree

class JunctionTree
{
    friend class ProbGraph;
    time_t quittingTime;
    // storing the data
    vector<TreeNode> tree;
    vector<Separator> sepVec;
    // storing the number of the smallest clique with a particular variable
    vector<int> smallestCliqueForVar;
    
    // which variables are the same in 2 nodes
    vector<int> jointVariables(int node1, int node2);

    // update the information in smallestCliqueForVar given a node
    void updateSmallestCliqueInfo(int node);
    
    // check if a tree node given by a pointer has all incoming mailboxes of children full
    bool childInboxFull(int node);
    // check if the mailbox of the parent is full
    bool parentInboxFull(int node);
    // check if this node has already done the collection
    bool nodeCollected(int node);
    // check if the node has already done the distribtion
    bool nodeDistributed(int node);

    // collect a single node
    void collectNode(int node);
    // distribute a single node
    void distributeNode(int node);

    // collect the evidence
    void collectEvidence();
    // distribute the evidence again
    void distributeEvidence();
    // take the inboxes of a particular node and mutliply them in
    void includeNodeInboxes(int node);
    // take all inboxes and multiply with clique potential
    void includeAllInboxes();
    // run the junction tree algorithm
    void runJTAlgorithm();

    // copies the tree x; used for copy and assignment constructor
    void copyTree(const JunctionTree& x);

    // check if the quitting time has been reached; if yes, emit an error in R and exit (only place where R is used here)
    void checkQuittingTime();

public:
    // constructor is only default
    JunctionTree();
    // copy constructor
    JunctionTree(const JunctionTree& x);
    
    // assignment operator
    JunctionTree& operator=(const JunctionTree& x);
    
    // destructor needed to release the separators
    ~JunctionTree();

    // set end time
    void setQuittingTime(time_t quittingTime);
    // insert a new root node (not connected to anything else); returns NULL if there is one already
    int addRoot(const vector<int> &vars, const Potential& rootPot = Potential());
    // add another node that has a connection to an existing node
    int addNode(int oldNode, const vector<int> &vars, const vector<int> &sepVars, const Potential& newNodePot=Potential(), const Potential& outSepPot=Potential(), const Potential& inSepPot=Potential());
    
    // extract the marginal probabilities for every variable; vector of bools inidicates which were
    // variables in the tree; other vector of distributions gives the one-dim result
    // already have to be of sufficient size
    // if element already taken when trying to write in, throws an error
    void getMarginalProbs(vector<bool> &varInTree, vector<Distribution>& margDists);
    
    // extract the second moment of every pair of variables present in cliques
    // smallest cliques will be treated first
    // the input variables have to be large enough to accomodate the data by themselves; 
    // no additional memory will be performed allocated
    // it also calculates the expectation of every node (the diagonal of the second moment matrix
    void getMarginalSecMom(vector<vector<bool> > &pairsInTree, vector<vector<double> > &secMomMat, const int graphSize, const vector<bool> &adjMat);
    
    // calculate just the normalization constant in the current tree
    double getNormalizationConstant();
    
    // enter the evidence into the appropriate clique
    void enterEvidence(const int var, const int value);

    // print the current junctionTree
    void printTree(ostream &out);
};

#endif
