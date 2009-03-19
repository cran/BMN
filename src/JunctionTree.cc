#include <set>
#include <cstdlib>
#include "JunctionTree.h"
#include <algorithm>
#include <map>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

// which variables are the same in 2 nodes
vector<int> JunctionTree::jointVariables(int node1, int node2)
{
    checkQuittingTime();
    set<int> foo;
    for(unsigned int i=0; i!=tree[node1].variables.size(); ++i)
    {
        foo.insert(tree[node1].variables[i]);
    }
    vector<int> joint;
    for(unsigned int i=0; i!=tree[node2].variables.size(); ++i)
    {
        if(foo.count(tree[node2].variables[i]))
        {
            joint.push_back(tree[node2].variables[i]);
        }
    }
    sort(joint.begin(), joint.end());
    return(joint);
}


// update the information in smallestCliqueForVar given a node
void JunctionTree::updateSmallestCliqueInfo(int node)
{
    checkQuittingTime();
    int newNodeSize = tree[node].variables.size();
    vector<int>::iterator vecIt;
    // go through all variables; check if the size of this new node is smaller than
    // the size of their old nodes; if vector not long enough, extend
    for(vecIt=tree[node].variables.begin(); vecIt!=tree[node].variables.end(); ++vecIt)
    {
        if(((int)smallestCliqueForVar.size())<(*vecIt)+1)
        {
            smallestCliqueForVar.resize((*vecIt)+1,noClique);
        }
        if(smallestCliqueForVar[*vecIt]==noClique)
        {
            smallestCliqueForVar[*vecIt]=node;
        }
        else
        {
            if(newNodeSize < ((int)tree[smallestCliqueForVar[*vecIt]].variables.size()))
            {
                smallestCliqueForVar[*vecIt] = node;
            }
        }
    }
}






// check if a tree node given by a pointer has all incoming mailboxes of children full
bool JunctionTree::childInboxFull(int node)
{
    checkQuittingTime();
    list<int>::iterator listIt= tree[node].in.begin();
    // first check if the node is the rootnode
    if(node!=0) // then the first mailbox goes to the parent
    {
        listIt++;
    }
    
    // check all remaining incoming mailboxes; return false if one is not full
    while(listIt!=tree[node].in.end())
    {
        if(!sepVec[*listIt].full)
        {
            return(false);
        }
        ++listIt;
    }
    return(true);
}

// check if the mailbox of the parent is full
bool JunctionTree::parentInboxFull(int node)
{
    checkQuittingTime();
    if(node==0) // it is the root node
    {
        return(true); // has no parent
    }
    // not the root node
    return(sepVec[tree[node].in.front()].full);
}

// check if this node has already done the collection
bool JunctionTree::nodeCollected(int node)
{
    checkQuittingTime();
    // if it is the root node, it does not need to be collected
    if(node == 0)
    {
        return(true);
    }
    // not the root node
    return(sepVec[tree[node].out.front()].full);
}

// check if the node has already done the distribtion
bool JunctionTree::nodeDistributed(int node)
{
    checkQuittingTime();
    list<int>::iterator listIt= tree[node].out.begin();
    // first check if the node is the rootnode
    if(node!=0) // then the first mailbox goes to the parent
    {
        listIt++;
    }
    
    // check all remaining incoming mailboxes; return false if one is not full
    while(listIt!=tree[node].out.end())
    {
        if(!sepVec[*listIt].full)
        {
            return(false);
        }
    }
    return(true);
}

// collect a single node
void JunctionTree::collectNode(int node)
{
    checkQuittingTime();
    Potential curPot = tree[node].pot;
    // check if it is the root node
    list<int>::iterator listIt=tree[node].in.begin();
    if(node!=0) // first link goes to the parent
    {
        ++listIt;
    }
    // go through all separators from children and include the potentials
    for(;listIt!=tree[node].in.end(); ++listIt)
    {
        curPot.multiply(sepVec[*listIt].pot);
    }
    // now marginalize to the variables of the parent separator, but only if not the root node
    if(node!=0)
    {
    // copy the potential of the current node
        curPot.marginalizeToVars(sepVec[tree[node].out.front()].variables);
        sepVec[tree[node].out.front()].full=true;
        sepVec[tree[node].out.front()].pot=curPot;
    }
    return;
}


// distribute a single node
void JunctionTree::distributeNode(int node)
{
    checkQuittingTime();
    // now send to every child; root only has children, no parents
    list<int>::iterator listIt = tree[node].out.begin();
    if(node!=0)
    {
        ++listIt;
    }
    Potential curPot;
    for(;listIt!=tree[node].out.end(); ++listIt)
    {
        curPot = tree[node].pot;
        // go through all inboxes and multiply the potentials except for the inbox
        // corresponding to the outbox I am currently working on
        list<int>::iterator listItInbox;
        for(listItInbox = tree[node].in.begin(); listItInbox!=tree[node].in.end(); ++listItInbox)
        {
            // check if the inbox corresponds to the current outbox
            if(sepVec[*listItInbox].from != sepVec[*listIt].to)
            {
                curPot.multiply(sepVec[*listItInbox].pot);
            }
        }
        curPot.marginalizeToVars(sepVec[*listIt].variables);
        sepVec[*listIt].full=true;
        sepVec[*listIt].pot = curPot;
    }
    return;
}



// collect the evidence
void JunctionTree::collectEvidence()
{
    checkQuittingTime();
    list<int> nodesToCollectFrom;
    // first walk through all nodes and identify the leaf node
    for(unsigned int node=1; node != tree.size(); ++node)
    {
        if(tree[node].out.size()==1) // only a parent node, therefroe a leaf node
        {
            nodesToCollectFrom.push_back(node);
        }
    }
    
    int curNode;
    // now collect from the nodes
    while(!nodesToCollectFrom.empty())
    {
        // get the front node and delete
         curNode = nodesToCollectFrom.front();
         nodesToCollectFrom.pop_front();
         
         // check if all children inboxes already full, otherwise put the node back into the queue
         if(!childInboxFull(curNode))
         {
            nodesToCollectFrom.push_back(curNode);
         }
         else if(!nodeCollected(curNode))// check if already collected, if not collect
         {
            collectNode(curNode);
            // add the parent 
            nodesToCollectFrom.push_back(sepVec[tree[curNode].out.front()].to);
         }
    }
    return;
}




// distribute the evidence again
void JunctionTree::distributeEvidence()
{
    checkQuittingTime();
    list<int> nodesToDistributeFrom;
    // start node is the root
    nodesToDistributeFrom.push_front(0);
    
    int curNode;
    list<int>::iterator listIt;
    // distribute to the nodes in the list
    while(!nodesToDistributeFrom.empty())
    {
        curNode = nodesToDistributeFrom.front();
        nodesToDistributeFrom.pop_front();
        distributeNode(curNode);
        // add all children to the list
        // go through all nodes; jump over first one unless root
        listIt=tree[curNode].out.begin();
        if(curNode!=0)
        {
            ++listIt;
        }
        while(listIt!=tree[curNode].out.end())
        {
            nodesToDistributeFrom.push_back(sepVec[*listIt].to);
            ++listIt;
        }
    }
}



// take the inboxes of a particular node and mutliply them in
void JunctionTree::includeNodeInboxes(int node)
{
    checkQuittingTime();
    // now walk through all inboxes and multiply the included potential with the node potential
    list<int>::iterator listIt;
    for(listIt=tree[node].in.begin(); listIt!=tree[node].in.end(); ++listIt)
    {
        tree[node].pot.multiply(sepVec[*listIt].pot);
    }
    return;
}


// take all inboxes and multiply with clique potential
void JunctionTree::includeAllInboxes()
{
    checkQuittingTime();
    // go through all nodes in the tree
    unsigned int treeIt;
    for(treeIt = 0; treeIt!=tree.size(); ++treeIt)
    {
        includeNodeInboxes(treeIt);
    }
    return;
}


// run the junction tree algorithm
void JunctionTree::runJTAlgorithm()
{
    collectEvidence();
    distributeEvidence();
    includeAllInboxes();
    return;
}



// copies the tree x; used for copy and assignment constructor
void JunctionTree::copyTree(const JunctionTree& x)
{
    quittingTime = x.quittingTime;
    tree=x.tree;
    sepVec = x.sepVec;
    smallestCliqueForVar = x.smallestCliqueForVar;
}

// check if the quitting time has been reached; if yes, emit an error in R and exit (only place where R is used here)
void JunctionTree::checkQuittingTime()
{
    R_CheckUserInterrupt();
    if(time(NULL) > quittingTime) {
        error("Timeout");
    }
}


JunctionTree::JunctionTree() : tree(0)
{
    quittingTime = -1; // indicates that no time was set
}


// copy constructor
// the copy constructor assumes that the object has been constructed correctly
JunctionTree::JunctionTree(const JunctionTree& x)
{
    copyTree(x);
}

// assignment operator
JunctionTree& JunctionTree::operator=(const JunctionTree& x)
{
    if(this != &x) // different object
    {
        copyTree(x);
    }
    return(*this);
}



// destructor not needed; everything called automatically
JunctionTree::~JunctionTree()
{
}

void JunctionTree::setQuittingTime(time_t quittingTime)
{
    this->quittingTime = quittingTime;
}


// insert a new root node (not connected to anything else)
int JunctionTree::addRoot(const vector<int> &vars, const Potential& rootPot)
{
    checkQuittingTime();
    if(tree.size()>0) // can only add 1 root
    {
        return(-1);
    }
    tree.resize(1);
    tree[0].variables = vars;
    tree[0].pot = rootPot; 
    updateSmallestCliqueInfo(0);
    return(0);
}

// add another node that has a connection to an existing node
int JunctionTree::addNode(int oldNode, const vector<int> &vars, const vector<int> &sepVars, const Potential& newNodePot, const Potential& outSepPot, const Potential& inSepPot)
{
    checkQuittingTime();
    int newNode = tree.size();
    tree.resize(newNode+1);
    tree[newNode].variables = vars;
    tree[newNode].pot = newNodePot;
    
    // generate the 2 separators and connect oldNode to newNode
    int outSep = sepVec.size();
    int inSep=outSep+1;
    sepVec.resize(sepVec.size()+2);
    
    // set the variables in the separator
    sepVec[outSep].variables = sepVars;
    sepVec[inSep].variables = sepVars;
    sepVec[outSep].pot = outSepPot;
    sepVec[inSep].pot = inSepPot;
    sepVec[outSep].full = (outSepPot!=Potential());
    sepVec[inSep].full = (inSepPot!=Potential());
    
    // set from and to nodes in separator
    sepVec[outSep].from = oldNode;
    sepVec[outSep].to = newNode;
    sepVec[inSep].from = newNode;
    sepVec[inSep].to = oldNode;
    
    // connect the nodes to the separators
    tree[oldNode].out.push_back(outSep);
    tree[oldNode].in.push_back(inSep);
    tree[newNode].out.push_back(inSep);
    tree[newNode].in.push_back(outSep);
    
    updateSmallestCliqueInfo(newNode);

    return(newNode);
}


// extract the marginal probabilities for every variable; vector of bools inidicates which were
// variables in the tree; other vector of distributions gives the one-dim result
void JunctionTree::getMarginalProbs(vector<bool> &varInTree, vector<Distribution>& margDists)
{
    // run the algorithm
    runJTAlgorithm();
    
    // walk through the pointers to the smallest cliques
    for(unsigned int i=0; i!=smallestCliqueForVar.size();++i)
    {
        // if clique given, marginalize out and include as distribution
        if(smallestCliqueForVar[i]!=noClique)
        {
            if(varInTree[i])
            {
                throw("Already filled");
            }
            Potential curPot(tree[smallestCliqueForVar[i]].pot);
            vector<int> curVar(1,i); // vector consisting only of the currently treated variable
            // marginalize to this one variable and include in the results
            curPot.marginalizeToVars(curVar);
            curPot.normalize();
            varInTree[i]=curPot.returnOneDimDist(margDists[i]);
        }
    }
    return;
}


// extract the second moment of every pair of variables present in cliques
// smallest cliques will be treated first
void JunctionTree::getMarginalSecMom(vector<vector<bool> > &pairsInTree, vector<vector<double> > &secMomMat, const int graphSize, const vector<bool> &adjMat)
{
    // run the algorithm
    runJTAlgorithm();
    
    // run through the nodes in the tree and for every node, file its number under its size
    // this way, we can start for the smallest clique and work upward
    multimap<int, int> cliqueSizes;
    vector<TreeNode>::iterator treeIt;
    int nodeNum;
    for(treeIt = tree.begin(), nodeNum=0; treeIt !=tree.end(); ++treeIt, ++nodeNum)
    {
        cliqueSizes.insert(pair<int,int>(treeIt->variables.size(), nodeNum));
    }
    
    // now start from the smallest clique
    multimap<int, int>::iterator cliqueIt;
    vector<int> vars;
    int size, index1, index2;
    for(cliqueIt = cliqueSizes.begin(); cliqueIt != cliqueSizes.end(); ++cliqueIt)
    {
        vars = tree[cliqueIt->second].variables;
        size = vars.size();
        
        // go through all pairs in the clique
        for(index1=0; index1<size; ++index1)
        {
            // work on the diagonal, afterwards all the other elements
            if(!pairsInTree[vars[index1]][vars[index1]])
            {
                pairsInTree[vars[index1]][vars[index1]]=true;
                    
                Potential curPot(tree[cliqueIt->second].pot);
                vector<int> margVars(1);
                margVars[0]=vars[index1];
                curPot.marginalizeToVars(margVars);
                curPot.normalize();
                Distribution dist;
                if(!curPot.returnOneDimDist(dist))
                {
                     throw("wasn't marginalized to 1 variable");
                }
                secMomMat[vars[index1]][vars[index1]]= dist[1];
            }
            for(index2 = index1+1; index2< size; ++index2)
            {
                int var1, var2;
                // save the numbers of the currently used variables in increasing order
                if(vars[index1]<vars[index2])
                {
                    var1=vars[index1];
                    var2=vars[index2];
                }
                else
                {
                    var1=vars[index2];
                    var2=vars[index1];
                }
                // if the two nodes have an edge between them in the original graph
                if(adjMat[var1*graphSize+var2])
                {
                    // if the pair has not been calculated before, calculate it now
                    if(!pairsInTree[var1][var2])
                    {
                        pairsInTree[var1][var2]=true;
                        pairsInTree[var2][var1]=true;
                    
                        Potential curPot(tree[cliqueIt->second].pot);
                        vector<int> margVars(2);
                        margVars[0]=var1;
                        margVars[1]=var2;
                        curPot.marginalizeToVars(margVars);
                        curPot.normalize();
                    
                        if(!curPot.returnSecondMoment(secMomMat[var1][var2]))
                        {
                            throw("wasn't marginalized to 2 variables");
                        }
                        secMomMat[var2][var1]=secMomMat[var1][var2];
                    }
                }
            }
        }
    }
}


// calculate just the normalization constant in the current tree
double JunctionTree::getNormalizationConstant()
{
    collectEvidence();
    // only include the inboxes at the root node
    includeNodeInboxes(0);
    // get the normalization constant of the root node and return
    return(tree[0].pot.totalSum());
}




// enter the evidence into the appropriate clique
void JunctionTree::enterEvidence(const int var, const int value)
{
    // check that no overrun
    if(var < (int) smallestCliqueForVar.size())
    {
        if(smallestCliqueForVar[var]!=noClique) // variable included in tree
        {
            tree[smallestCliqueForVar[var]].pot.enterEvidence(var, value);
        }
    }
    return;
}




// print the current junctionTree
void JunctionTree::printTree(ostream &out)
{
    // necessary iterators
    list<int>::iterator sepIt;
    vector<int>::iterator vecIt;
    
    // go through all nodes; print the included variables and the separators
    // also print the potential if it exists
    out << "Junction Tree Size: " << tree.size() << endl;
    for(unsigned int treePos=0; treePos!=tree.size(); ++treePos)
    {
        out << "--------------------- New Node --------------------------" << endl;
        out << "Number:  " << treePos << endl;
        out << "Variables:";
        for(vecIt = tree[treePos].variables.begin(); vecIt!=tree[treePos].variables.end(); ++vecIt)
        {
            out << *vecIt << " ";
        }
        out << endl;
        
        out << "Separators Out:";
        for(sepIt = tree[treePos].out.begin(); sepIt!=tree[treePos].out.end(); ++sepIt)
        {
            out << *sepIt << " ";
        }
        out << endl;
            
        out << "Separators In:";
        for(sepIt = tree[treePos].in.begin(); sepIt!=tree[treePos].in.end(); ++sepIt)
        {
            out << *sepIt << " ";
        }
        out << endl;
            
//        tree[treePos].pot.normalize();
        tree[treePos].pot.print(out);
    }
    // print out which cliques are the smallest cliques for a given variable
    
    out << "---------------------- Separators --------------------------" << endl;
    for(unsigned int sepPos=0; sepPos!=sepVec.size(); ++sepPos)
    {
        out << "Separator Number: " << sepPos << endl;
        if(sepVec[sepPos].full){out << "FULL ";}else{out << "EMPTY";}
        for(vecIt = sepVec[sepPos].variables.begin(); vecIt!=sepVec[sepPos].variables.end(); ++vecIt)
        {
            out << *vecIt << " ";
        }
//        sepVec[sepPos].pot.normalize();
        sepVec[sepPos].pot.print(cout);
        out << endl;
    }
    
    
    out << "----------- Smallest Clique for a given variable ------------" << endl;
    for(unsigned int i=0; i < smallestCliqueForVar.size(); ++i)
    {
        out << i << ": ";
        if(smallestCliqueForVar[i]!=noClique)
        {
            for(unsigned int j=0; j< tree[smallestCliqueForVar[i]].variables.size(); ++j)
            {
                out << tree[smallestCliqueForVar[i]].variables[j] << " ";
            }
        }
        out << endl;
    }
    return;
}



