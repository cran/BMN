#include "ProbGraph.h"
#include <algorithm>

// reset the adjList according to what is in the adjMat
void ProbGraph::setAdjList()
{
    for(int i=0; i!=size; ++i)
    {
        adjList[i].clear();
    }
    for(int i=0; i!=size; ++i)
    {
        for(int j=i+1; j<size; ++j)
        {
            if(adjMat[i*size+j]>0)
            {
                adjList[i].push_back(j);
                adjList[j].push_back(i);
            }
        }
    }
    return;
}


// initialize eliminated to false
void ProbGraph::setEliminated()
{
    for(int i=0; i!=size; ++i)
    {
        eliminated[i]=false;
    }
    return;
}






// checks if a node is a simplical node
bool ProbGraph::isSimplical(int nodeNum)
{
    list<int>::iterator listItInner, listItOuter;
    
    // check that the node still exists
    if(eliminated[nodeNum])
    {
        throw("Node already eliminated in isSimplical");
    }
    
    // walk through all neighbours of nodeNum and compare if they are neighbours of each other
    // note that we assume an undirected graph
    for(listItOuter = adjList[nodeNum].begin(); listItOuter!=adjList[nodeNum].end(); ++listItOuter)
    {
        for((listItInner = listItOuter)++; listItInner!=adjList[nodeNum].end(); ++listItInner)
        {
            if(!adjMat[*listItInner * size + *listItOuter])
            {
                return(false);
            }
        }
    }
    return(true);
}



// find the size of the family of the node
int ProbGraph::familySize(int nodeNum)
{
    // check that the node still exists
    if(eliminated[nodeNum])
    {
        throw("Node already eliminated in familySize");
    }
    return(adjList[nodeNum].size()+1);
}


// return the current family as a vector (for potentials later)
vector<int> ProbGraph::family(int nodeNum)
{
    vector<int> foo = neighbours(nodeNum);
    foo.push_back(nodeNum);
    return(foo);
}


// return the current family as a vector (for potentials later)
vector<int> ProbGraph::neighbours(int nodeNum)
{
    vector<int> foo;
    list<int>::iterator listIt;
    for(listIt = adjList[nodeNum].begin(); listIt!=adjList[nodeNum].end(); ++listIt)
    {
        foo.push_back(*listIt);
    }
    return(foo);
}





// fill in all edges necessary to make nodeNum's family a clique
void ProbGraph::fillIn(int nodeNum)
{
    list<int>::iterator listItInner, listItOuter;
    
        // check that the node still exists
    if(eliminated[nodeNum])
    {
        throw("Node already eliminated in fillIn");
    }

    // walk through all neighbours of nodeNum and compare if they are neighbours of each other
    // note that we assume that we have an undirected graph; fill in missing edges to make clique
    for(listItOuter = adjList[nodeNum].begin(); listItOuter!=adjList[nodeNum].end(); ++listItOuter)
    {
        for((listItInner = listItOuter)++; listItInner!=adjList[nodeNum].end(); ++listItInner)
        {
            if(!adjMat[*listItInner * size + *listItOuter])
            {
                // fill in edges in adjacency matrix
                adjMat[*listItInner * size + *listItOuter]=true;
                adjMat[*listItOuter * size + *listItInner]=true;
                // fill in edges in adjacency list
                adjList[*listItInner].push_back(*listItOuter);
                adjList[*listItOuter].push_back(*listItInner);
            }
        }
    }
    return;
}


// deletes node nodeNum from the adjList and notes in eliminated that it has been eliminated
void ProbGraph::eliminateNode(int nodeNum)
{
    if(!eliminated[nodeNum])
    {
        eliminated[nodeNum]=true;
        
        // go through all neighbours of the node that is to be eliminated
        list<int>::iterator listItOuter, listItInner;
        for(listItOuter = adjList[nodeNum].begin(); listItOuter!=adjList[nodeNum].end(); ++listItOuter)
        {
            // find the edge to nodeNum and eliminate
            for(listItInner = adjList[*listItOuter].begin(); listItInner!=adjList[*listItOuter].end(); ++listItInner)
            {
                if(*listItInner==nodeNum)
                {
                    adjList[*listItOuter].erase(listItInner);
                    break;
                }
            }
        }
    }
    return;
}

// check if all nodes have been eliminated
bool ProbGraph::allEliminated()
{
    for(int i=0; i!=size; ++i)
    {
        if(!eliminated[i])
        {
            return(false);
        }
    }
    return(true);
}



// checks if the second sorted vector is a subset of the first sorted vector
bool ProbGraph::isSubset(const vector<int> &x, const vector<int> &y)
{
    vector<int>::const_iterator vecItX = x.begin();
    vector<int>::const_iterator vecItY;
    
    for(vecItY=y.begin(); vecItY!=y.end(); ++vecItY)
    {
        if(vecItX == x.end()) // no items in x to search through left
        {
            return(false);
        }
        while(*vecItY != *vecItX) // search for current item in Y
        {
            ++vecItX; // not found, increase X
            if(vecItX == x.end()) // at end, item not included
            {
                return(false);
            }
        }
        // item found
        ++vecItX;
    }
    return(true);
}







// find the next node to eliminate
int ProbGraph::findNextElimNode()
{
    int curFamSize = size; // larger than maximal possible size
    int curElimNode = -1;
    for(int i=0; i!=size; ++i)
    {
        if(!eliminated[i])
        {
            if(isSimplical(i))
            {
                return(i);
            }
            else
            {
                int foo = familySize(i);
                if(foo < curFamSize)
                {
                    curFamSize = foo;
                    curElimNode = i;
                }
            }
        }
    }
    return(curElimNode);
}



// see if there is a node in the set with neighbours only within the set itself
// return -1 if no such node is found
int ProbGraph::findNodeWithCertainNeighbours(vector<int> nodeVec)
{
    vector<bool> nodeIndVec(size); // vector that states true if the node is included in the set
    vector<int>::iterator vecIt;
    for(vecIt = nodeVec.begin(); vecIt!=nodeVec.end(); ++vecIt)
    {
        nodeIndVec[*vecIt]=true;
    }
    
    // go through the nodes in nodeVec and see if one has only neighbours that are also in nodeVec
    list<int>::iterator listIt;
    bool found;
    for(vecIt = nodeVec.begin(); vecIt!=nodeVec.end(); ++vecIt)
    {
        found = true; // set to false if a neighbour is not in the set
        for(listIt = adjList[*vecIt].begin(); listIt!=adjList[*vecIt].end(); ++listIt)
        {
            if(!nodeIndVec[*listIt])
            {
                found = false;
                break;
            }
        }
        if(found) // all neighbours in the set
        {
            return(*vecIt);
        }
    }
    // nothing found
    return(-1);
}




// generate a list of nodes/sep that represent the elimination sequence in the 
// graph. This list will later be used to generate the junction tree.
list<CliqueSepPair> ProbGraph::generateCliqueSepList()
{
    list<CliqueSepPair> res; // saves the results
    vector<int> clique, sep;
    
    while(!allEliminated())
    {
        // find the next node to eliminate
        int elimNode = findNextElimNode();
        clique = family(elimNode);
        sep = neighbours(elimNode);
        // check if simplical; if not fill in
        if(!isSimplical(elimNode))
        {
            fillIn(elimNode);
        }
        // eliminate the node
        eliminateNode(elimNode);
        
        // now check if others from the same clique can be eliminated
        while((elimNode=findNodeWithCertainNeighbours(sep))!=-1)
        {
            sep = neighbours(elimNode);
            eliminateNode(elimNode); // filling in not necessary by construction
        }
        
        // sort the entries
        sort(clique.begin(), clique.end());
        sort(sep.begin(), sep.end());
        
        res.push_front(CliqueSepPair(clique, sep));
    }
    return(res);
}




//  initialize the potentials in the tree;
void ProbGraph::initializePotentials()
{
    // turn the thetaMatrix into a list
    vector<list<pair<int,double> > > thetaList(size);
    for(int i=0; i!=size; ++i)
    {
        for(int j=i; j!=size; ++j)
        {
            if(thetaMat[i*size+j]!=0)
            {
                thetaList[i].push_back(pair<int,double>(j,thetaMat[i*size+j]));
            }
        }
    }

    // now walk through all nodes in the tree and initialize the potentials, include the thetaMat
    list<JunctionTree>::iterator jtIt;
    vector<TreeNode>::iterator treeIt;
    for(jtIt=jtl.begin(); jtIt!=jtl.end(); ++jtIt)
    {
        for(treeIt=jtIt->tree.begin(); treeIt!=jtIt->tree.end(); ++treeIt)
        {
            // initialize the potential
            vector<int> curVars = treeIt->variables;
            vector<int>::iterator vecIt, vecItInner;
            list<pair<int,double> >::iterator listIt;
            treeIt->pot = Potential(curVars);
        
            // now walk through all variables and see if there is any thetas to include in the potential
            for(vecIt=curVars.begin(); vecIt!=curVars.end(); ++vecIt)
            {
                vecItInner=vecIt;
                for(listIt = thetaList[*vecIt].begin(); listIt!=thetaList[*vecIt].end(); )
                {
                    while(vecItInner!=curVars.end() && listIt->first > *vecItInner){++vecItInner;} // increase sorted variables until at least as large as currently treated variable
                    if(vecItInner!=curVars.end() && listIt->first == *vecItInner)
                    {
                        // add the theta
                        if(*vecIt==*vecItInner) // for single variable
                        {
                            treeIt->pot.multSingleExp(listIt->second, *vecIt);
                        }
                        else // for two variables
                        {
                            treeIt->pot.multDoubleExp(listIt->second, *vecIt, *vecItInner);
                        }
                        ++vecItInner;
                        listIt = thetaList[*vecIt].erase(listIt);
                    }
                    else
                    {
                        ++listIt;
                    }
                }
            }
        }
    }
    return;
}








// initialize the object and generate the adjacency list
ProbGraph::ProbGraph(const int numNodes, const int *graphAdjMat, const double *graphThetaMat, time_t quittingTime) : adjMat(numNodes*numNodes,0), thetaMat(numNodes*numNodes,0), eliminated(numNodes,0), adjList(numNodes,list<int>()), size(numNodes)
{
    // copy the data
    this->quittingTime = quittingTime;
    int vecSize = size*size;
    for(int i=0; i!=vecSize; ++i)
    {
        adjMat[i] = graphAdjMat[i];
        thetaMat[i] = graphThetaMat[i];
    }
    makeJoinTrees();
    initializePotentials();
}


// copy constructor
ProbGraph::ProbGraph(const ProbGraph& x):adjMat(x.adjMat), thetaMat(x.thetaMat), eliminated(x.eliminated), adjList(x.adjList), size(x.size), jtl(x.jtl)
{
    quittingTime = x.quittingTime;
}

// destructor
ProbGraph::~ProbGraph()
{
    // not necessary to call anything
}




// generate a list of junction trees (more than one if the graph is separable)
void ProbGraph::makeJoinTrees()
{
    // first initialize eliminated and the adjacency list
    setAdjList();
    setEliminated();
    
    // now get the list with the cliques and the separators
    list<CliqueSepPair> cliqueSepList = generateCliqueSepList();
    
    // convert these into a list of junction trees
    list<int> activeTreeNodeNums;
    
    while(!cliqueSepList.empty())
    {
        // take the first element the cliqueSepList; start a new tree and make the clique active
        JunctionTree newTree;
        newTree.setQuittingTime(quittingTime);
        int curTreeNodeNum = newTree.addRoot(cliqueSepList.front().first);
        activeTreeNodeNums.push_front(curTreeNodeNum);
        cliqueSepList.pop_front();
        
        // run through the active list one at a time until all elemnts are eliminated
        // for each active clique, find separators that are subsets and add them to the joinTree
        while(!activeTreeNodeNums.empty())
        {
            // get the next currently active treeNode
            curTreeNodeNum = activeTreeNodeNums.front();
            activeTreeNodeNums.pop_front();
            
            // go through the cliqueSepList and find any with a seperator that is a subset of the current node variables
            // if found add it to the tree and the active list; delete from cliqueSepList and go on
            list<CliqueSepPair>::iterator cliqueIt;
            for(cliqueIt = cliqueSepList.begin(); cliqueIt!=cliqueSepList.end(); )
            {
                // if no separator, then root node and not eligible for comparison; ignore
                if((cliqueIt->second != vector<int>(0)) && isSubset(newTree.tree[curTreeNodeNum].variables, cliqueIt->second))
                {
                    int newTreeNodeNum = newTree.addNode(curTreeNodeNum, cliqueIt->first, cliqueIt->second);
                    activeTreeNodeNums.push_back(newTreeNodeNum);
                    cliqueIt = cliqueSepList.erase(cliqueIt);
                }
                else
                {
                    ++cliqueIt;
                }
            }
        }
        jtl.push_front(newTree);
    }
    return;
}


// collect the marginal probabilities from all trees in the object
void ProbGraph::marginalProbs(vector<bool> &varInTree, vector<Distribution>& margDists)
{
    // initialize the variables
    varInTree.clear();
    margDists.clear();
    varInTree.resize(size);
    margDists.resize(size);
    
    // walk through all trees
    list<JunctionTree>::iterator listIt;
    for(listIt = jtl.begin(); listIt!=jtl.end(); ++listIt)
    {
        listIt->getMarginalProbs(varInTree, margDists);
    }
    return;
}


// writes the expectation and covariance matrix. Given pointers have to have sufficient
// space allocated for them
void ProbGraph::getExpectSecMom(double *expectation, double* secMomMat)
{
    // first make a copy of the current junction trees
    list<JunctionTree> jtlCopy=jtl;
    
    // calculate the expection
    vector<bool> varInTree;
    vector<Distribution> margDists;
    marginalProbs(varInTree, margDists);
    // transfer to the double vector
    for(int i=0; i!=size; ++i)
    {
        if(!varInTree[i]) // variable not processed
        {
            throw("Variable not processed");
        }
        expectation[i] = margDists[i][1];
    }

    // now calculate the covariance matrix
    for(int i=0; i!=size; ++i)
    {
        // return the copied trees
        jtl=jtlCopy;
        // initialize varInTree and margDists
        varInTree.clear(); varInTree.resize(size);
        margDists.clear(); margDists.resize(size);
        
        // include the necessary evidence
        enterEvidence(i,1);
        // get marginal probs
        marginalProbs(varInTree, margDists);
        
        // transfer to covMat
        for(int j=0; j!=size; ++j)
        {
            secMomMat[i*size+j] = margDists[j][1]*expectation[i];
        }
    }
    return;
}



void ProbGraph::getExpectSecMomActive(double* secMomMat)
{
    // initialize pairsInTree and secondMomentMatrix
    vector<vector<bool> > pairsInTree; pairsInTree.resize(size);
    vector<vector<double> > secondMoments; secondMoments.resize(size);

    for(int i=0; i!=size; ++i)
    {
        pairsInTree[i].resize(size, false);
        secondMoments[i].resize(size,0);
    }

    // get second moments    // walk through all trees
    list<JunctionTree>::iterator listIt;
    for(listIt = jtl.begin(); listIt!=jtl.end(); ++listIt)
    {
        listIt->getMarginalSecMom(pairsInTree, secondMoments, size, adjMat);
    }
        
    // transfer to secMomMat
    for(int i=0; i!=size; ++i)
    {
        for(int j=0; j!=size; ++j)
        {
            secMomMat[i*size+j] = secondMoments[i][j];
        }
    }
    return;
}



// writes the expectation and covariances for one variable. Given pointers have to have sufficient
// space allocated for them
void ProbGraph::getExpectSecMomSingleVar(int var, double *expectation, double* secMomVec)
{
    // first make a copy of the current junction trees
    list<JunctionTree> jtlCopy=jtl;
    
    // calculate the expection
    vector<bool> varInTree;
    vector<Distribution> margDists;
    marginalProbs(varInTree, margDists);
    // transfer to the double vector
    for(int i=0; i!=size; ++i)
    {
        if(!varInTree[i]) // variable not processed
        {
            throw("Variable not processed");
        }
        expectation[i] = margDists[i][1];
    }

    // return the copied trees
    jtl=jtlCopy;
    // initialize varInTree adn margDists
    varInTree.clear(); varInTree.resize(size);
    margDists.clear(); margDists.resize(size);
        
    // include the necessary evidence
    enterEvidence(var,1);
    // get marginal probs
    marginalProbs(varInTree, margDists);
        
    // transfer to covMat
    for(int j=0; j!=size; ++j)
    {
        secMomVec[j] = margDists[j][1]*expectation[var];
    }
    return;
}


// calculate the normalization constant for the model
double ProbGraph::getNormalizationConstant()
{
    double normConst=1;
    // walk through all trees
    list<JunctionTree>::iterator listIt;
    for(listIt = jtl.begin(); listIt!=jtl.end(); ++listIt)
    {
        normConst*=listIt->getNormalizationConstant();
    }
    return(normConst);

}




// enter evidence into all trees
void ProbGraph::enterEvidence(const int var, const int value)
{
    list<JunctionTree>::iterator listIt;
    for(listIt = jtl.begin(); listIt!=jtl.end(); ++listIt)
    {
        listIt->enterEvidence(var, value);
    }
    return;
}



// print the probGraph object
void ProbGraph::print(ostream &out)
{
    out << "---------------------- ProbGraph Object -----------------------------" << endl;
    // print the adjacency matrix
    out << "Adjacency Matrix" << endl;
    for(int i=0; i!=size; ++i)
    {
        for(int j=0; j!=size; ++j)
        {
            out << adjMat[i*size+j] << " ";
        }
        out << endl;
    }
    
    // print the adjacency list
    out << "Adjacency List" << endl;
    for(int i=0; i!=size; ++i)
    {
        out << "Node " << i << ":";
        list<int>::iterator listIt;
        for(listIt = adjList[i].begin(); listIt!=adjList[i].end(); ++listIt)
        {
            out << *listIt << " ";
        }
        out << endl;
    }
    
    // print the theta matrix
    out << "Theta Matrix" << endl;
    for(int i=0; i!=size; ++i)
    {
        for(int j=0; j!=size; ++j)
        {
            out << thetaMat[i*size+j] << " ";
        }
        out << endl;
    }
    
    // print the vector of eliminated nodes
    out << "Eliminated Nodes" << endl;
    for(int i=0; i!=size; ++i)
    {
        out << eliminated[i] << " ";
    }
    out << endl;
    
    out << "------------------- TREES ----------------------" << endl;
    list<JunctionTree>::iterator listIt;
    for(listIt=jtl.begin(); listIt!=jtl.end(); ++listIt)
    {
        out << "---------- NEW TREE ----------------" << endl;
        listIt->printTree(out);
    }
    
    
    return;
}

