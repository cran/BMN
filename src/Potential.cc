#include "Potential.h"
#include <algorithm>
#include <set>

int Potential::findLocation(const int var)
{
    int location = -1;
    for(int i=0; i< (int) variables.size(); ++i)
    {
        if(variables[i]==var)
        {
            location = i;
            break;
        }
    }
    return(location);
}





// build a vector with the old and newVariables
// also build up vectors mappingThis and mappingNew with the information
// at which position of the joint Variables this and the new Variables are
void Potential::splitVariables(const vector<int>& newVars, vector<int>& jointVars, vector<int>& mappingThis, vector<int>& mappingNew)
{
    set<int> jointSet, thisSet, newSet;
    
    // generate all these sets
    for(unsigned int i=0; i!=variables.size(); ++i)
    {
        jointSet.insert(variables[i]);
        thisSet.insert(variables[i]);
    }
    for(unsigned int i=0; i!=newVars.size(); ++i)
    {
        jointSet.insert(newVars[i]);
        newSet.insert(newVars[i]);
    }

    // transfer jointset into a vector
    set<int>::iterator setIt;
    jointVars.clear();
    for(setIt = jointSet.begin(); setIt!=jointSet.end(); ++setIt)
    {
        jointVars.push_back(*setIt);
    }
    
    // find the mapping
    mappingThis.clear();
    mappingNew.clear();
    for(unsigned int jointPos=0; jointPos!=jointVars.size(); ++jointPos)
    {
        if(thisSet.count(jointVars[jointPos]))
        {
            mappingThis.push_back(jointPos);
        }
        if(newSet.count(jointVars[jointPos]))
        {
            mappingNew.push_back(jointPos);
        }
    }
    return;
}





// given a integer and a mapping of the bits
// find a new integer by mapping the old bits to their new positions
// the value at mapping[i] describes that the mapping[i]-th bit of x should be at
// position i of y
unsigned int Potential::mapIndex(const unsigned int x, const vector<int>& mapping)
{
    int len = mapping.size();
    unsigned int y=0;
    
    for(int i=0; i!=len; ++i)
    {
        if( x & (1 << mapping[i])) // is the ith bit set, set the mapped bit in y
        {
            y |= (1 << i);
        }
    }
    return(y);
}




// standard constructor
Potential::Potential(const vector<int>& vars) : variables(vars), pot(1 << vars.size(),1)
{
    sort(variables.begin(), variables.end());
    if(vars.size()>20)
    {
        throw("Too many variables in one clique. Must be 20 or less.");
    }
}

Potential::Potential()
{
}


// copy constructor
Potential::Potential(const Potential& p) : variables(p.variables), pot(p.pot)
{
}


    // destroy the potential
Potential::~Potential()
{
}

// equality and inequality operator
bool Potential::operator==(const Potential& x) const
{
    return((this->variables==x.variables) && (this->pot==x.pot));
}


bool Potential::operator!=(const Potential& x) const
{
    return(!(*this==x));
}




// marginalize up to the given set of variables
void Potential::marginalizeToVars(const vector<int> &vars)
{
    vector<int> potVars(variables); // need to be copied as they are changed by marginalizing
    vector<int>::const_iterator vecIt, vecItVars;
    bool isInVars=false;
    for(vecIt = potVars.begin(); vecIt!=potVars.end(); ++vecIt)
    {
        // see if the variable is in vars
        isInVars=false;
        for(vecItVars=vars.begin(); vecItVars!=vars.end(); ++vecItVars)
        {
            if(*vecIt==*vecItVars){isInVars=true;}
        }
        
        // if no, marginalize
        if(!isInVars){marginalize(*vecIt);}
    }
    return;
}


// return the potential as a one dimensional distribution, if only one variable left
// return value says if operation succeded. Result is written into dist
bool Potential::returnOneDimDist(Distribution& dist) const
{
    if(variables.size()!=1)
    {
        return(false);
    }
    dist = pot;
    return(true);
}

// return the second moment of the potential when only two variables are left
bool Potential::returnSecondMoment(double& secMom) const
{
    if(variables.size()!=2)
    {
        return(false);
    }
    secMom = pot[3];
    return(true);
}



// marginalize over the variable var and return a new potential
void Potential::marginalize(const int var)
{
    int location = findLocation(var);
    
    // initialize the new variables
    vector<int> newVars = variables;
    newVars.erase(newVars.begin() + location);
    int newLength = 1 << newVars.size();
    vector<double> newPot(newLength);
    
    // if variable exists, marginalize it out
    if(location !=-1)
    {
        // walk through everything before and after the location
        // need 2 loops for this
        unsigned int innerLength = 1 << location;
        unsigned int outerLength = 1 << (variables.size() - 1 - location);
        unsigned int curIndex, curIndexNew, outerIndex;
        unsigned int maskLocation = 1 << location;
        for(unsigned int outer = 0; outer!=outerLength; ++outer)
        {
            outerIndex = outer << (location + 1);
            for(unsigned int inner = 0; inner!=innerLength; ++inner)
            {
                // sum over the variable at location
                curIndexNew = (outerIndex >> 1) | inner;
                curIndex = outerIndex | inner;
                newPot[curIndexNew] = pot[curIndex] + pot[curIndex | maskLocation];
            }
        }
        // copy all the new elements; free memory where needed
        pot = newPot;
        variables = newVars;
    }
    return;
}



// find the sum of all values in the potential
double Potential::totalSum()
{
    double sum=0;
    for(unsigned int i=0; i!=pot.size(); ++i)
    {
        sum+=pot[i];
    }
    return(sum);
}


// normalize so that the sum of all values is 1
void Potential::normalize()
{
    double sum = totalSum();
    for(unsigned int i=0; i!=pot.size(); ++i)
    {
        pot[i]/=sum;
    }
    return;
}


//


//multiply this potential with another potential
void Potential::multiply(const Potential& factor)
{
    // join the variables and return vector with the positions
    // in the joint set
    vector<int> jointVars, mappingThis, mappingFactor;
    splitVariables(factor.variables, jointVars, mappingThis, mappingFactor);
    
    // run through all indices combinations in the joint variables
    // map it onto the indices of the variables in this and in factor
    unsigned int jointLength = 1 << jointVars.size();
    unsigned int thisIndex, factorIndex;
    // reserve new space
    vector<double> jointPot(jointLength);
    for(unsigned int jointIndex = 0; jointIndex!=jointLength; ++jointIndex)
    {
        thisIndex = mapIndex(jointIndex, mappingThis);
        factorIndex = mapIndex(jointIndex, mappingFactor);
        jointPot[jointIndex] = pot[thisIndex] * factor.pot[factorIndex];
    }
    // free the old space and set the new
    pot=jointPot;
    variables = jointVars;
    return;
}


// multiply the potential with a single marginal exponential
void Potential::multSingleExp(const double theta, const int var)
{
    // the multiplier for var = 1
    double multiple = exp(theta);

    int location = findLocation(var);
    
    // if variable exists, multiply
    if(location !=-1)
    {
        // walk through everything before and after the location
        // need 2 loops for this
        unsigned int innerLength = 1 << location;
        unsigned int outerLength = 1 << (variables.size() - 1 - location);
        unsigned int curIndex, outerIndex;
        unsigned int maskLocation = 1 << location;
        for(unsigned int outer = 0; outer!=outerLength; ++outer)
        {
            outerIndex = outer << (location + 1);
            for(unsigned int inner = 0; inner!=innerLength; ++inner)
            {
                // multiply
                curIndex = outerIndex | inner | maskLocation;
                pot[curIndex] *= multiple;
            }
        }
    }
    return;
}



// multiply the potential with a single or double marginal exponential
void Potential::multDoubleExp(const double theta, const int var1, const int var2)
{
    // find the location; if the variable is not included, do nothing
    int location1 = findLocation(var1);
    int location2 = findLocation(var2);
    if((location1==-1) || (location2==-1))
    {
        return;
    }
    
    // if both variables the same, singleExp
    if(location1==location2)
    {
        multSingleExp(theta, var1);
        return;
    }
    
    // sort location1 and 2 by size
    if(location1 > location2)
    {
        int foo = location2;
        location2 = location1;
        location1 = foo;
    }
    // the multiplier for var = 1
    double multiple = exp(theta);
    
    // walk through everything before and after the location
    // need 3 loops for this
    unsigned int innerLength = 1 << location1;
    unsigned int middleLength = 1 << (location2-location1-1);
    unsigned int outerLength = 1 << (variables.size() - 1 - location2);
    unsigned int curIndex, outerIndex, middleIndex;
    unsigned int maskLocation = (1 << location1) | (1<<location2);
    
    for(unsigned int outer = 0; outer!=outerLength; ++outer)
    {
        outerIndex = outer << (location2 + 1);
        for(unsigned int middle = 0; middle!=middleLength; ++middle)
        {
            middleIndex = middle << (location1 + 1);
            for(unsigned int inner = 0; inner!=innerLength; ++inner)
            {
                // multiply
                curIndex = outerIndex | middleIndex | inner | maskLocation;
                pot[curIndex] *= multiple;
            }
        }
    }
    return;
}



void Potential::enterEvidence(const int var, const int value)
{
    if(value >1 || value < 0) // value has to be either 0 or 1
    {
        return;
    }
    
    int location = findLocation(var);
    
    // if variable exists, multiply
    if(location !=-1)
    {
        // walk through everything before and after the location
        // need 2 loops for this
        unsigned int innerLength = 1 << location;
        unsigned int outerLength = 1 << (variables.size() - 1 - location);
        unsigned int curIndex, outerIndex;
        unsigned int maskLocation;
        if(value == 0) // set all 1's to 0
        {
            maskLocation = 1 << location  ;
        }
        else // set all 0's to 0
        {
            maskLocation = 0 ;
        }
        for(unsigned int outer = 0; outer!=outerLength; ++outer)
        {
            outerIndex = outer << (location + 1);
            for(unsigned int inner = 0; inner!=innerLength; ++inner)
            {
                // multiply
                curIndex = outerIndex | inner | maskLocation;
                pot[curIndex] = 0;
            }
        }
    }
    return;
}





// print the current potential
void Potential::print(ostream& out)
{
    // first print the variables that were used
    out << "Variables:" << endl;
    for(unsigned int i=0; i!= variables.size(); ++i)
    {
        out << variables[i] << " ";
    }
    out << endl;
    
    // now print the actual potential; 4 on every line
    for(unsigned int i=0; i!=pot.size(); ++i)
    {
        if(i % 4 ==0)
        {
            out.setf(ios::oct,ios::basefield);
            out << "Oct: " << i << " :";
            out.setf(ios::dec, ios::basefield);
        }
        out << pot[i] << "  ";
        if(i % 4 ==3)
        {
            out << endl;
        }
    }
    return;
}


