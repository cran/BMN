/*********************************************
***
*** class that encodes a potential of binary variables
*** possible to marginalize over a varible
*** and multiply with another potential
***
*********************************************/

#ifndef _POTENTIAL_
#define _POTENTIAL_
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

typedef vector<double> Distribution;

class Potential
{
    vector<int> variables; // which variables are in the potential
    vector<double> pot; // stores the actual potential
    
    // identify at which position a variable is saved
    // returns -1 if variale is not included at all
    int findLocation(const int var);

    // build a vector with the old and newVariables
    // also build up vectors mappingThis and mappingNew with the information
    // at which position of the joint Variables this and the new Variables are
    void splitVariables(const vector<int>& newVars, vector<int>& jointVars, vector<int>& mappingThis, vector<int>& mappingNew);
    
    // given an integer and a mapping of the bits
    // find a new integer by mapping the old bits to their new positions
    unsigned int mapIndex(const unsigned int x, const vector<int>& mapping);
    
    // marginalize over the variable var and return a new potential
    void marginalize(const int var);

public:
    // initialize the potential
    Potential(const vector<int>& vars);
    Potential();
    // copy constructor
    Potential(const Potential& p);
    
    // destroy the potential
    ~Potential();
    
    // equality and inequality operator
    bool operator==(const Potential& x) const;
    bool operator!=(const Potential& x) const;
    
    // marginalize up to the given set of variables
    void marginalizeToVars(const vector<int> &vars);
    // return the potential as a one dimensional distribution, if only one variable left
    // return value says if operation succeded. Result is written into dist
    bool returnOneDimDist(Distribution& dist) const;
    
    // return the second moment (i.e. P(x1=1,x2=1)) if only 2 variables are left
    bool returnSecondMoment(double& secMom) const;

    // what is the total sum of all values in the potential
    double totalSum();
    
    // normalize so that the sum of all values is 1
    void normalize();
    
    //multiply this potential with another potential
    void multiply(const Potential& factor);
    
    // multiply the potential with a single or double marginal exponential
    void multSingleExp(const double theta, const int var);
    void multDoubleExp(const double theta, const int var1, const int var2);
    // enter evidence into a potential; value can only be either 0 or 1
    void enterEvidence(const int var, const int value);
    
    // print the current potential
    void print(ostream& out);
};


#endif
