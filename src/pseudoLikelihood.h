#include <iostream>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <math.h>
#include <Rmath.h>

const int CHECKFREQUENCY = 50;

class PseudoLikelihood
{
    double *Theta;
    const double *X;
    const double *XTX;
    const double *rhoMat;
    const double *Delta;
    const double *ThetaStart;
    const double thr;
    const int maxOuterIter;
    const int maxInnerIter;
    const int maxVarAdd;
    const int numObs;
    const int numVars;
    double stepSize;
    bool performLineSearch;
    bool success;
    // iteration value for the inner and outer loops
    int innerIter;
    int outerIter;

    // arrays that store which variables are active
    int numActive;
    int* activeC; 
    int* activeR;
    // prepare the matrix that stores the current gradient
    double* gradient;
    double *etaHat;
    double penLogLik;

    // helper functions
    void identifyCurrentlyActive();
    void identifyFutureActiveWithSort();

    void copyThetaActive(const double* x, double* y);
    double deltaThetaActive(const double* x, const double* y);
    void calcGradientAll();
    void calcGradientActive();

    void calcModelEtas(double* Theta);
    void calcModelProbs(double* pHat);
    void calculatePHatTXAll(double* pHatTX, const double* pHat);
    void calculatePHatTXActive(double* pHatTX, const double* pHat);
    double calcPenLogLik(double* Theta);
    double minQuadL1(double a, double b, double rho);

    // functions necessary for the line search
    void adjustGradientForPenalty(double* adjGradient, const double* Theta);
    double calcGradientTimesAscentDirec(const double* gradient, const double *oldTheta, const double *newTheta);
    void calcBacktrackTheta(double *newTheta, const double* oldTheta, const double beta);
    double lineSearch(const double* oldTheta, double* newTheta, const double alpha, const double beta);

    void printMatrixActive(double* Mat);


    // runs the inner loop
    void pseudoInnerLoop();
    // runs the outer loop
    void pseudoOuterLoop();


public:

    PseudoLikelihood(const double* X, const double* XTX, const double* rhoMat, const double* Delta, const double* ThetaStart, const double thr, const int maxOuterIter, const int maxInnerIter, const int maxVarAdd, const int numObs, const int numVars, double stepSize, bool performLineSearch);
    
    ~PseudoLikelihood();
    
    double* returnTheta();
    bool returnSuccess();
};

