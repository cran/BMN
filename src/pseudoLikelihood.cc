#include "pseudoLikelihood.h"
using namespace std;

// identify the active variable from the Theta matrix;
// make sure that every off-diagonal element exists only once
// return value is the number of active variables
// activeR are the active rownumbers, activeC is the columnumbers
void PseudoLikelihood::identifyCurrentlyActive()
{
    numActive = 0;
    for(int i=0; i<numVars; ++i)
    {
        for(int j=i; j<numVars; ++j)
        {
            if(Theta[j*numVars+i]!=0)
            {
                activeR[numActive]=i;
                activeC[numActive]=j;
                ++numActive;
            }
        }
    }
}

// identify variables that will be active in the next step; however, at most add maxVarAdd previously inactive ones
void PseudoLikelihood::identifyFutureActiveWithSort()
{
    numActive=0;
    int index;
    double* gradDiffInactive = new double[numVars*(numVars+1)/2];
    double* gradDiff = new double[numVars*numVars];
    int numInactive=0;
    for(int i=0; i<numVars; ++i)
    {
        for(int j=i; j<numVars; ++j)
        {
            if(Theta[j*numVars+i]!=0)
            {
                activeR[numActive]=i;
                activeC[numActive]=j;
                ++numActive;
            }
            else 
            {
                index = j*numVars+i;
                gradDiff[index] = fabs(gradient[j*numVars+i]) - fabs(rhoMat[j*numVars+i]);
                if(gradDiff[index]>0)
                {
                    gradDiffInactive[numInactive] = gradDiff[index];
                    ++numInactive;
                }
            }
        }
    }
    
    // now sort the difference between gradient and penalty for the inactive variables;
    // only a partial sort is necessary, since we only need a cutoff
    int k = imax2(1, numInactive - maxVarAdd+1); 
    rPsort(gradDiffInactive, numInactive, k);
    double cutoff = gradDiffInactive[k-1];
    // now go through all variables again and add the one above the cutoff that are also currently theta=0
    for(int i=0; i<numVars; ++i)
    {
        for(int j=i; j<numVars; ++j)
        {
            if(Theta[j*numVars+i]==0 && gradDiff[j*numVars+i]>=cutoff)
            {
                activeR[numActive]=i;
                activeC[numActive]=j;
                ++numActive;
            }
        }
    }
    delete [] gradDiffInactive;
    delete [] gradDiff;
}

// copies the upper triangular part of x onto y using the active data stored currently in the object
void PseudoLikelihood::copyThetaActive(const double* x, double* y)
{
    for(int i=0; i<numActive; ++i)
    {
        y[activeC[i]*numVars+activeR[i]]=x[activeC[i]*numVars+activeR[i]];
    }
}

// finds the maximal difference in the upper-diagonal matrix between matrix x and y
double PseudoLikelihood::deltaThetaActive(const double* x, const double* y)
{
    double delta = 0;
    for(int i=0; i<numActive; ++i)
    {
        delta = fmax2(delta, fabs(y[activeC[i]*numVars+activeR[i]]-x[activeC[i]*numVars+activeR[i]]));
    }
    return(delta);
}


// calculate the model probabilities pHat
// will be saved under the pointer pHat, which has to have enough space
// all matrices are saved internally as columnwise vectors
void PseudoLikelihood::calcModelEtas(double *Theta)
{
    // first initialize the pHat matrix
    for(int i=0; i<numVars; ++i)
    {
        // initialize to the diagonal elements
        double foo = Theta[i*numVars+i];
        int index = i*numObs;
        for(int j=0; j<numObs; ++j)
        {
            etaHat[index]=foo;
            ++index;
        }
    }
    // go through all the active variables
    for(int i=0; i<numActive; ++i)
    {
        // use the variable if it is not on the diagonal (already processed)
        if(activeR[i]!=activeC[i])
        {
            int inc = 1;
            // use DAXPY which replaces y by alpha*x + y
            F77_CALL(daxpy)(&numObs, &Theta[activeC[i]*numVars+activeR[i]], X+activeC[i]*numObs, &inc, etaHat+activeR[i]*numObs, &inc);
            F77_CALL(daxpy)(&numObs, &Theta[activeC[i]*numVars+activeR[i]], X+activeR[i]*numObs, &inc, etaHat+activeC[i]*numObs, &inc);
        }
    }
}

// calculate the model probabilities from the model etas
void PseudoLikelihood::calcModelProbs(double* pHat)
{
    int index = 0;
    for(int i=0; i<numVars; ++i)
    {
        for(int j=0; j<numObs; ++j)
        {
            pHat[index] = 1/(1+exp(-etaHat[index]));
            ++index;
        }
    }
}

void PseudoLikelihood::calculatePHatTXAll(double* pHatTX, const double* pHat)
{
    // go over the active variables
    for(int i=0; i<numVars; ++i)
    {
        for(int j=0; j< numVars; ++j)
        {
            if(i!=j)
            {
                int incr = 1;
                // ddot is the BLAS inner product function with double precision
                pHatTX[j*numVars+i] = F77_CALL(ddot)(&numObs, pHat+i*numObs,&incr, X+j*numObs, &incr);
            }
            else // computation only on the diagonal; take the sum over the pHat
            {
                int index = i*numVars + i;
                int indexInner = i*numObs;
                // ddot is the BLAS inner product function with double precision
                pHatTX[index] = 0;
                for(int k=0; k< numObs; ++k)
                {
                    pHatTX[index] += pHat[indexInner];
                    ++indexInner;
                }
            }
        }
    }
}


// calculate only the active pHat^T %*% X; however on the diagonal, do not use X[,i], but a vector consisting of 1s
void PseudoLikelihood::calculatePHatTXActive(double* pHatTX, const double* pHat)
{
    // go over the active variables
    for(int i=0; i<numActive; ++i)
    {
        if(activeR[i]!=activeC[i]) // on the off diagonal also the symmetric elements have to be computed
        {
            int incr = 1;
            // ddot is the BLAS inner product function with double precision
            pHatTX[activeC[i]*numVars+activeR[i]] = F77_CALL(ddot)(&numObs, pHat+activeR[i]*numObs,&incr, X+activeC[i]*numObs, &incr);
            pHatTX[activeR[i]*numVars+activeC[i]] = F77_CALL(ddot)(&numObs, pHat+activeC[i]*numObs,&incr, X+activeR[i]*numObs, &incr);
        }
        else // computation only on the diagonal; take the sum over the pHat
        {
            int index = activeC[i]*numVars + activeR[i];
            int indexInner = activeR[i]*numObs;
            // ddot is the BLAS inner product function with double precision
            pHatTX[index] = 0;
            for(int k=0; k< numObs; ++k)
            {
                pHatTX[index] += pHat[indexInner];
                ++indexInner;
            }
        }
    }
}



void PseudoLikelihood::calcGradientAll()
{
    // first calculate the etas
    double* pHat = new double[numVars*numObs];
    double* pHatTX = new double[numVars*numVars];
    calcModelEtas(Theta);
    // get the probabilities from the etas
    calcModelProbs(pHat);
    // now get pHatTX for all variables
    calculatePHatTXAll(pHatTX, pHat);
    
    for(int i=0; i<numVars; ++i)
    {
        for(int j=i; j<numVars; ++j)
        {
            // check if we are on the diagonal
            if(i!=j)
            {
                int index1=j*numVars+i;
                int index2=i*numVars+j;
                gradient[index1] = 2*XTX[index1] - pHatTX[index1] - pHatTX[index2] - Delta[index1];
            }
            else // now we are on the diagonal
            {
                int index = i*numVars + i;
                gradient[index] = XTX[index] - pHatTX[index] - Delta[index];
            }
        }
    }
    delete [] pHat;
    delete [] pHatTX;
}

// given a matrix of thetas and the active variables, calculate the gradients of the active variables to use it 
// in the inner loop
void PseudoLikelihood::calcGradientActive()
{
    // first calculate the etas
    double* pHat = new double[numVars*numObs];
    double* pHatTX = new double[numVars*numVars];
    // get the probabilities from the etas
    calcModelProbs(pHat);
    // now get pHatTX
    calculatePHatTXActive(pHatTX, pHat);
    
//    Rprintf("pHatTX\n");
//    printMatrixDouble(pHatTX, numVars, numVars);
//    Rprintf("XTX\n");
//    printMatrixDouble(XTX, numVars, numVars);

    
    // now use this data to calculate the gradient, but only over the active variables
    for(int i=0; i<numActive; ++i)
    {
        // check if we are on the diagonal
        if(activeR[i]!=activeC[i])
        {
            int index1=activeC[i]*numVars+activeR[i];
            int index2=activeR[i]*numVars+activeC[i];
            gradient[index1] = 2*XTX[index1] - pHatTX[index1] - pHatTX[index2] - Delta[index1];
        }
        else // now we are on the diagonal
        {
            int index = activeC[i]*numVars + activeR[i];
            gradient[index] = XTX[index] - pHatTX[index] - Delta[index];
        }
    }
    delete [] pHat;
    delete [] pHatTX;
}


double PseudoLikelihood::calcPenLogLik(double* Theta)
{
    double result=0;
    int index;
    // first add up the normalization constants
    for(int i=0; i< numObs; ++i)
    {
        for(int j=0; j < numVars; ++j)
        {
//            result -= logspace_add(etaHat[j*numObs+i],0);
            result -= log1p(exp(etaHat[j*numObs+i]));
        }
    }
    
    // now add the linear part and subtract the penalty, also adjust for Delta
    for(int i=0; i<numActive; ++i)
    {
        if(activeR[i]==activeC[i]) // on the diagonal
        {
            index = activeR[i]*numVars+activeR[i];
            result += (XTX[index] -Delta[index])* Theta[index] - rhoMat[index]*fabs(Theta[index]);
        }
        else // on the off-diagonal need corresponding elements on both sides, which are equal, so times 2
        {
            index = activeC[i]*numVars+activeR[i];
            result += (2*XTX[index]-Delta[index])*Theta[index] - rhoMat[index]*fabs(Theta[index]);
        }
    }
    return(result);
}

// function that find the minimum for an L1-penalized quadratic function
// a>0 assumed but not checked by the function for performance reason 
double PseudoLikelihood::minQuadL1(double a, double b, double rho)
{
    if(b+rho<0) // x has to be positive
    {
        return(-(b+rho)/a);
    }
    else if(b-rho <0) // x is 0
    {
        return(0);
    }
    else
    {
        return((rho-b)/a);
    }
}



// adjust the gradient for the penalty
void PseudoLikelihood::adjustGradientForPenalty(double* adjGradient, const double* Theta)
{
    int index;
    for(int i=0; i<numActive; ++i)
    {
        index = activeC[i]*numVars + activeR[i];
        if(Theta[index]!=0)
        {
            adjGradient[index] = gradient[index] - rhoMat[index]*sign(Theta[index]);
        }
        else
        {
            adjGradient[index] = gradient[index] - sign(gradient[index])*rhoMat[index];
        }
    }
}


// calculate the descent direction times the gradient and return the value
double PseudoLikelihood::calcGradientTimesAscentDirec(const double* gradient, const double *oldTheta, const double *newTheta)
{
    int index;
    double result=0;
    for(int i=0; i<numActive; ++i)
    {
        index = activeC[i]*numVars+activeR[i];
        result += (newTheta[index]-oldTheta[index])*gradient[index];
    }
    return(result);
}

// shortens the difference between oldTheta and newTheta by a factor beta and writes it back into
// newTheta
void PseudoLikelihood::calcBacktrackTheta(double *newTheta, const double* oldTheta, const double beta)
{
    int index;
    for(int i=0; i<numActive; ++i)
    {
        index = activeC[i]*numVars + activeR[i];
        newTheta[index] = oldTheta[index] + beta *(newTheta[index]-oldTheta[index]);
    }
}


// performs a line search, writes the used Theta value into newTheta, overwriting the guess from before;
// as a return value, the new logLikelihood is stored in penLogLik
double PseudoLikelihood::lineSearch(const double* oldTheta, double* newTheta, const double alpha, const double beta)
{
    double* adjGradient = new double[numVars*numVars];
    calcModelEtas(newTheta);
    double newPenLogLik=calcPenLogLik(newTheta);
    // adjust the gradient for the penalty
    adjustGradientForPenalty(adjGradient, oldTheta);
    // get the descent direction
    double lowerBoundGrad = calcGradientTimesAscentDirec(adjGradient, oldTheta, newTheta);
    
    // now do the backtracking line search
    double t=1;
//    Rprintf("Lowerbound %f\n",lowerBoundGrad);
//    Rprintf("PenLogLik %f\n",(*penLogLik));
//    Rprintf("NewPenLogLik %f\n",newPenLogLik);
    while(newPenLogLik < penLogLik + alpha * t * lowerBoundGrad)
    {
        t = t*beta;
        // make sure to stop if t gets too small; use the last Theta and eta computed (which are almost the same as the starting one
        if(t < 1e-5)
        {
            break;
        }
        calcBacktrackTheta(newTheta, oldTheta, beta);
        calcModelEtas(newTheta);
        newPenLogLik = calcPenLogLik(newTheta);
    }
    // save the current penalize logLikelihood
    penLogLik = newPenLogLik;
    delete [] adjGradient;
    return(t);
}

void PseudoLikelihood::printMatrixActive(double* Mat)
{
    cout << "Printing matrix:" << endl;
    for(int i=0; i<numActive; ++i)
    {
        cout << activeR[i] << "," << activeC[i] << ":" << Mat[activeC[i]*numVars+activeR[i]] << endl;
    }
}


// runs the inner loop of the pseudo algorithm which optimizes the pseudo-likelihood over a given set of variables
// return value indicates if it was successful; 0 for success; 1 for failure
void PseudoLikelihood::pseudoInnerLoop()
{
    double maxChange; // variable tracks the maximal change for all Theta values; if under thr, Theta is considered to be the solution
    double* oldThetaInner = new double[numVars*numVars];
    double* savedTheta = new double[numVars*numVars];
    double savedPenLogLik;
    bool performCheck;
    
    // necessary calculations for the first run of the inner loop
    calcModelEtas(Theta);
    copyThetaActive(Theta, savedTheta);
    penLogLik = calcPenLogLik(Theta);
    savedPenLogLik = penLogLik;

    // run the inner loop
    for(innerIter=0; innerIter < maxInnerIter; ++innerIter)
    {
        printMatrixActive(Theta);
        printMatrixActive(gradient);
        // check if in this iteration a line search is performed
        if(innerIter % CHECKFREQUENCY==0)
        {
            performCheck = true;
        }
        else
        {
            performCheck = false;
        }
        // initialize maxChange
        maxChange = 0;
        // calculate the gradient for the currently active variables
        calcGradientActive();
        // if line search performed in this iteration, save theta
        if(performLineSearch)
        {
            copyThetaActive(Theta, oldThetaInner);
        }

        // go through all currently active variables
        for(int i=0; i< numActive; ++i)
        {
            int index = activeR[i] + activeC[i] * numVars;
            double theta0 = Theta[index],a;
            // calculate the factor of (theta-theta0)
            if(activeR[i]==activeC[i])
            {
                a = numObs/4/(stepSize); // on the diagonal vector of 1's instead of column of X used (same as in pHatTX)
            }
            else // not on the diagonal
            {
                a = (XTX[activeR[i]*numVars+activeR[i]] + XTX[activeC[i]*numVars+activeC[i]])/4/(stepSize);
            }
            
            double b = -gradient[index]-theta0*a;
            
            cout << activeR[i] << "," << activeC[i] << ":" << a << "," << b << endl;
            
            Theta[index] = minQuadL1(a,b,rhoMat[index]);
            maxChange = fmax2(maxChange, fabs(Theta[index]-theta0));
        }
        cout << "MaxChange" << maxChange << endl;
        cout << "StepSize" << stepSize << endl;
        cout << "Iteration: " << innerIter << endl;
//        maxChange *=t;
//        Rprintf("Step length: %f\n", t);
//        Rprintf("Max change: %f\n", maxChange);
        if(maxChange <= thr)
        {
            delete [] oldThetaInner;
            delete [] savedTheta;
            return;
        }
        if(performLineSearch)
        {
            // line search also calculated model etas
            lineSearch(oldThetaInner, Theta, 0.1, 0.5);
        }
        else // no line searches
        {
            if(performCheck)
            {
                // check if the penalized Log likelihood has improved;
                calcModelEtas(Theta);
                penLogLik = calcPenLogLik(Theta);
//                Rprintf("InnerIter: %d, SavedPenLogLik %f, PenLogLik %f, stepSize %f\n", innerIter, savedPenLogLik, penLogLik, *stepSize);
                if(penLogLik >= savedPenLogLik - thr) // did improve
                {
                    copyThetaActive(Theta, savedTheta);
                    savedPenLogLik = penLogLik;
                }
                else // did not improve, revert to old setting
                {
                    // half the number of line steps and restore the old value for theta
                    innerIter -= CHECKFREQUENCY;
                    copyThetaActive(savedTheta, Theta);
                    stepSize /= 2;
                    if(stepSize < 0.001) // stepSize too small; switch to a line search
                    {
                        performLineSearch=true;
                        stepSize=1;
                    }
 //                   Rprintf("StepSize reduced to %f\n", *stepSize);
                }
            }
            calcModelEtas(Theta);
        }
//        Rprintf("Loglik %f\n", penLogLik);
    }
    delete [] oldThetaInner;
    delete [] savedTheta;

//    Rprintf("InnerIter: %d\n", innerIter);
//    Rprintf("Active Inner: %d\n", numActive);
}




void PseudoLikelihood::pseudoOuterLoop()
{
    // initialize variables
    double maxThetaChange=0;
    identifyCurrentlyActive();
    double* oldThetaOuter = new double[numVars*numVars];
    
    // now run the outer loop
    for(outerIter=0; outerIter<maxOuterIter; ++outerIter)
    {
        // in the outer loop, a set of active variables is set
        // then the inner loop is run until convergence
        // then a new active set is defined
        // set the maximum number of inner loops to run to MAX_INNER_ITER
        
        // copy the current theta vector to oldTheta)
        copyThetaActive(Theta, oldThetaOuter);
        
        pseudoInnerLoop();
        if(innerIter==maxInnerIter) // check if the number of iterations in the inner loop was at maximum
        {// error occured
            outerIter = maxOuterIter; 
            delete [] oldThetaOuter;
            return; // will signal an error
        }
        
        // find out how much Theta changed from the last iteration
        maxThetaChange = deltaThetaActive(Theta, oldThetaOuter);
//        Rprintf("OuterMaxThetaChange %f\n", maxThetaChange);
        if(maxThetaChange < thr && outerIter > 0) // only break if at least one iteration has been completed
        {
            success = true;
            delete [] oldThetaOuter;
            return;
        }
        
        // find the new active variables
        calcGradientAll();
//        numActive=identifyFutureActive(activeR, activeC, Theta, gradient, rhoMat, numVars);
        identifyFutureActiveWithSort();
    }
}





// sets all the variables and starts the calculations
PseudoLikelihood::PseudoLikelihood(const double* X, const double* XTX, const double* rhoMat, const double* Delta, const double* ThetaStart, const double thr, const int maxOuterIter, const int maxInnerIter, const int maxVarAdd, const int numObs, const int numVars, double stepSize, bool performLineSearch) : X(X), XTX(XTX), rhoMat(rhoMat), Delta(Delta), thr(thr), maxOuterIter(maxOuterIter), maxInnerIter(maxInnerIter), maxVarAdd(maxVarAdd), numObs(numObs), numVars(numVars)
{
    // initialize all the variables
    this->stepSize = stepSize;
    this->performLineSearch = performLineSearch;
    this->success = false;

    // arrays that store which variables are active
    numActive=0;
    activeC = new int[numVars*numVars];
    activeR = new int[numVars*numVars];
    // prepare the matrix that stores the current gradient
    gradient = new double[numVars*numVars];
    etaHat = new double[numVars*numObs];

    // Theta has to be copied over from ThetaStart
    int incr=1;
    int ThetaLen = numVars*numVars;
    this->Theta = new double[numVars*numVars];
    F77_CALL(dcopy)(&ThetaLen, ThetaStart, &incr, Theta, &incr);
    pseudoOuterLoop();
}


PseudoLikelihood::~PseudoLikelihood()
{
    delete[](Theta);
    delete[](activeC);
    delete[](activeR);
    delete[](gradient);
    delete[](etaHat);
}


double* PseudoLikelihood::returnTheta()
{
    return(Theta);
}

bool PseudoLikelihood::returnSuccess()
{
    return(success);
}

// take the upper triangular matrix in Theta and copy it onto the 
// lower triangular matrix
void symmetrizeTheta(double* Theta, const int numVars)
{
    for(int i=0; i<numVars; ++i)
    {
        for(int j=i+1; j<numVars; ++j)
        {
            Theta[i*numVars+j]=Theta[j*numVars+i];
        }
    }
}


extern "C" {
// write the interface function to R
// Parameters:  X is an n by p matrix
//              XTX is the inner product of X with itself
//              rho is a p by p matrix of the penalty parameters
//              Delta is the adjustment matrix for the gradients, p by p matrix
//              ThetaStart is an initialization value for Theta, p by p matrix
//              thr is the maximum change of Theta allowed at convergence
//              maxIter is the maximum number of iterations allowed in the algorithms outer loop before it stops
SEXP pseudoLikelihood(SEXP R_X, SEXP R_XTX, SEXP R_rhoMat, SEXP R_Delta, SEXP R_ThetaStart, SEXP R_thr, SEXP R_maxOuterIter, SEXP R_maxInnerIter, SEXP R_maxVarAdd, SEXP R_stepSize, SEXP R_performLineSearch)
{
    // protect the input items
    PROTECT(R_X);
    PROTECT(R_XTX);
    PROTECT(R_rhoMat);
    PROTECT(R_Delta);
    PROTECT(R_ThetaStart);
    PROTECT(R_thr);
    PROTECT(R_maxOuterIter);
    PROTECT(R_maxInnerIter);
    PROTECT(R_maxVarAdd);
    PROTECT(R_stepSize);
    PROTECT(R_performLineSearch);
    
    // get a few necessary parameters (consistency check is done in the R wrapper function)
    int numObs = INTEGER(GET_DIM(R_X))[0];
    int numVars = INTEGER(GET_DIM(R_X))[1];
    
    // get space for the output Theta matrix
    SEXP R_Theta, R_success;
    PROTECT(R_Theta = allocMatrix(REALSXP, numVars, numVars));
    PROTECT(R_success = allocVector(LGLSXP, 1));
    
    // now run the actual pseudolikelihood code
    PseudoLikelihood pseudoObject(REAL(R_X), REAL(R_XTX), REAL(R_rhoMat), REAL(R_Delta), REAL(R_ThetaStart), REAL(R_thr)[0], INTEGER(R_maxOuterIter)[0], INTEGER(R_maxInnerIter)[0], INTEGER(R_maxVarAdd)[0], numObs, numVars, REAL(R_stepSize)[0], LOGICAL(R_performLineSearch)[0]);
    // only upper triangular used; make matrix symmetric
    double* Theta = pseudoObject.returnTheta();
    symmetrizeTheta(Theta, numVars);

    // now copy the data onto the R_object
    int ThetaLen = numVars*numVars;
    int incr=1;
    F77_CALL(dcopy)(&ThetaLen, Theta, &incr, REAL(R_Theta), &incr);

    LOGICAL(R_success)[0]=pseudoObject.returnSuccess();
    
    // generate a list for the return values
    SEXP retList, dimnames;
    PROTECT(retList = allocVector(VECSXP,2));
//    SET_VECTOR_ELT(retList,0,R_Theta);
    SET_VECTOR_ELT(retList,0,R_Theta);
    SET_VECTOR_ELT(retList,1,R_success);
    PROTECT(dimnames = allocVector(STRSXP,2));
    SET_STRING_ELT(dimnames, 0, mkChar("Theta"));
    SET_STRING_ELT(dimnames, 1, mkChar("success"));
    setAttrib(retList,R_NamesSymbol, dimnames); 
    
    // unprotect them again
    UNPROTECT(15);
    return(retList);
}
}

