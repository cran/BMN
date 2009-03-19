/***************************************************************
***
*** write a wrapper function for probgraph to R
***
***
***************************************************************/


#include "ProbGraph.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <time.h>

extern "C"
{

SEXP runJTAlgSecMomMat(SEXP adjMat, SEXP thetaMat, SEXP maxRuntime)
{
    // protect the items
    PROTECT(adjMat);
    PROTECT(thetaMat);
    PROTECT(maxRuntime);
    
    int size = INTEGER(GET_DIM(adjMat))[0];
    time_t quittingTime = time(NULL) + INTEGER(maxRuntime)[0];
    
    // get space for the expectation and covariance matrix
    SEXP covMat, expec;
    PROTECT(covMat = allocMatrix(REALSXP, size, size));
    PROTECT(expec = allocVector(REALSXP, size));

    double *covMatPtr = REAL(covMat);
    double *expecPtr = REAL(expec);
    
    // get the probgraph object
    try
    {
        ProbGraph graph(size, INTEGER(adjMat), REAL(thetaMat), quittingTime);
        graph.getExpectSecMom(expecPtr, covMatPtr);
    }
    catch(char const* errorMsg)
    {
        error(errorMsg);
    }
    
    // generate a list for the return values
    SEXP retList, dimnames;
    PROTECT(retList = allocVector(VECSXP,2));
    SET_VECTOR_ELT(retList,0,expec);
    SET_VECTOR_ELT(retList,1,covMat);
    PROTECT(dimnames = allocVector(STRSXP,2));
    SET_STRING_ELT(dimnames, 0, mkChar("Expectation"));
    SET_STRING_ELT(dimnames, 1, mkChar("SecondMomentMatrix"));
    setAttrib(retList,R_NamesSymbol, dimnames); 
    
    // unprotect them again
    UNPROTECT(7);
    return(retList);
}





// get the second moments but only for the edges in the graph
SEXP runJTAlgSecMomMatActive(SEXP adjMat, SEXP thetaMat, SEXP maxRuntime)
{
    // protect the items
    PROTECT(adjMat);
    PROTECT(thetaMat);
    
    int size = INTEGER(GET_DIM(adjMat))[0];
    time_t quittingTime = time(NULL) + INTEGER(maxRuntime)[0];
    
    // get space for the expectation and covariance matrix
    SEXP secMomMat;
    PROTECT(secMomMat = allocMatrix(REALSXP, size, size));

    double *secMomMatPtr = REAL(secMomMat);
    // initialize the matrix
    for(int i=0; i<size; ++i)
    {
        for(int j=0; j<size; ++j)
        {
            secMomMatPtr[i*size + j] = 0;
        }
    }

    // get the probgraph object
    try
    {
        ProbGraph graph(size, INTEGER(adjMat), REAL(thetaMat), quittingTime);
        graph.getExpectSecMomActive(secMomMatPtr);
    }
    catch(char const* errorMsg)
    {
        error(errorMsg);
    }
    
    // unprotect them again
    UNPROTECT(4);
    return(secMomMat);
}



// calculates the whole expectation vector but only the 
// covariances for one variables
SEXP runJTAlgSecMomVec(SEXP adjMat, SEXP thetaMat, SEXP varR, SEXP maxRuntime)
{
    // protect the items
    PROTECT(adjMat);
    PROTECT(thetaMat);
    PROTECT(varR);
    
    int size = INTEGER(GET_DIM(adjMat))[0];
    int var = INTEGER(varR)[0];
    time_t quittingTime = time(NULL) + INTEGER(maxRuntime)[0];

    
    // get space for the expectation and covariance matrix
    SEXP covVec, expec;
    PROTECT(covVec = allocVector(REALSXP, size));
    PROTECT(expec = allocVector(REALSXP, size));

    double *covVecPtr = REAL(covVec);
    double *expecPtr = REAL(expec);
    
    // get the probgraph object
    try
    {
        ProbGraph graph(size, INTEGER(adjMat), REAL(thetaMat), quittingTime);
        graph.getExpectSecMomSingleVar(var,expecPtr, covVecPtr);
    }
    catch(char const* errorMsg)
    {
        error(errorMsg);
    }
    
    // generate a list for the return values
    SEXP retList, dimnames;
    PROTECT(retList = allocVector(VECSXP,2));
    SET_VECTOR_ELT(retList,0,expec);
    SET_VECTOR_ELT(retList,1,covVec);
    PROTECT(dimnames = allocVector(STRSXP,2));
    SET_STRING_ELT(dimnames, 0, mkChar("Expectation"));
    SET_STRING_ELT(dimnames, 1, mkChar("SecondMomentVector"));
    setAttrib(retList,R_NamesSymbol, dimnames); 
    
    // unprotect them again
    UNPROTECT(8);
    return(retList);
}




// calculates the whole expectation vector but only the 
// covariances for one variables
SEXP runJTAlgNormalizationConstant(SEXP adjMat, SEXP thetaMat, SEXP maxRuntime)
{
    // protect the items
    PROTECT(adjMat);
    PROTECT(thetaMat);
    
    int size = INTEGER(GET_DIM(adjMat))[0];
    time_t quittingTime = time(NULL) + INTEGER(maxRuntime)[0];
    
    // get space for the expectation and covariance matrix
    SEXP normConst;
    PROTECT(normConst = allocVector(REALSXP, 1));

    double *normConstPtr = REAL(normConst);
    
    // get the probgraph object
    try
    {
        ProbGraph graph(size, INTEGER(adjMat), REAL(thetaMat), quittingTime);
        *normConstPtr = graph.getNormalizationConstant();
    }
    catch(char const* errorMsg)
    {
        error(errorMsg);
    }
    
    
    // unprotect them again
    UNPROTECT(4);
    return(normConst);
}




}
