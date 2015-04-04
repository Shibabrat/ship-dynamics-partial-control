/*
    Driver to solve ODES using GSL 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <matrix.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    
    double *phiMesh,*pPhiMesh;
    size_t N;
    mxArray *timeArray, *statesArray;
    int i,j;
/*   
    Step 1: Error Checking Step 1a: is nlhs 1?  If not,
    generate an error message and exit mexample (mexErrMsgTxt
    does this for us!) 
*/
    if (nlhs!=2)
        mexErrMsgTxt("Function requires two output argument.");
      
/*  Step 1b: is nrhs 2? */
    if (nrhs!=4)
        mexErrMsgTxt("Function requires four input arguments.");
      
/*   
    Step 1c: Is input array a numeric (nonstring), real (noncomplex) matrix? 
*/
//     if (mxIsChar(ARRAY_IN) || mxIsComplex(ARRAY_IN))
//         mexErrMsgTxt("First argument to mexample must be \
//         a numeric, real-valued array.");
      
  /*   Step 1d: Is VECTOR_IN a 2x1 or 1x2 numeric, full, real-
       valuedvector? */
    N=mxGetM(prhs[2]);  /* Assigning result of mxGetM and
                            mxGetN to variables. Assuming no of 
                            rows and columns are equal */
      
    long double initialTime, finalTime;   
    initialTime = mxGetScalar(prhs[0]);
    finalTime = mxGetScalar(prhs[1]);
   
    phiMesh = mxGetPr(prhs[2]);
    pPhiMesh = mxGetPr(prhs[3]);
    N = mxGetM(prhs[2]);
    if (N != mxGetM(prhs[3]))
    {
        mexErrMsgTxt("Size of mesh for both variables should be equal");
    }    
    if (N != mxGetN(prhs[3]))
    {
        mexErrMsgTxt("Size of mesh for both variables should be equal");
    }
    
    /*Begin solving the ODES*/
    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= N; j++)
        {
            
        }     
            
    }
    
    
    
}























