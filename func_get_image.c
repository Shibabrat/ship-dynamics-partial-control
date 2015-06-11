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
    int i,j,Nx,Ny;

    if (nlhs!=1)
        mexErrMsgTxt("Function requires two output argument.");
      
    if (nrhs!=2)
        mexErrMsgTxt("Function requires five input arguments.");

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
    
    /*Populate the domain inside the box with points*/
    long double deltaX = 
//     for (i = 1; i <= N; i++)
//     {
//         for (j = 1; j <= N; j++)
//         {
//             
//         }     
//             
//     }
    
    
    
}























