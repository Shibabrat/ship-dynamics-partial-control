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
#include "parameters.h"
#include "integration.h"

/* Input Arguments */

#define	tk      prhs[0]
#define	xIn     prhs[1]


/* Output Arguments */

#define	xOut    plhs[0]


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
                 /* variable declarations here */

{
//     int numPts, numDim;
    int i,j;
    double *currPos, *iterPos, *timeSpan;
    size_t m,n;
    double params = 0.1;
    double currPt[2], iterPt[2];
    
/* verify the number of input and output */
	if (nrhs != 2){
		mexErrMsgIdAndTxt("MyToolbox:mex_integration:nrhs",
                "This function takes 2 arguments.");
	}

	if (nlhs != 1){
		mexErrMsgIdAndTxt("MyToolbox:mex_integration:nlhs",
                "This function gives 1 output.");
	}

/* verify the type of arguments: array or not */
//     numPts = mxGetScalar(prhs[0]);
//     numDim = mxGetScalar(prhs[1]);
//     printf("Size of input array: %d x %d\n",numPts, numDim);

    m = mxGetM(xIn);
    n = mxGetN(xIn);
    printf("Size of input array: %d x %d\n",m, n);
    
    /* 3rd argument of time span is row vector of initial and final time */
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:mex_integration:notRowVector",
                      "Time span must be a row vector.");
    }

    /* 4th argument of grid points is an array of size #pts x #dim */
    /* make sure the first input argument is scalar */
    if( mxIsComplex(prhs[1]) ||
        mxGetNumberOfElements(prhs[1]) == 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:mex_integration:notArray",
                      "2nd input must be an array.");
    }
    
    if(mxGetM(prhs[1])!= m){
        mexErrMsgIdAndTxt("MyToolbox:mex_integration:notCorrectRows",
                      "Array has incorrect number of rows.");     
    }
    
    if(mxGetN(prhs[1]) != n){
        mexErrMsgIdAndTxt("MyToolbox:mex_integration:notCorrectColumns",
                      "Array has incorrect number of columns.");
    }

/* Extract time span data */
    timeSpan = mxGetPr(tk);
    double currT = timeSpan[0];
    double tau = timeSpan[1] - timeSpan[0];
    // printf("%lf\n",tau);

/* Get the pointer to the initial condition */
    currPos = mxGetPr(xIn);
    
/* Create a matrix for the return argument */
    xOut = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
       
/* Assign pointers to each input and output */
    iterPos = mxGetPr(xOut);

/* code here */
    
    /* Segment for integration grid of points */ 
    for (i = 0; i < m; i++){
        // for (j = 0; j < n; j++){
            currPt[0] = currPos[i];
            currPt[1] = currPos[i + m];
            printf("%lf \t %lf \n",currPt[0],currPt[1] );
            evolve_pt(currT, tau, currPt, iterPt, params);
            iterPos[i] = iterPt[0];
            iterPos[i + m] = iterPt[1];
            // printf("%lf \t %lf \n",iterPt[0],iterPt[1] );
        // }
    }

    return;
}
    
    
    
    
    