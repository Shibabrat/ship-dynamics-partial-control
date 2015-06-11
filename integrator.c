#include "mex.h"
#include "math.h"

/* This is the C routine that generates the grid of points and */
/* that carries out the integration of that grid of points. It */
/* uses a 4th order Runge-Kutta method that is applied to every */
/* point of the grid of points. This routine also discriminates */
/* which regions of the of the initial grid are going to be in */
/* the forbidden regions and which not (for example in the case */
/* of the Duffing oscillator the regions that are around the */
/* periodic attractors). */



// void flow(long double xi[],long double xo[])    /* This function contains the flow that we are going to use to find the safe sets. */
// {                                               /* In this example contains the Duffing oscillator but it could be any other system. */ 
// long double c1,c2,c3,rho,phi;
// 
// long double vi1,vi2,vi3,vf1,vf2,vf3;
// 
// c1=0.15;                                       /* The are the parameters of the Duffing oscillator. */
// c2=1;
// c3=1;
// rho=0.245;
// phi=1;
// 
// vi1=xi[0];
// vi2=xi[1];
// vi3=xi[2];
// 
// vf1=vi2;                                       /* Here is the definition of the Duffing oscillator. */
// vf2=-c1*vi2+c2*vi1-c3*pow(vi1,3)+rho*sin(vi3);
// vf3=1;
// 
// xo[0]=vf1;
// xo[1]=vf2;
// xo[2]=vf3;
// 
// return;   
// }

void flow(long double xi[],long double xo[])    /* This function contains the flow that we are going to use to find the safe sets. */
{                                               
long double c2,c3,c4,c5,b1,b2,HBar,omegaN;
long double alpha0, H, lambda, omegaBar, omegaE;

long double vi1,vi2,vi3,vf1,vf2,vf3;

/* Parameters for the ship roll model: Soliman,Thompson-1991 */
omegaN = 0.62;
omegaE = 0.527;
alpha0 = 0.73;
H = 18.054;
lambda = 221.94;
c2=0.1296/pow(omegaN,2);
c3=1.0368/pow(omegaN,2);
c4 = -4.059/pow(omegaN,2);
c5 = 2.4052/pow(omegaN,2);
b1 = 0.0043/omegaN;
b2 = 0.0225;
HBar = alpha0*M_PI*(H/lambda);
omegaBar = omegaE/omegaN;

vi1=xi[0];
vi2=xi[1];
vi3=xi[2];

vf1=vi2;                                      
vf2=-vi1 - c2*fabs(vi1)*vi1 - c3*pow(vi1,3) - c4*fabs(vi1)*pow(vi1,3) -c5*pow(vi1,5) - b1*vi2 - fabs(vi2)*vi2 + HBar*sin(vi3); 
vf3= omegaBar;

xo[0]=vf1;
xo[1]=vf2;
xo[2]=vf3;

return;   
}

void rkutta_safe_flag(int N,long double coord_input[], long double coord_output[], int *safeFlag)    /* This is the function that carries out the integration. */
{                                                                 
int steps,i;

long double x,y,h, omegaBar, omegaE,omegaN;

long double vn[3];

long double k1[3];
long double k2[3];
long double k3[3];
long double k4[3];

long double store[3];

steps=200;    /* This is the number of steps of our integrator. */
omegaE = 0.527;
omegaN = 0.62;
omegaBar = omegaE/omegaN;

h=(2*M_PI)/(omegaBar*steps);    /* This is the size of the step that we will use for the integration. */
*safeFlag = 1;                  

x=coord_input[0];
y=coord_input[1];

vn[0]=x;
vn[1]=y;
vn[2]=0;

for(i=0;i<N*steps;i++){    /* Here is the loop for the integration. In this case we take the stroboscopic map between 0 and 2*pi. */

flow(vn,k1);
 
store[0]=k1[0];
store[1]=k1[1];
store[2]=k1[2];
store[0]=store[0]*0.5*h;
store[1]=store[1]*0.5*h;
store[2]=store[2]*0.5*h;
store[0]=store[0]+vn[0];
store[1]=store[1]+vn[1];
store[2]=store[2]+vn[2];
 
flow(store,k2);
 
store[0]=k2[0];
store[1]=k2[1];
store[2]=k2[2];
store[0]=store[0]*0.5*h;
store[1]=store[1]*0.5*h;
store[2]=store[2]*0.5*h;
store[0]=store[0]+vn[0];
store[1]=store[1]+vn[1];
store[2]=store[2]+vn[2];
 
flow(store,k3);
 
store[0]=k3[0];
store[1]=k3[1];
store[2]=k3[2];
store[0]=store[0]*h;
store[1]=store[1]*h;
store[2]=store[2]*h;
store[0]=store[0]+vn[0];
store[1]=store[1]+vn[1];
store[2]=store[2]+vn[2];
 
flow(store,k4);
 
k2[0]=k2[0]*2;
k2[1]=k2[1]*2;
k2[2]=k2[2]*2;
 
k3[0]=k3[0]*2;
k3[1]=k3[1]*2;
k3[2]=k3[2]*2;
 
k1[0]=(k1[0]+k2[0]+k3[0]+k4[0])*h/6;
k1[1]=(k1[1]+k2[1]+k3[1]+k4[1])*h/6;
k1[2]=(k1[2]+k2[2]+k3[2]+k4[2])*h/6;
 
vn[0]=vn[0]+k1[0];
vn[1]=vn[1]+k1[1];
vn[2]=vn[2]+k1[2];

/*
 * Check if the trajectory has landed in the unsafe(forbidden) region, very rudimentary 
 * approach to event function in MATLAB.
 * If it is the unsafe(forbidden) region we set safe flag to '0'. If it is 
 * not in the forbidden region we set flag to '1'. 
*/  
if ((fabs(vn[0]) >= 0.88) && (*safeFlag != 0))
{
    *safeFlag = 0;
    // printf("%d\n",safeFlag);
//     printf("%Lf \t %Lf\n",vn[0],vn[1]);
}
else
{
    continue;
}

}

// printf("%d\n",safeFlag);
x=vn[0];
y=vn[1];

coord_output[0]=x;
coord_output[1]=y;

return;

}

void grid_generator(double N, double xi, double yi, double xf, double yf, double grid_points_base[], double grid_points_image[])   /* This is the function that generates the grids and that carries out the */ 
{                                                                                                                                  /* integration of every point of the grid.It also checks which points are */
                                                                                                                                   /* in the forbidden region which not of the initial set. */
int j,k,resolution,lock, safeFlag;
long double d1,d2,d3,d4,d5;
long double coord_input[2],coord_output[2];
long double t1[2],t2[2],t3[2],t4[2],t5[2];

resolution=N;  /* We pass as the argument N the resolution of the grid that will have N*N points. */

for(j=0;j<=resolution;j++){
     for(k=0;k<=resolution;k++){
            grid_points_base[5*(resolution+1)*k+j*5]=xi + (xf-xi)*j/N;                   /* Here we generate the grid of points. */
            grid_points_base[5*(resolution+1)*k+j*5+1]=yi + (yf-yi)*k/N;
            coord_input[0]= xi + (xf-xi)*j/N;
            coord_input[1]= yi + (yf-yi)*k/N;
            coord_output[0]=0;
            coord_output[1]=0;
            safeFlag = 1;
            rkutta_safe_flag(1,coord_input,coord_output, &safeFlag);                                   /* Here we carry out the integration of every point of the grid. */
//             printf("%Lf \t %Lf \n",coord_input[0],coord_input[1]);
//             printf("%Lf \t %Lf \n",coord_output[0],coord_output[1]);
//             printf("%d\n",safeFlag);
            
            grid_points_image[5*(resolution+1)*k+j*5]=coord_output[0];
            grid_points_image[5*(resolution+1)*k+j*5+1]=coord_output[1];
            
            lock=1;

//             t1[0]=0.815;
//             t1[1]=0.242;
//            
//             t2[0]=-0.933;
//             t2[1]=0.299;
//            
//             t3[0]=-1.412;
//             t3[1]=-0.137;
//            
//             t4[0]=-0.354;
//             t4[1]=-0.614;
//            
//             t5[0]=0.645;
//             t5[1]=-0.464;
//  
//            
//             t1[0]=t1[0]-coord_input[0];
//             t1[1]=t1[1]-coord_input[1];
//            
//             t2[0]=t2[0]-coord_input[0];
//             t2[1]=t2[1]-coord_input[1];
//            
//             t3[0]=t3[0]-coord_input[0];
//             t3[1]=t3[1]-coord_input[1];
//            
//             t4[0]=t4[0]-coord_input[0];
//             t4[1]=t4[1]-coord_input[1];
//            
//             t5[0]=t5[0]-coord_input[0];
//             t5[1]=t5[1]-coord_input[1];
//             
//            
//             d1=sqrt(t1[0]*t1[0]+t1[1]*t1[1]);
//             d2=sqrt(t2[0]*t2[0]+t2[1]*t2[1]);
//             d3=sqrt(t3[0]*t3[0]+t3[1]*t3[1]);
//             d4=sqrt(t4[0]*t4[0]+t4[1]*t4[1]);
//             d5=sqrt(t5[0]*t5[0]+t5[1]*t5[1]);
// 
// 
// 	        if (d1>0.2 && d2>0.2 && d3>0.2 && d4>0.2 && d5>0.2){     /* Here we check which point of the grid is in the forbidden region and which not. */      
//                 lock=1;                                              /* If it is the forbidden region we set lock to '0'. If it is not in the forbidden */
//             }                                                        /* region we set lock to '1'. */
//             else
//             {
//                 lock=0;
//             }
            
//             if (fabs(coord_input[0]) < 0.88)
//             {
//                 lock = 1;
//             }
//             else
//             {
//                 lock = 0;
//             }
//             
            grid_points_base[5*(resolution+1)*k+j*5+2]= safeFlag;                  
     }
}

}

void mexFunction( int nlhs, mxArray *plhs[],                        /* This is the function that sets up the interface between Matlab and our C code. */
                  int nrhs, const mxArray *prhs[] )
{
    
  double *xi,*yi,*xf,*yf,*grid,*grid_imagen;
  mwSize mrows,ncols,dimen[4];
  double *N;  
  int resolution;
  
  
  /* The input must be a noncomplex scalar double. */
  
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar.");
  }
  
  N = mxGetPr(prhs[0]);
  xi = mxGetPr(prhs[1]);
  yi = mxGetPr(prhs[2]);
  xf = mxGetPr(prhs[3]);
  yf = mxGetPr(prhs[4]);
  resolution=*N;
  
  dimen[0]=1;
  dimen[1]=5;
  dimen[2]=resolution;
  dimen[3]=resolution;
  
  /* Create matrixes for the return argument */
  
 
  plhs[0] = mxCreateNumericArray(4, dimen, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(4, dimen, mxDOUBLE_CLASS, mxREAL);
  
  /* Assign pointers to each input and output. */
  
  grid = mxGetPr(plhs[0]);
  grid_imagen = mxGetPr(plhs[1]);
 
  grid_generator((*N)-1,*xi,*yi,*xf,*yf,grid,grid_imagen);


  
}



