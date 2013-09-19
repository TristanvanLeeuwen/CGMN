#include "math.h"
#include "mex.h"   
/*
Kaczmarz sweep on band-matrix. The essential loop is:
for each row i
	x = x + w(i)*(b(i) - A(i,:)*x)*A(i,:)'
end
The matrix is given in banded storage format, where each column of the array S
stores a band of the matrix such that A(i,i+idx(j)) = S(i,j).

use:
	y = sweepR_mex(S,idx,x,b,w,dir)
 
 Tristan van Leeuwen, 2013
 Tristan.vanLeeuwen@gmail.com
*/

void init_zero(double *x, int n){
	int i;
	
	for(i=0;i<n;i++){
		x[i] = 0;
	}
}

void init_array(double *y, double *x, int n){
	int i;
	
	for(i=0;i<n;i++){
		y[i] = x[i];
	}
}

void do_sweep(int nrow, int ncol, int nx, double *Sr, double *Si, double *idx, double *yr, double *yi, double *br, double *bi, double *w, int dir){
	double cr, ci;
	int i,j,k;
	/*actual sweep, forward direction*/
	if(dir>0){
        /*i loops over the rows*/
		for(i=0;i<nrow;i++){
			/* first calculate the inner products*/
			cr = br[i];
			ci = bi[i];
			/*j loops over bands*/
			for(j=0;j<ncol;j++){
				/*k is the column index*/
				k = i + (int)idx[j];
				if((k>=0)&(k<nx)){
					cr -= Sr[i + j*nrow]*yr[k] - Si[i + j*nrow]*yi[k];
					ci -= Sr[i + j*nrow]*yi[k] + Si[i + j*nrow]*yr[k];
				}
			}
			/* now, update the vector */
            cr *= w[i];
            ci *= w[i];
			for(j=0;j<ncol;j++){
				k = i + (int)idx[j];
				if((k>=0)&(k<nx)){
					yr[k] +=   cr*Sr[i + j*nrow] + ci*Si[i + j*nrow];
					yi[k] +=  -cr*Si[i + j*nrow] + ci*Sr[i + j*nrow];
				}
			}
		}
	}
	/* same, in reverse direction*/
	else{
		for(i=nrow-1;i>=0;i--){
			/* first calculate the inner products*/
			cr = br[i];
			ci = bi[i];
			/*j loops over bands*/
			for(j=0;j<ncol;j++){
				/*k is the column index*/
				k = i + (int)idx[j];
				if((k>=0)&(k<nx)){
					cr -= Sr[i + j*nrow]*yr[k] - Si[i + j*nrow]*yi[k];
					ci -= Sr[i + j*nrow]*yi[k] + Si[i + j*nrow]*yr[k];
				}
			}
			/* now, update the vector */
            cr *= w[i];
            ci *= w[i];
			for(j=0;j<ncol;j++){
				k = i + (int)idx[j];
				if((k>=0)&(k<nx)){
					yr[k] +=   cr*Sr[i + j*nrow] + ci*Si[i + j*nrow];
					yi[k] +=  -cr*Si[i + j*nrow] + ci*Sr[i + j*nrow];
				}
			}
		}
	}  
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double *Sr,*Si,*xr,*xi,*br,*bi,*yr,*yi;
    double *idx, *w;
	int dir, nrow, ncol, nx;
	
	/* read input, initialize complex part to zero if input is real.*/
    ncol = mxGetN(prhs[0]);
	nrow = mxGetM(prhs[0]);
	Sr  = mxGetPr(prhs[0]);
    if(mxIsComplex(prhs[0])){
        Si  = mxGetPi(prhs[0]);
    }
    else{
        Si = mxCalloc(nrow*ncol,sizeof(double));
		init_zero(Si,nrow*ncol);
    }
	idx = mxGetPr(prhs[1]);	
    nx  = mxGetM(prhs[2]);
	xr  = mxGetPr(prhs[2]);
    if(mxIsComplex(prhs[2])){
        xi  = mxGetPi(prhs[2]);
    }
    else{
        xi = mxCalloc(nx,sizeof(double));
		init_zero(xi,nx);
    }
	br  = mxGetPr(prhs[3]);
    if(mxIsComplex(prhs[3])){
        bi  = mxGetPi(prhs[3]);
    }
    else{
        bi = mxCalloc(nrow,sizeof(double));
		init_zero(bi,nrow);
    }
	w   = mxGetPr(prhs[4]);
	dir = mxGetScalar(prhs[5]);
	
	/* define output vector y and initialize with input vector x.*/
	plhs[0]  = mxCreateDoubleMatrix(nx, 1, mxCOMPLEX); 
	yr       = mxGetPr(plhs[0]);
    yi       = mxGetPi(plhs[0]);
	init_array(yr,xr,nx);
	init_array(yi,xi,nx);
	
	/*sweep*/
	do_sweep(nrow, ncol, nx, Sr,Si,idx,yr,yi,br,bi,w,dir);
	
    return;
}



       
