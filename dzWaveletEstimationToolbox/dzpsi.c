

#include "mex.h"
#include "matrix.h"
#include <math.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	int N=0;
	int r=0, c=0;
	int rr=0,cc=0;
	int i=0;
	double *x, *epsilon;
	double *Psi;
	double z=0;

    // Usage Check
	if(nlhs!=1||nrhs!=2)
	{
		mexPrintf("Usage: data=dzpsi(data,epsilon)\n");
		return;
	}
    
    // Input
    x=mxGetPr(prhs[0]);
	epsilon=mxGetPr(prhs[1]);

	// Input Size
	r=mxGetM(prhs[0]);
    c=mxGetN(prhs[0]);

    // Output
    plhs[0]=mxCreateDoubleMatrix(r,c,mxREAL);
    Psi=mxGetPr(plhs[0]);

    // Calc
    N=252*pow(*epsilon,-1.0/6.0)+1; //ceil
	for(rr=0;rr<r;++rr)
		for(cc=0;cc<c;++cc)
		{
			z=*(x+r*cc+rr);
			if(z>N)
				*(Psi+r*cc+rr) =log(z) - 1.0/2.0/z - 1.0/12.0/pow(z,2.0) + 1.0/120.0/pow(z,4.0) - 1.0/252.0/pow(z,6.0);
			else
			{
				*(Psi+r*cc+rr) =log((z+N)) - 1.0/2.0/(z+N) - 1.0/12.0/pow(z+N,2.0) + 1.0/120.0/pow(z+N,4.0) - 1.0/252.0/pow(z+N,6.0);
				for(i=0;i<N;++i)
					*(Psi+r*cc+rr)-=pow(z+i,-1.0);
			}
		}
  //fprintf(1,'%12.5f  %8d %24.10e   %24.10e \n',x(i,j), N, log(z)     - 1/2/z     - 1/12/z^2     + 1/120/z^4     - 1/252/z^6, Psi(i,j));
}