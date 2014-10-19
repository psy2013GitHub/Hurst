

#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int n=0;
  double n_double=0;
  int r=0, c=0;
  int rr=0,  cc=0;
  int poffset=0;
  double *q, *epsilon;
  mxArray *N_mx;
  double  *N_dp;
  double *zeta;
  
  // Usage Check
  if(nlhs!=1||nrhs!=2)
  {
  	mexPrintf("Usage: zeta=dzrzeta(q,epsilon)");
  	return;
  }

  // Input
  q=mxGetPr(prhs[0]);
  epsilon=mxGetPr(prhs[1]);
  
  // Input Size
  r=mxGetM(prhs[0]);
  c=mxGetN(prhs[0]);
  mexPrintf("size(zeta)=[%d %d]\n",r,c);
  
  // Output
  plhs[0]=mxCreateDoubleMatrix(r,c,mxREAL);
  zeta=mxGetPr(plhs[0]);

  // Calc  
  // 1, m: N=max(round(epsilon^(-1)-q-1)+1,0);
  N_mx=mxCreateDoubleMatrix(r,c,mxREAL);
  N_dp=mxGetPr(N_mx);
  for(cc=0;cc<c;++cc)
  	for(rr=0;rr<r;++rr)
  	{
  		poffset=cc*r+rr;
  		n_double=pow(*epsilon,-1.0)-*(q+poffset)-1; //
  		n=(int)n_double; //
  		n=(n_double-n<n+1-n_double)?n:(n+1); // round
  		*(N_dp+poffset)=(n+1>0)?(n+1):0; //maximum
  	}
  
  // 2, m: zeta(i,j)=sum(q(i,j)+(0:N(i,j).^(-2))+(q(i,j)+N(i,j)+1)^(-1)
  for(cc=0;cc<c;++cc)
  	for(rr=0;rr<r;++rr)
  	{
  		poffset=cc*r+rr;
  		*(zeta+poffset)=0;
  		for(n=0;n<*(N_dp+poffset)+1;++n)
  		{
  			*(zeta+poffset)+=pow(*(q+poffset)+(double)(n),-2.0);
    //%fprintf(1,'%12.0f  %8d %24.10e   %24.10e \n',q(i,j), N(i,j), 1/q(i,j), zeta(i,j));
  		}
  		*(zeta+poffset)+=pow(*(q+poffset)+*(N_dp+poffset)+1,-1.0);
  	}

  	//
    mxDestroyArray(N_mx);

}