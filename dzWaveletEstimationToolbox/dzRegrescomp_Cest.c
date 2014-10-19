


#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	int len=0,  c=0;
	int l=0, cc=0;
    double *fr, *Calpha, *ft_r, *ft_i;
    double *output;
    double abs_ft2=0;

    // Usage
    if(nrhs!=3)
    	mexErrMsgTxt("Uage: dzRegrescomp(fr,Calpha,ft)\n");

    // Input
    fr=mxGetPr(prhs[0]);
    Calpha=mxGetPr(prhs[1]);
	ft_r=mxGetPr(prhs[2]);
	ft_i=mxGetPi(prhs[2]);


    // Input Size
    len=mxGetN(prhs[0]); //size(fr,2)
    c=mxGetN(prhs[1]); // size(Calpha,2)

    //Output
    plhs[0]=mxCreateDoubleMatrix(1,c,mxREAL);
    output=mxGetPr(plhs[0]);
    
    // Calc
    for(cc=0;cc<c;++cc)
    {
    	for(l=0;l<len;++l)
    	{
    		abs_ft2=pow(*(ft_r+l),2.0)+pow(*(ft_i+l),2.0);
    		*(output+cc)+=pow(*(fr+l),-1.0 * (*(Calpha+cc))) * abs_ft2;
    	}
    }

}