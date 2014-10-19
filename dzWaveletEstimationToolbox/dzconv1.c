
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /*
         data, coloum as variable;
         filter, column vector;
      optional:
         shape, "full" as default

    */

    int r_indata=0,  c_indata=0, len_filter=0;
    int r_outdata=0;
    int c=0, r=0, cc=0, rr=0;
    mwSize lenstring=0;
    double *indata, *filter, *outdata;
    char *shape;
    mxArray *tmp_out_mx;
    double *tmp_out_dp;
    mxArray *vin_mx;
    double *vin_dp;
    mxArray *rhs[3];

    // Usage Check
	if(nlhs!=1||!(nrhs==2||nrhs==3))
    {
		mexPrintf("Usage1: data=dzconv1(data,filter)\nUsage2: data=dzconv1(data,filter,'full/same/'')\n");
        return;
    }
    if(nrhs==3&&!mxIsChar(prhs[2]))
    {
        mexPrintf("The second parameter must be a string");
        return;
    }

	// Input
    indata=mxGetPr(prhs[0]);
    filter=mxGetPr(prhs[1]);
    if(nrhs==3)
    {
        lenstring=mxGetN(prhs[2])*sizeof(mxChar)+1;
        shape=(char *)mxMalloc(lenstring);
        mxGetString(prhs[2],shape,lenstring);
    }
    // fix me
    else if(nrhs==2)
    {
        shape=(char *) malloc(sizeof(char)*5); //4 +'\0'
        *shape='f'; *(shape+1)='u'; *(shape+2)='l'; *(shape+3)='l'; *(shape+4)='\0'; // default full
    }

    // Input size
    r=mxGetM(prhs[0]);
    c=mxGetN(prhs[0]);
    len_filter=(mxGetM(prhs[1])>mxGetN(prhs[1]))?mxGetM(prhs[1]):mxGetN(prhs[1]);


    // Ouput size, max{r_indata+len_filter-1,r_indata,len_filter}
    r_outdata=r+len_filter-1;
    if(r_outdata<r)
        r_outdata=r;
    if(r_outdata<len_filter)
        r_outdata=len_filter;
    
    // Output
    plhs[0]=mxCreateDoubleMatrix(r_outdata,c,mxREAL);
    outdata=mxGetPr(plhs[0]);
    //tmp_out_mx=mxCreateDoubleMatrix(r_outdata,1,mxREAL);
    
    vin_mx=mxCreateDoubleMatrix(r,1,mxREAL);
    vin_dp=mxGetPr(vin_mx);
    rhs[0]=vin_mx;
    rhs[1]=prhs[1];
    if(nrhs==3)
        rhs[2]=prhs[2];

    // Conv Calc
    for(cc=0;cc<c;++cc)
    {

        // copy 1 column
    	for(rr=0;rr<r;++rr)
        {
            *(vin_dp+rr)=*(indata+cc*r+rr);
        }
        
        // conv1
        if(nrhs==2)
            mexCallMATLAB(1,&tmp_out_mx,2,rhs,"conv");
        else
            mexCallMATLAB(1,&tmp_out_mx,3,rhs,"conv");
        
        // cp m:outdata(1,:)=tmp_out
        tmp_out_dp=mxGetPr(tmp_out_mx);
    	for(rr=0;rr<r_outdata;++rr)
    		*(outdata+cc*r_outdata+rr)=*(tmp_out_dp+rr);

        mxDestroyArray(tmp_out_mx);
    }

    // free memory
    mxDestroyArray(vin_mx);
    mxFree(shape);

}
