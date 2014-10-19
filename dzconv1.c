
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
         shape, /f/

    */

    int r_indata=0,  c_indata=0, r_filter=0;
    int r_outdata=0;
    int cc=0, rr=0;
    double *filter,*outputdata;
    char *shape;
    //
	if(nlhs!=1||!(nrhs==2||nrhs==3))
    {
		mexPrintf("Usage: data=dzconv1(data,filter)");
        return -1;
    }

	//
    indata=mxGetPr(prhs[0]);
    filter=mxGetPr(prhs[1]);
    if(nrhs==3)
        shape=*(mxGetPr(plhs[2]));
    else if(prhs==2)
        shape=(char *) malloc(sizeof(char)*5); //4+'\0'
        shape="full"; // default full


    r_indata=mxGetM(prhs[0]);
    c_indata=mxGetN(prhs[0]);
    r_filter=mxGetM(prhs[1]);


    // r
    r=r_indata+r_filter-1;
    if(r_outdata<r_indata)
        r_outdata=r_indata;
    if(r_outdata<r_filter)
        r_outdata=r_filter;
    
    //create matrix
    plhs[0]=mxCreateDoubleMatrix();
    outdata=mxGetPr(plhs[0]);
    tmp_out=mxCreateDoubleMatrix(r_outdata,c_indata,mxREAL);
    tmp_vec=mxCreateDoubleMatrix(r_outdata,1,mxREAL);

    //conv calc
    for(cc=0;cc<c;++cc)
    {
        //cp
    	for(rr=0;rr<r;++rr)
    		*(tmp_vec+rr)=*(indata+cc*r+rr);

        //conv1
        mexCallMatlab(1,tmp_out,3,tmp_vec,filter,shape,'conv2');

        //cp
    	for(rr=0;rr<r;++rr)
    		*(outputdata+cc*r+rr)=*(tmp_output+rr);
    }

    // free memory
    mxDestroyArray(tmp_out);
    mxDestroyArray(tmp_vec);
    if(nrhs==2)
        free(shape);

}
