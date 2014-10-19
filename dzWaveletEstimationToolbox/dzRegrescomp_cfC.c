

#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int r=0,  c=0;
    int rr=0, cc=0;
    int present_len=0, cuml_len=0;
	double *njj, *vjj,  *jj_safe, *gj, *jj_safe_num; // input
	double *output; //output
	double *tmp_njj, *tmp_vjj;
    mxArray *gammaln_mx;
    double  *gammaln_dp;
    mxArray *out1_mx, *out2_mx, *out3_mx;
    double  *out1_dp, *out2_dp, *out3_dp;
    int tmp_jj_safe=0;

    // Usage Check
	if(nlhs!=1||nrhs!=4)
	{
		mexPrintf("Usage: output=dzRegrescomp_cfC(njj,vjj,jj_safe,jj_safe_num)\n");
		return;
	}

    // Input
	njj=mxGetPr(prhs[0]);
	vjj=mxGetPr(prhs[1]);
	jj_safe=mxGetPr(prhs[2]);
	jj_safe_num=mxGetPr(prhs[3]);

    // Input Size
	r=mxGetM(prhs[1]); //vjj
	c=mxGetN(prhs[1]); //vjj

    // Output
	plhs[0]=mxCreateDoubleMatrix(1,c,mxREAL);
	output=mxGetPr(plhs[0]);
   

    // Calc
    for(cc=0;cc<c;++cc)
    {

    	present_len=(int)*(jj_safe_num+cc);

        // out1
        gammaln_mx=mxCreateDoubleMatrix(present_len,1,mxREAL);
        gammaln_dp=mxGetPr(gammaln_mx);
        for(rr=0;rr<present_len;++rr)
        {
            tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
            if(tmp_jj_safe<0)
            {
                mexPrintf("Error: The fourth parameter have zero or negative value\n");
                return;
            }
            *(gammaln_dp+rr)=*(vjj+cc*r+tmp_jj_safe) * 2 + *(njj+tmp_jj_safe)/2; // 2*vjj(jj_safe)+njj(jj_safe)/2
        }
        mexCallMATLAB(1,&out1_mx,1,&gammaln_mx,"gammaln");

        // out2
        for(rr=0;rr<present_len;++rr)
        {
        	tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
            *(gammaln_dp+rr)=*(njj+tmp_jj_safe) / 2; // njj(jj_safe)/2
        }
        mexCallMATLAB(1,&out2_mx,1,&gammaln_mx,"gammaln");

        // out3
        for(rr=0;rr<present_len;++rr)
        {
            tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
            *(gammaln_dp+rr)=*(vjj+cc*r+tmp_jj_safe) + *(njj+tmp_jj_safe)/2; // vjj(jj_safe)+njj(jj_safe)/2
        }
        mexCallMATLAB(1,&out3_mx,1,&gammaln_mx,"gammaln");

    	mxDestroyArray(gammaln_mx);
        
        out1_dp=mxGetPr(out1_mx);
        out2_dp=mxGetPr(out2_mx);
        out3_dp=mxGetPr(out3_mx);
    	for(rr=0;rr<present_len;++rr)
    	{
    		tmp_jj_safe=(int)*(jj_safe+cuml_len+rr);
            *(output+cc)+=*(out1_dp+rr) + *(out2_dp+rr) - *(out3_dp+rr)*2;
        }
        *(output+cc)=exp(*(output+cc))-1;
        
        //cumulate
        cuml_len+=present_len;

        //
        mxDestroyArray(out1_mx);
        mxDestroyArray(out2_mx);
        mxDestroyArray(out3_mx);

    }
}
