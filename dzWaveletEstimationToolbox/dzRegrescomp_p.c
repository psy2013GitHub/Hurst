

#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int r=0,  c=0;
    int rr=0, cc=0;
    int tmp_length=0, cuml_length=0;
	double *njj, *vjj,  *jj_safe, *gj, *jj_safe_num; // input
	double *p; //output
	double *tmp_njj, *tmp_vjj;
    mxArray *in1_gammaln_mx, *in2_gammaln_mx;
    double  *in1_gammaln_dp, *in2_gammaln_dp;
    mxArray *out1_mx, *out2_mx;
    double  *out1_dp, *out2_dp;
	double *tmp1, *tmp2, *tmp3;
    int tmp_jj_safe=0;

	if(nlhs!=1||nrhs!=5)
	{
		mexPrintf("Usage: p=dzRegrescomp_p(njj,vjj,gj,jj_safe,jj_safe_num)\n");
		return;
	}

	njj=mxGetPr(prhs[0]);
	vjj=mxGetPr(prhs[1]);
	gj =mxGetPr(prhs[2]);
	jj_safe=mxGetPr(prhs[3]);
	jj_safe_num=mxGetPr(prhs[4]);

	r=mxGetM(prhs[1]); //vjj
	c=mxGetN(prhs[1]); //vjj

	plhs[0]=mxCreateDoubleMatrix(1,c,mxREAL);
	p=mxGetPr(plhs[0]);

    cuml_length=0;
    for(cc=0;cc<c;++cc)
    {
    	tmp_length=(int)*(jj_safe_num+cc);
        // tmp_njj
    	in1_gammaln_mx=mxCreateDoubleMatrix(tmp_length,1,mxREAL);
    	in1_gammaln_dp=mxGetPr(in1_gammaln_mx);
        for(rr=0;rr<tmp_length;++rr)
        {
        	tmp_jj_safe=(int)*(jj_safe+cuml_length+rr)-1;
            if(tmp_jj_safe<0)
            {
                mexPrintf("Error: The fourth parameter have zero or negative value\n");
                return;
            }
            *(in1_gammaln_dp+rr)=*(njj+tmp_jj_safe)/2; // njj(jj_safe)/2
        }
        jj_safe=mxGetPr(prhs[3]); //reset
      
        // tmp_vjj
        in2_gammaln_mx=mxCreateDoubleMatrix(tmp_length,1,mxREAL);
        in2_gammaln_dp=mxGetPr(in2_gammaln_mx);
        for(rr=0;rr<tmp_length;++rr)
        {
        	tmp_jj_safe=(int)*(jj_safe+cuml_length+rr)-1;
        	*(in2_gammaln_dp+rr)=*(vjj+cc*r+tmp_jj_safe)+*(njj+tmp_jj_safe)/2; // njj(jj_safe)
        }
        jj_safe=mxGetPr(prhs[3]); //reset
        

    	mexCallMATLAB(1,&out1_mx,1,&in1_gammaln_mx,"gammaln");
        mxDestroyArray(in1_gammaln_mx);
        mexPrintf("%d\n",tmp_length);
    	mexCallMATLAB(1,&out2_mx,1,&in2_gammaln_mx,"gammaln");
        mxDestroyArray(in2_gammaln_mx);


        out1_dp=mxGetPr(out1_mx);
        out2_dp=mxGetPr(out2_mx);
        // p = exp( sum( gammaln(njj(jj_safe)/2) + vjj(jj_safe).*log(njj(jj_safe)./(2.^(1-gj(jj_safe)))) - gammaln(vjj(jj_safe) + njj(jj_safe)/2) ));
    	for(rr=0;rr<tmp_length;++rr)
    	{
    		tmp_jj_safe=(int)*(jj_safe+cuml_length+rr)-1;
            *(p+cc)+=*(out1_dp+rr)+*(vjj+cc*r+tmp_jj_safe)*log(*(njj+tmp_jj_safe)/pow(2.0,1.0-*(gj+tmp_jj_safe)))-*(out2_dp+rr);
        }
        *(p+cc)=exp(*(p+cc));
        
        //cumulate
        cuml_length+=tmp_length;

    }
}
