


#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
	int r=0,  c=0;
	int rr=0, cc=0;
	double *wjj, *njj, *vjj, *jj_safe, *jj_safe_num, *precis_psi;
	mxArray *tmp_mx;
    double *tmp_dp;
	double *output;
    int present_len=0, cuml_len=0;
    int tmp_jj_safe=0;
    mxArray *rhs[2];
    mxArray *psi1_mx;
    double *psi1_dp;
    mxArray *psi2_mx;
    double *psi2_dp;
    
    if(nlhs!=1||nrhs!=6)
	{
		mexPrintf("Usage: output=dzRegrescomp_CoValphacfC(wjj,njj,vjj,jj_safe,jj_safe_num,precis_psi)\n");
		return;
	}

    wjj=mxGetPr(prhs[0]);
	njj=mxGetPr(prhs[1]);
	vjj=mxGetPr(prhs[2]);
	jj_safe=mxGetPr(prhs[3]);
	jj_safe_num=mxGetPr(prhs[4]);
	precis_psi=mxGetPr(prhs[5]);

	r=mxGetM(prhs[2]);
	c=mxGetN(prhs[2]);

	plhs[0]=mxCreateDoubleMatrix(c,1,mxREAL);
	output=mxGetPr(plhs[0]);

    for(cc=0;cc<c;++cc)
    {
    	present_len=(int)*(jj_safe_num+cc);
        
        tmp_mx=mxCreateDoubleMatrix(present_len,1,mxREAL);
        tmp_dp=mxGetPr(tmp_mx);
        for(rr=0;rr<present_len;++rr)
        {
        	tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
            if(tmp_jj_safe<0)
            {
                mexPrintf("Error: The fourth parameter have zero or negative value\n");
                return;
            }            
        	*(tmp_dp+rr)=*(njj+tmp_jj_safe)/2 + *(vjj+cc*r+tmp_jj_safe);
        }
    
        rhs[0]=tmp_mx;
        rhs[1]=prhs[5];
        mexCallMATLAB(1,&psi1_mx,2,rhs,"dzpsi");


        for(rr=0;rr<present_len;++rr)
        {
        	tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
        	*(tmp_dp+rr)=*(njj+tmp_jj_safe)/2;
        }
        
        rhs[0]=tmp_mx;
        rhs[1]=prhs[5];
        mexCallMATLAB(1,&psi2_mx,2,rhs,"dzpsi");

        mxDestroyArray(tmp_mx);
        
        psi1_dp=mxGetPr(psi1_mx);
        psi2_dp=mxGetPr(psi2_mx);
    	for(rr=0;rr<present_len;++rr)
    	{
    		tmp_jj_safe=(int)*(jj_safe+cuml_len+rr)-1;
    		*(output+cc)+=*(wjj+cc*r+tmp_jj_safe) * (*(psi1_dp+rr) - *(psi2_dp+rr));
    	}
        mxDestroyArray(psi1_mx);
        mxDestroyArray(psi2_mx);
        
        
        cuml_len+=present_len;


    	/* Matlab 
    	sum(wjj(jj_safe) .* (psi(njj(jj_safe)/2+vjj(jj_safe),precis_psi)-psi(njj(jj_safe)/2,precis_psi)))
    	*/

    }



}