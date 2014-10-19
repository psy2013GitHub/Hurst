



#include "mex.h"
#include "matrix.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	double *raw_tr,*window_size,*cuml_tr,*FA,*Hurst;
	int r=0,c=0,rr=0,cc=0,bb=0,w=0,ww=0,www=0,rrr=0;
	int tmp_raw_offset=0,tmp_culm_offset=0,tmp_FA_offset;
    int tmp_window_size=0;
    int nblock=0;
	double tmp_X=0, tmp_Y=0, tmp_FA=0;
	double sum_X=0, sum_Y=0, sum_XY=0, sum_XX=0, sum_YY=0;
	double a=0.0, b=0.0;

    if(nrhs!=2)
    {
        mexPrintf("Usage: ...=Hurst_DDFA(raw_timeseries,window_size)");
        return;
    }


	raw_tr=mxGetPr(prhs[0]);
	window_size=mxGetPr(prhs[1]);
	w=mxGetM(prhs[1]);
	r=mxGetM(prhs[0]);
	c=mxGetN(prhs[0]);
    //mexPrintf("r=%d, c=%d\n",r,c);
    plhs[0]=mxCreateDoubleMatrix(1,c,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(w,c,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(r,c,mxREAL);
	Hurst=mxGetPr(plhs[0]);
    FA=mxGetPr(plhs[1]);
    cuml_tr=mxGetPr(plhs[2]);

    // Culmulated Sum
	for(cc=0;cc<c;++cc)
	{
        //average
        sum_Y=0;
        for(rr=0;rr<r;++rr)
        {
            tmp_raw_offset=cc*r+rr;
            sum_Y+=*(raw_tr+tmp_raw_offset);
        }
        sum_Y=sum_Y/((double)(r));
        //mexPrintf("average %lf\n",sum_Y);
		for(rr=0;rr<r;++rr)
		{
		 tmp_culm_offset=cc*r+rr;
         *(cuml_tr+tmp_culm_offset)=0;
		 for(rrr=0;rrr<=rr;++rrr)
		 {
		    tmp_raw_offset=cc*r+rrr;
            *(cuml_tr+tmp_culm_offset)+=*(raw_tr+tmp_raw_offset)-sum_Y; // X-mean(X)
		 }
         //mexPrintf("%d_%d: tmp_cuml=%f\n",cc,rr,*(cuml_tr+tmp_culm_offset));
		}
	}

    // detrend ---> Hurst calc
    for(ww=0;ww<w;++ww)
    {
    	tmp_window_size=*window_size++;
        if(tmp_window_size<=1)
        {
            mexPrintf("\n\nwindow_size must greater than 1\n\n");
            return;
        }
    	if(r%tmp_window_size)
    	{
    		mexPrintf("window_size not dividable: %d\n",tmp_window_size);
    		return;
    	}
        //mexPrintf("\n * window_size: %d *\n",tmp_window_size);
        nblock=r/tmp_window_size;
    	// detrend & FA Calc & Hurst Calc
    	for(cc=0;cc<c;++cc)
    	{
            tmp_FA=0;
            for(bb=0;bb<nblock;++bb)
            {
                //mexPrintf("\n bb= %d\n",bb);
            // detrend
                // In thoery, cuml_sum better caculated from bottom to top, for memory purpose
                sum_X=0;
                sum_Y=0;
                sum_XY=0;
                sum_XX=0;
                for(www=0;www<tmp_window_size;++www)
                {
                    tmp_culm_offset=cc*r+bb*tmp_window_size+www;
                    //mexPrintf("tmp_culm_offset=%d\n",cc*r+bb*tmp_window_size+www);
                    tmp_X=(double)(www+1);
                    tmp_Y=*(cuml_tr+tmp_culm_offset);
                    //mexPrintf("tmp_X=%f ,tmp_Y=%f\n",tmp_X,tmp_Y);
                    sum_Y+=tmp_Y;
                    sum_X+=tmp_X;
                    sum_XY+=tmp_X*tmp_Y;
                    sum_XX+=tmp_X*tmp_X;
                }
                //mexPrintf("sum_Y=%f, sum_X=%f, sum_XY=%f, sum_XX=%f\n",sum_Y,sum_X,sum_XY,sum_XX);
                b=(double)tmp_window_size*sum_XY-sum_X*sum_Y;
                b/=((double)(tmp_window_size)*sum_XX - sum_X*sum_X);
                a = sum_Y / ((double)(tmp_window_size)) - b * sum_X / ((double)(tmp_window_size));
                //mexPrintf("a=%f b=%f\n",a,b);
            // FA Calc
                sum_Y=0;
                sum_YY=0;
                for(www=0;www<tmp_window_size;++www)
                {
                    tmp_culm_offset=cc*r+bb*tmp_window_size+www;
                    tmp_Y=*(cuml_tr+tmp_culm_offset)-b*(www+1)-a;
                    sum_Y+=tmp_Y;
                    sum_YY+=tmp_Y*tmp_Y;
                }
                //mexPrintf("sum_Y=%f, sum_YY=%f, tmp_window_size=%d\n",sum_Y,sum_YY,tmp_window_size);
                tmp_FA+=sum_YY/((double)(tmp_window_size))-sum_Y*sum_Y/((double)(tmp_window_size*tmp_window_size));
            }
            tmp_FA/=(double)(nblock);
            tmp_FA=pow(tmp_FA,0.5);
        *(FA+cc*w+ww)=tmp_FA;
        //mexPrintf("cc=%d w=%d ww=%d :%d %.2f\n",cc,w,ww,cc*w+ww,tmp_FA);
        }
    }
    if(nlhs<2)
        mxDestroyArray(plhs[2]);
    
    // Hurst Calc
    window_size=mxGetPr(prhs[1]);
    for(cc=0;cc<c;++cc)
    {
        sum_X=0;
        sum_Y=0;
        sum_XY=0;
        sum_XX=0;
        for(ww=0;ww<w;++ww)
    		{
    			tmp_FA_offset=cc*w+ww;
    			tmp_X=*(window_size+ww);
                tmp_X=log(tmp_X);
    			tmp_Y=*(FA+tmp_FA_offset);
                tmp_Y=log(tmp_Y);
                //mexPrintf("tmp_X=%lf ,tmp_Y=%lf\n",tmp_X,tmp_Y);
    			sum_Y+=tmp_Y;
    			sum_X+=tmp_X;
    			sum_XY+=tmp_X*tmp_Y;
    			sum_XX+=tmp_X*tmp_X;
    		}

        *(Hurst+cc) =((double)(w))*sum_XY-sum_X*sum_Y;
        *(Hurst+cc)/=((double)(w))*sum_XX-sum_X*sum_X;
        //mexPrintf("Hurst: (%d*%f-%f*%f)/(%d*%f-%f*%f)=%lf\n",w,sum_XY,sum_X,sum_Y,w,sum_XX,sum_X,sum_X,(((double)(w))*sum_XY-sum_X*sum_Y)/(((double)(w))*sum_XX-sum_X*sum_X));
    	//a=sum_Y/w-b*sum_X/((double)w);
    }
    if(nlhs<3)
        mxDestroyArray(plhs[1]);
}