#include <mex.h>
#include "matrix.h"

int64_T size, mz_max, tot_counts;
double *dati;
double *prior;
double *stmp;
double *tmp;
int64_T mz, t ,ct1, ct2,i, counts;

void mexFunction(int nlhs, mxArray *plhs[ ],
int nrhs, const mxArray *prhs[ ])

{
    double *stmp;
    dati = mxGetPr(prhs[0]);
    size = mxGetM(prhs[0]);
    tot_counts=0;
    
    for(i=0;i<size; i++)
    {
        if (dati[i]>1)
            tot_counts+=abs((int64_T) (dati[i]+0.5));
         //counts = (int64_T)(dati[i] + 0.5 );
    }
     mexPrintf("Tot counts= %d\n", tot_counts);
    
    mz_max=mxGetScalar(prhs[1]);
    //tot_counts=mxGetScalar(prhs[2]);
    if (nrhs>2)
    {
        prior= mxGetPr(prhs[2]);
        stmp=(double*)mxCalloc(tot_counts, sizeof(double));
    }
    
    
    tmp=(double*)mxCalloc(tot_counts*2, sizeof(double));
    
    mz = 0;
    t = 1;
    ct1 = 0;
    
    for(i = 0; i < size; i++ )
    {
        mz++;
        if (mz>mz_max) // # rows matrix int
        {
            mz = 1;
            t++;
        }
        //mexPrintf("Dati[%d]=%f, t=%d, mz=%d\n", i, dati[i], t, mz);
        if ( dati[i]>1 )
        {
            
            counts = (int64_T)(dati[i] + 0.5 );
            for (ct2 = 0; ct2 < counts; ct2++)
            {
                tmp[ct1+ct2] = (double)t;
                tmp[ct1+ct2+tot_counts] = (double) mz;
                if (nrhs>2) //if there is a prior
                {
                    stmp[ct1+ct2] = prior[i];
                }
            }
            ct1 = ct1 + counts;
        }
    }
    plhs[0]=mxCreateDoubleMatrix(tot_counts,2,mxREAL);
    
    mxSetPr(plhs[0], tmp);
    if (nlhs>1)
    {
        plhs[1]=mxCreateDoubleMatrix(tot_counts,1,mxREAL);
        if (nrhs>2)  //if there is a prior
        {
            mxSetPr(plhs[1], stmp);
        }
    }
    //delete[] tmp;
    
}