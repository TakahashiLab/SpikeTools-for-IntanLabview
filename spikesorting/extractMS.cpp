#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

void makeOutputs(double *wv,double *raw,int step,double *seq,int X,int seqLen){
	int i,j;	
	int ind=0;
	int center=(step/2)*X;	
	int base;
	int counter=0;
	for(i=0;i<seqLen;i++){

	  ind=seq[i];

	  base=ind*X;
	  //for (j=(ind-1)*step*X;j<ind*step*X;j++){
	  for (j=base-center;j<base+center;j++){
	    wv[counter]=(double)raw[j];
            counter++;
	    printf("%d\n",counter);
	  }
	}

	return;	
}

	
void mexFunction(int nOUT, mxArray *pOUT[],int nINP, const mxArray *pINP[])
{
    	/* check number of arguments */
       if (nINP != 3 || nOUT !=1)
           mexErrMsgTxt("y=getSpikes(x,step,cluster)\n");
       
	/* read inputs */
	double *raw=mxGetPr(pINP[0]);
	mwSize X=mxGetM(pINP[0]);
	mwSize Y=mxGetN(pINP[0]);

	int step=(int)mxGetScalar(pINP[1]);

	double *seq= mxGetPr(pINP[2]);
	mwSize seqX=mxGetM(pINP[2]);
	mwSize seqY=mxGetN(pINP[2]);


	/* create outputs */
	if (seqY==0){
      	  	const mwSize wvDims[] = {X,seqY*step};
		pOUT[0] = mxCreateNumericArray(2, wvDims, mxDOUBLE_CLASS, mxREAL);
		double *wv = (double *)mxGetPr(pOUT[0]);
	}
	else{
	  const mwSize wvDims[] = {X,seqY*step};
	  pOUT[0] = mxCreateNumericArray(2, wvDims, mxDOUBLE_CLASS, mxREAL);
	  double *wv = (double *)mxGetPr(pOUT[0]);

          makeOutputs(wv,raw,step,seq,X,seqY);
	}

}







