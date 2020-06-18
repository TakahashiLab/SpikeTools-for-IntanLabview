 #include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

void makeOutputs(unsigned int *sq,short *raw,short *wv,int len,int step,int X){
  int i,j;
  int ind;
  int count=0;
	
  for(i=0;i<len;i++){
    ind=sq[i]-(step/4);
    
    for (j=0;j<X*step;j++){	  
      wv[count++]=raw[ind*X+j];
      //      printf("count=%d\n",count);
      //      printf("%d\n",ind*X+j);
    }
  }
  return;	
}

void ExistOrNot(short *raw,unsigned int *sp,int step,double* th,int len,int X,short *raw2){
  int i,j,k,l,m;
  int minP,maxP,minP2;
  short minV,maxV;
  int count=0;
  int Cross=0;
  int startPoint=0;
  int endPoint=0;
  int maxCh,minNum,maxNum;


  for(i=0;i<len;i++){

    Cross=0;
    //	  printf("i=%d\n",i);
    for (j=0;j<X;j++){
      /* threshold crossing for fish */
      if (raw[i*X+j] > th[j]){
	Cross=1;
	maxV=raw[i*X+j];
	break;
      }
    }

    
    if (Cross==1){

      /* search negative peak point */

      minV=32767;
      //	    endPoint=(step*X+(i+1)*X+j);
      startPoint=((-step/2)*X+(i+1)*X+j);

      if(startPoint<step/4){
	startPoint=step/4;
	i=step/4+step/2;
      }
      endPoint=startPoint+step*X;
	
      i=i-(step/2);
      for (k=startPoint;k<endPoint;k=k+X){
	if (raw2[k] < minV){
	  minP=i;
	  minV=raw2[k];
	}
	i++;
	if (i> len)
	  return;
      }
      //recheck
      sp[count++]=minP;
      i=minP+52;//plus 1msec
    }
  }

  return;	
}
	
void mexFunction(int nOUT, mxArray *pOUT[],int nINP, const mxArray *pINP[])
{
  	/* check number of arguments */
       if (nINP != 4 || nOUT !=2){
           mexErrMsgTxt("y=extractSpikesP(raw,x,xp,step,stdV,NoE)\n");
	   return;
       }

       int i;
       mwSize outlength;
       
	/* read inputs */
	short *raw=(short *)mxGetPr(pINP[0]);
        mwSize X=mxGetM(pINP[0]);
	mwSize Y=mxGetN(pINP[0]);

	double *th=mxGetPr(pINP[1]);

	mwSize step=(int)mxGetScalar(pINP[2]);
	//	int EstStep=(int)mxGetScalar(pINP[3]);


	short *raw2=(short *)mxGetPr(pINP[3]);

	int len=(int)floor((double)(Y/step));
	unsigned int *sq=(unsigned int *)malloc(sizeof(int)*len*3);
	memset(sq,-1,sizeof(unsigned int)*len*3);

	ExistOrNot(raw,sq,step,th,Y,X,raw2);

	/* create outputs */
	for (i=0,outlength=0;i<len*3;i++){
		if(sq[i]==-1){
			break;	
		}
		outlength++;
	}


	if(outlength==0){
          pOUT[0] = mxCreateDoubleMatrix(1,1, mxREAL);
          pOUT[1] = mxCreateDoubleMatrix(1,1, mxREAL);
          *mxGetPr(pOUT[0])=(double)-1;
          *mxGetPr(pOUT[1])=(double)-1;

  	  free(sq);
	  return;
	}

	//	printf("%d\n",outlength);
	
	
	mwSize wvDims[] = {X,outlength*step};
	pOUT[0] = mxCreateNumericArray(2, wvDims, mxINT16_CLASS, mxREAL);
	short *wv =(short*) mxGetPr(pOUT[0]);


	makeOutputs(sq,raw2,wv,outlength,step,X);	

        mwSize wvDimsSq[] = {1,outlength};
	pOUT[1] = mxCreateNumericArray(2, wvDimsSq, mxUINT32_CLASS, mxREAL);
	unsigned int *seq = (unsigned int *)mxGetPr(pOUT[1]);

	for (i=0;i<outlength;i++){
	  seq[i]=sq[i];
	}

	free(sq);
	
	return;	

}




