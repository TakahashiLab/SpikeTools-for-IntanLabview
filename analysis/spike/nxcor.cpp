#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

void XcorrCore(double *A,int lenA,double *B,int lenB,double *event,int
lenE,int dispRange,int dataRange,double *output,double *outputS,
double *count,int methodN){
	int i;	
	double **tmpA;
	double **tmpB;
	double *rawXcorr;
	double *shiftXcorr;	
	int j=0,k=0;
	double currentP;
	int shift;
	int tmp;
	
	//initialize paramters to calculate the cross correlation

	int displength=(dispRange*2+1);
	int datalength=(dataRange+1);

	tmpA=(double **)malloc(sizeof(double *)*(lenE));
	tmpB=(double **)malloc(sizeof(double *)*(lenE));
	for (i=0;i<lenE;i++){	
		tmpA[i]=(double *)malloc(sizeof(double)*(datalength));
		tmpB[i]=(double *)malloc(sizeof(double)*(datalength+dispRange*2));
		memset(tmpA[i],0,sizeof(double)*(datalength));
		memset(tmpB[i],0,sizeof(double)*(datalength+dispRange*2));	
	}
	rawXcorr=(double *)malloc(sizeof(double)*(displength));
	shiftXcorr=(double *)malloc(sizeof(double)*(displength));
	memset(rawXcorr,0,sizeof(double)*(displength));
	memset(shiftXcorr,0,sizeof(double)*(displength));
	count[0]=0.0;
	count[1]=0.0;	
	
	for (i=0;i<lenE;i++){
		currentP=event[i];
		
		//for spikeA (Basis)
		for(;j<lenA;j++){
			if ( A[j] >= (currentP) && A[j] <=(currentP+dataRange)){
				tmpA[i][(int)(A[j]-(currentP))]=tmpA[i][(int)(A[j]-(currentP))]+1;
				count[0]=count[0]+1.0;
			}
			else if( A[j] > (currentP+dataRange) ){
				break;
			}	
		}

		// for spikeB (Reference)
		for(;k<lenB;k++){
		
			if ( B[k] >= (currentP-dispRange) && B[k] <= (currentP+dataRange+dispRange)){
				tmpB[i][(int)(B[k]-(currentP-dispRange))]=tmpB[i][(int)(B[k]-(currentP-dispRange))]+1;
				count[1]=count[1]+1.0;
			}
			else if( B[k] > (currentP+dataRange+dispRange) ){
				break;
			}	
		}
	}	

	
	//calculate the raw cross correlogram
	for( i=0;i<lenE;i++){
		for (j=0;j<datalength;j++){
		    if (tmpA[i][j] != 0){
			   for (k=0 ; k < datalength+dispRange*2; k++){
			    if(tmpB[i][k]!=0){
			        tmp=k-j;
				    if (tmp >= 0 && tmp <=dispRange*2){
						rawXcorr[tmp]=rawXcorr[tmp]+tmpA[i][j]*tmpB[i][k];
					}	
				}
			   }
			}  
		}
	}

	//calculate the shift cross correlogram
	for( i=0;i<lenE;i++){
		if(i==lenE-1){
			shift=0;
		}
		else{
			shift=i+1;
		}

		for (j=0;j<datalength;j++){
		    if (tmpA[i][j] != 0){
			   for (k=0 ; k < datalength+dispRange*2; k++){
			    if(tmpB[shift][k]!=0){
			        tmp=k-j;
				    if (tmp >= 0 && tmp <=dispRange*2){
						shiftXcorr[tmp]=shiftXcorr[tmp]+tmpA[i][j]*tmpB[shift][k];
					}	
				}
			   }
			}  
		}
	}		
	



	for(i=0;i<displength;i++){
		switch(methodN){
			case 1:
				output[i]=rawXcorr[i];
				break;
			case 2:
				output[i]=shiftXcorr[i];
				break;
			case 3:
				output[i]=rawXcorr[i]-shiftXcorr[i];
				break;
		}
		outputS[i]=shiftXcorr[i];
	}


	//free paramters	
	for (i=0;i<lenE;i++){	
		free(tmpA[i]);
		free(tmpB[i]);
	}
	free(tmpA);
	free(tmpB);
	free(rawXcorr);
	free(shiftXcorr);

	return;
}

void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
double *output;
double *outputS;
double *count;
	/* read inputs */
	double *A=mxGetPr(pINP[0]);
	int lenA=mxGetN(pINP[0]);

	double *B=mxGetPr(pINP[1]);
	int lenB=mxGetN(pINP[1]);

	double *event=mxGetPr(pINP[2]);
	int lenE=mxGetN(pINP[2]);

//	printf("%d %d %d\n",lenA,lenB,lenE);	
	int dispRange=(int)mxGetScalar(pINP[3]);

	int dataRange=(int)mxGetScalar(pINP[4]);

	int methodN=(int)mxGetScalar(pINP[5]);

	int length=dispRange*2+1;

	pOUT[0] = mxCreateDoubleMatrix(1,length,mxREAL);
	output = mxGetPr(pOUT[0]);

	pOUT[1] = mxCreateDoubleMatrix(1,length,mxREAL);
	outputS = mxGetPr(pOUT[1]);

	pOUT[2] = mxCreateDoubleMatrix(1,2,mxREAL);
	count = mxGetPr(pOUT[2]);

	XcorrCore(A,lenA,B,lenB,event,lenE,dispRange,dataRange,output,outputS,count,methodN);

	return;	
}
