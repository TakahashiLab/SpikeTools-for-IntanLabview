#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

void RXcorrCore(double *A,int lenA,double *B,int lenB,int dispRange,double *rawAcorr){
	int i,j;	
	double reference,delays;
	int first,firstP;	
	
	//initialize paramters to calculate the cross correlation
	int displength=(dispRange*2+1);
	first=0;
	firstP=0;	
	for (i=0;i<lenA;i++){
		reference=A[i];
		for (j=firstP;j<lenB;j++){
			delays=B[j]-reference;
			if ( delays >= -dispRange && delays <=dispRange){

				if (first==0 && delays < -dispRange){
					first=j;	
				}
				rawAcorr[(int)(delays+dispRange)]=rawAcorr[(int)(delays+dispRange)]+1.0;
			}
			if (delays > dispRange){
				break;
			}	
		}
		firstP=first;
		first=0;	
	}	

	return;
}

void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	double *output;

	/* read inputs */
	double *A=mxGetPr(pINP[0]);
	int lenA=mxGetN(pINP[0]);

	double *B=mxGetPr(pINP[1]);
	int lenB=mxGetN(pINP[1]);

	int dispRange=(int)mxGetScalar(pINP[2]);
	int length=dispRange*2+1;	

	pOUT[0] = mxCreateDoubleMatrix(1,length,mxREAL);
	output = mxGetPr(pOUT[0]);
	memset(output,0,sizeof(double)*(length));

	RXcorrCore(A,lenA,B,lenB,dispRange,output);

	return;	
}
