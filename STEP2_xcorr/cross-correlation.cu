#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <math.h>

#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas_v2.h>

#include "segy_PASSCAL.h"
#include "../params.h"
#include "sub_crosscorrelation.cu"
#include "xcorKernel.cu"

int main(int argc, char *argv[]){
	struct segy segy;
	struct filelist *fl;
	int deviceid;
	int i, j, k;
	int index, ichan, iseg;
	int ichan_array1, ichan_array2; 
	int segyfile_number = 0;	
	int nt_array1, nt_array2;
	int nt;
	int nstat_array1, nstat_array2;
	int nstat_currentloop_array1, nstat_currentloop_array2;		// add for loop
	int nloop_array1, nloop_array2;					// add for loop
	int iloop_array1, iloop_array2;					// add for loop
	int iBatch_array1_begin, iBatch_array1_end;			// add for loop
	int Batch_currentloop_array1, Batch_currentloop_array2;		// add for loop
	int fd_array1, fd_array2;
	int fd_output;
	double delta_array1, delta_array2;
	double delta;
	int npts_lag, npts_cc;
	int npts_segment, seg_number;
	int npts_fft, npts_fft2;
//	int Batch_array1, Batch_array2;
	int npts_halfram, npts_whiten;
	float *waveform_temp_array1, *waveform_temp_array2;
//	float *cc;
	float **cc_stack;
	float *norm_array1, *norm_array2;

	float *weight_array1, *weight_array2;
	float df;
	int   Nlow, Nhigh;

	char filename[256];
	char filename_array1[256],filename_array2[256];
	char filename_output[256];
	char buf[3600];
	FILE *fp;

	int block_num = 1024;
	int grid_num_array1, grid_num_array2;
	float *dnorm_array1, *dnorm_array2;
        cublasHandle_t handle;
        cufftHandle plan_array1, plan_array2;
	cufftHandle plan_whiten_array1, plan_whiten_array2;
	cufftComplex *waveform_array1_orig, *waveform_array2_orig;
	cufftComplex *waveform_array2_orig_temp;			// add for loop
	cufftComplex *dwaveform_array1_orig, *dwaveform_array2_orig;
	cufftComplex *waveform_array1, *waveform_array2;
        cufftComplex *dwaveform_array1, *dwaveform_array2;
        cufftComplex *dxcorr;
//	cufftComplex *xcorr;
	float *dcc;
	

//-------------------------------------------------------------------------
// read segy file names from SEGY_LIST
	if((fp=fopen(SEGY_LIST,"r")) == NULL) perror(SEGY_LIST);
        while(fscanf(fp, "%s", filename) != EOF) segyfile_number += 1;
        rewind(fp);
        fl = (struct filelist *)malloc(sizeof(struct filelist) * segyfile_number);
        for(i=0;i<segyfile_number;i++){
           fscanf(fp, "%s", filename);
           strcpy(fl[i].file, filename);
	   convert_datetime2ut(fl[i].file,&fl[i].ut);
        }
        printf("  INPUT: %d HOUR SEGY FILES (SORTED BY TIME)\n", segyfile_number);
        fclose(fp);
	if(segyfile_number >= 2){
	   qsort(fl,segyfile_number,sizeof(struct filelist),filecomp);
           printf("\t\t%s: %.2f --> %s: %.2f\n\n",fl[0].file,fl[0].ut,fl[segyfile_number-1].file,fl[segyfile_number-1].ut);
	}
	else printf("\tSEGY FILES: \n\t%s: %.2f\n\n",fl[0].file,fl[0].ut);

//-------------------------------------------------------------------------
// read segy information from SEGY_INFO 
        printf("  SEGY INFO from %s\n", SEGY_INFO);
	if((fp=fopen(SEGY_INFO,"r")) == NULL) perror(SEGY_INFO);
	while(fgetc(fp) != '\n') i++;
	fscanf(fp, "%s%d%d%lf", buf, &nstat_array1, &nt_array1, &delta_array1);
	fscanf(fp, "%s%d%d%lf", buf, &nstat_array2, &nt_array2, &delta_array2);
	printf("\t        nstat   nt    delta\n");
	printf("\tArray1: %5d %d %.5f\n",nstat_array1, nt_array1, delta_array1);
	printf("\tArray2: %5d %d %.5f\n\n",nstat_array2, nt_array2, delta_array2);
	fclose(fp);

//-------------------------------------------------------------------------
// final check segy files before cross-correlation
        printf("  CHECKING SEGY FILES ...\n");
	if(fabs(delta_array1-delta_array2) > 1e-5){
	   printf("Error: delta between array1 (%f) and array2 (%f) are different!\n",delta_array1, delta_array2);
	   exit(0);
	}
	else{
	   delta = delta_array1;
	}

//	for(i=0;i<segyfile_number;i++){
//	   sprintf(filename_array1, "%s/%s", ARRAY1_DOWNSAMPLE_DIR,fl[i].file);
//	   checksegyfiles(filename_array1, delta_array1, nt_array1, nstat_array1);
//	   sprintf(filename_array2, "%s/%s", ARRAY2_DOWNSAMPLE_DIR,fl[i].file);
//	   checksegyfiles(filename_array2, delta_array1, nt_array1, nstat_array2);
//	   printf("\t%s okay!\n",fl[i].file);
//	}

	if(nt_array1 != nt_array2){
	   nt = nt_array1 < nt_array2 ? nt_array1 : nt_array2;
	   printf("Warning: number of points between array1 (%d) and array2 (%d) are different! NPTS %d will be used for Xcorr!\n",nt_array1, nt_array2, nt);
	}
	else{ 
	   nt = nt_array1;
	}
	printf("\tAll SEGY FILES LOOK OKAY!\n\tPROCEED TO CROSS-CORRELATION!\n\n");

//-------------------------------------------------------------------------
// start cross-correlation loop
        printf("  START CROSS-CORRELATION ...\n");
	if(argc >= 2){deviceid = atoi(argv[1]);}
	else{deviceid = 0;}  // default GPU ID: 0
	if(mkdir(XCORR_OUTPUT_DIR,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) != 0) printf("\t%s already exist or Creating Wrong !\n",XCORR_OUTPUT_DIR);

	npts_lag = floor(MAX_LAG / delta);
	npts_cc = 2 * npts_lag + 1;
	npts_segment = floor(XCOR_SEGMENT / delta) + 1;
	npts_whiten = 2;
        while(npts_whiten < nt) npts_whiten *= 2;
	npts_fft = 2;
        while(npts_fft < npts_segment) npts_fft *= 2;
        npts_fft2 = npts_fft;
	npts_fft *= 2;
	seg_number = floor((nt-1) / npts_segment);
	if(seg_number == 0){
	   printf("Error: XCOR_SEGMENT TIME IN PARAMS.H %d is longer than 1 hour\n",XCOR_SEGMENT);
	   exit(0);
	}
	nloop_array1 = ceil((double)nstat_array1/(double)ARRAY1_CHANNELPERLOOP);
       	nloop_array2 = ceil((double)nstat_array2/(double)ARRAY2_CHANNELPERLOOP); 

	printf("\tXCOR OUTPUT LENGTH: %.2f ~ %.2fs, npts_cc: %d\n",-1.0*MAX_LAG, 1.0*MAX_LAG, npts_cc);
	printf("\tXCOR SEGMENT WINDOW LENGTH: %.2fs, npts_seg: %d, npts_fft: %d\n", XCOR_SEGMENT*1.0, npts_segment, npts_fft);
	printf("\tTOTAL TIME PER FILE: %.2fs, SEGMENT NUMBER PER CHANNEL: %d\n", (nt-1)*delta, seg_number);
	printf("\tTOTAL NPTS PER CHANNEL: %d, npts_whiten: %d\n", nt, npts_whiten);
	printf("\tXCORR LOOPS FOR ARRAY1: %d \n", nloop_array1);
	printf("\tXCORR LOOPS FOR ARRAY2: %d \n", nloop_array2);
	printf("\tGPU DEVICE ID: %d\n", deviceid);

//	Batch_array1 = nstat_array1 * seg_number;
//	Batch_array2 = nstat_array2 * seg_number;
//	cc = (float *)malloc(sizeof(float) * nstat * npts_cc);
//	xcorr = (cufftComplex *)malloc(sizeof(cufftComplex) * npts_fft * Batch);
	waveform_array1_orig = (cufftComplex *)malloc(sizeof(cufftComplex) * npts_whiten * nstat_array1);
	waveform_array2_orig = (cufftComplex *)malloc(sizeof(cufftComplex) * npts_whiten * nstat_array2);
	if(waveform_array1_orig == NULL) printf("Error in mallocing waveform_array1_orig\n");
	if(waveform_array2_orig == NULL) printf("Error in mallocing waveform_array2_orig\n");

	cc_stack = (float **) malloc(sizeof(float *) * nstat_array1);
	for(i=0; i<nstat_array1; i++){	
	   cc_stack[i] = (float *)malloc(sizeof(float) * nstat_array2 * npts_cc);
	   if(cc_stack[i] == NULL) printf("Error in ccstack %d\n",i);
	   for(k=0; k<nstat_array2*npts_cc; k++) cc_stack[i][k] = 0.0;
	}
	

	cudaSetDevice(deviceid);
//	cufftCreateCheck(cufftPlan1d(&plan, npts_fft, CUFFT_C2C, Batch));

	sprintf(filename_output, "%s/stack.segy", XCORR_OUTPUT_DIR);
	if((fd_output = creat(filename_output,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0){
           printf("cannot creat %s\n",filename_output);
	   exit(0);}

	for(i=0;i<segyfile_number;i++){
	   printf("\tlooping segyfile : %d %s\n",i, fl[i].file);
// read waveform data from segy files
// copy original downsampled data
	   waveform_temp_array1 = (float *)malloc(sizeof(float) * nt);
	   sprintf(filename_array1, "%s/%s", ARRAY1_DOWNSAMPLE_DIR,fl[i].file);
           fd_array1 = open(filename_array1, 0);
	   read(fd_array1, buf, 3600);
	   if(i==0) write(fd_output, buf, 3600);
	   for(ichan=0; ichan<nstat_array1; ichan++){
	     read(fd_array1, &segy, SEGY_HEAD_SIZE);
	     read(fd_array1, waveform_temp_array1, sizeof(float)*nt);
	     for(j=0; j<npts_whiten; j++){
	       if(j<nt){
		 waveform_array1_orig[ichan*npts_whiten+j].x = waveform_temp_array1[j];
		 waveform_array1_orig[ichan*npts_whiten+j].y = 0.0;
	       }
	       else{
		 waveform_array1_orig[ichan*npts_whiten+j].x = 0.0;
		 waveform_array1_orig[ichan*npts_whiten+j].y = 0.0;
	       }
	     }
	   }
           close(fd_array1);
	   free(waveform_temp_array1);

	   waveform_temp_array2 = (float *)malloc(sizeof(float) * nt);
	   sprintf(filename_array2, "%s/%s", ARRAY2_DOWNSAMPLE_DIR,fl[i].file);
           fd_array2 = open(filename_array2, 0);
	   read(fd_array2, buf, 3600);
	   for(ichan=0; ichan<nstat_array2; ichan++){
	     read(fd_array2, &segy, SEGY_HEAD_SIZE);
	     read(fd_array2, waveform_temp_array2, sizeof(float)*nt);
	     for(j=0; j<npts_whiten; j++){
	       if(j<nt){
		 waveform_array2_orig[ichan*npts_whiten+j].x = waveform_temp_array2[j];
		 waveform_array2_orig[ichan*npts_whiten+j].y = 0.0;
	       }
	       else{
		 waveform_array2_orig[ichan*npts_whiten+j].x = 0.0;
		 waveform_array2_orig[ichan*npts_whiten+j].y = 0.0;
	       }
	     }
	   } // end for ichan
           close(fd_array2);
	   free(waveform_temp_array2);

//////////////////////////////////
// start xcorr loop for array1
	   for(iloop_array1=0; iloop_array1<nloop_array1; iloop_array1++){
	     if(iloop_array1 != (nloop_array1-1)){
	       nstat_currentloop_array1 = ARRAY1_CHANNELPERLOOP;
	     }
	     else{
	       nstat_currentloop_array1 = nstat_array1 - ARRAY1_CHANNELPERLOOP * (nloop_array1-1);
	     }

	     
// STEP 1. for array1
// Copy Original data To GPU for Temporal normalization and Spectral whitening
	     if(cudaMalloc((void **)&dwaveform_array1_orig,sizeof(cufftComplex)*npts_whiten*nstat_currentloop_array1) != cudaSuccess) printf("Cuda error: Failed to allocate for dwaveform_array1_orig\n");
	     if(cudaMemcpy(dwaveform_array1_orig, (waveform_array1_orig+npts_whiten*iloop_array1*ARRAY1_CHANNELPERLOOP), (nstat_currentloop_array1*npts_whiten*sizeof(cufftComplex)), cudaMemcpyHostToDevice) != cudaSuccess) printf("Error in copying waveform_array1_orig\n");
	     grid_num_array1 = floor(1.0 * nstat_currentloop_array1 * npts_whiten / block_num) + 1;
	     checkNanKernel<<<grid_num_array1,block_num>>>(dwaveform_array1_orig, npts_whiten, nstat_currentloop_array1);


// STEP 2. for array1
// TEMPORAL NORMALIZATION: onebit or running-absolute-mean
	     if(TEMPORAL_NORMALIZATION == 1){
	       Onebit<<<grid_num_array1,block_num>>>(dwaveform_array1_orig, npts_whiten, nstat_currentloop_array1);
	     }
	     else if(TEMPORAL_NORMALIZATION == 2){
	       npts_halfram = floor(1.0 / freql / delta / 2.0);
	       if(cudaMalloc((void **)&weight_array1,sizeof(float)*npts_whiten*nstat_currentloop_array1) != cudaSuccess) printf("Cuda error: Failed to allocate for weight_array1\n");
	       RunabsmeanKernel<<<grid_num_array1,block_num>>>(dwaveform_array1_orig, weight_array1, npts_whiten, nstat_currentloop_array1, nt, npts_halfram);
	       cudaSyncCheck(cudaDeviceSynchronize());
	       divweight<<<grid_num_array1,block_num>>>(dwaveform_array1_orig, weight_array1, npts_whiten, nstat_currentloop_array1, nt, npts_halfram);
	       cudaSyncCheck(cudaDeviceSynchronize());
	       cudaFree(weight_array1);
	     }
	     else if(TEMPORAL_NORMALIZATION != 0){
	       printf("Unknow TEMPORAL_NORMALIZATION defined in params.h!\n");
	       exit(0);
	     }


// STEP 3. for array1
// SPECTRAL WHITHENING 
	     if(SPECTRAL_WHITENING){
	       df = (1.0 / delta) / npts_whiten;
	       Nlow = floor((double)freql / (double)df);
	       Nhigh = ceil((double)freqh / (double)df);
	       
	       cufftCreateCheck(cufftPlan1d(&plan_whiten_array1, npts_whiten, CUFFT_C2C, nstat_currentloop_array1));
	       cufftExeInverseCheck(cufftExecC2C(plan_whiten_array1, dwaveform_array1_orig, dwaveform_array1_orig, CUFFT_FORWARD));
	       cudaSyncCheck(cudaDeviceSynchronize());
	       SpectralWhitenKernel<<<grid_num_array1,block_num>>>(dwaveform_array1_orig, Nlow, Nhigh, npts_whiten, nstat_currentloop_array1);
	       cufftExeInverseCheck(cufftExecC2C(plan_whiten_array1, dwaveform_array1_orig, dwaveform_array1_orig, CUFFT_INVERSE));
	       cudaSyncCheck(cudaDeviceSynchronize());
	       cufftDestroy(plan_whiten_array1);
	     } 
	     cudaMemcpy((waveform_array1_orig+npts_whiten*iloop_array1*ARRAY1_CHANNELPERLOOP), dwaveform_array1_orig, (npts_whiten*nstat_currentloop_array1*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
	     cudaFree(dwaveform_array1_orig);

	     
// STEP 4. for array1
// RESHAPE WAVEFORM TO SEGMENTED DATA IN CPU FOR LATER FFT AND XCOR
	     waveform_array1 = (cufftComplex *)malloc(sizeof(cufftComplex)*npts_fft*nstat_currentloop_array1*seg_number);
	     iBatch_array1_begin = iloop_array1*ARRAY1_CHANNELPERLOOP*seg_number;
	     iBatch_array1_end = (iloop_array1*ARRAY1_CHANNELPERLOOP+nstat_currentloop_array1)*seg_number;
	     for(j=iBatch_array1_begin; j<iBatch_array1_end; j++){
	       ichan = floor(1.0 * j / seg_number);
	       iseg = j - ichan * seg_number;
      	       index = ichan * npts_whiten + iseg * npts_segment;
	       for(k=0; k<npts_segment; k++){
		 waveform_array1[(j-iBatch_array1_begin)*npts_fft+k].x = waveform_array1_orig[index+k].x;
	         waveform_array1[(j-iBatch_array1_begin)*npts_fft+k].y = 0.0;
	       }
	       for(k=npts_segment; k<npts_fft; k++){
	         waveform_array1[(j-iBatch_array1_begin)*npts_fft+k].x = 0.0;
	         waveform_array1[(j-iBatch_array1_begin)*npts_fft+k].y = 0.0;
	       } // end for k
	     } // end for j


// STEP 5. for array1
// COPY TO GPU FOR XCORR
	     Batch_currentloop_array1 = nstat_currentloop_array1 * seg_number;
	     if(cudaMalloc((void **)&dwaveform_array1,sizeof(cufftComplex) * npts_fft * Batch_currentloop_array1) != cudaSuccess) printf("Cuda error: Failed to allocate for dwaveform_array1\n");
	     cudaMemcpy(dwaveform_array1, waveform_array1, (Batch_currentloop_array1 * npts_fft * sizeof(cufftComplex)), cudaMemcpyHostToDevice);
	     grid_num_array1 = floor(1.0 * Batch_currentloop_array1 * npts_fft / block_num) + 1;


// STEP 6. for array1
// CALCULATE L2 NORM FOR EACH SEGMENT FOR LATER NORMALIZED XCORR
	   norm_array1 = (float *)malloc(sizeof(float) * Batch_currentloop_array1);
	   if(norm_array1 == NULL) printf("Error in mallocaing norm_array1\n");
	   if(cudaMalloc((void **)&dnorm_array1,sizeof(float)*Batch_currentloop_array1) != cudaSuccess) printf("Cuda error: Failed to allocate for dnorm_array1\n");
	   cublasCreateCheck(cublasCreate(&handle));
	   for(j=0; j<Batch_currentloop_array1; j++)  cublasSnrmCheck(cublasSnrm2(handle, npts_segment, (cufftReal *)(dwaveform_array1+j*npts_fft), 2, &norm_array1[j]));
	   cublasDestroyCheck(cublasDestroy(handle));
	   cudaMemcpy(dnorm_array1, norm_array1, (Batch_currentloop_array1*sizeof(float)), cudaMemcpyHostToDevice);
//	   for(j=0; j<10; j++) printf("ARRAY1: %f \n",norm_array1[j]);
	   free(norm_array1);


// STEP 7. for array1
// FAST FOURIER TRANSFORM FOR EACH SEGMENT
	   cufftCreateCheck(cufftPlan1d(&plan_array1, npts_fft, CUFFT_C2C, Batch_currentloop_array1));
	   cufftExeForwardCheck(cufftExecC2C(plan_array1, dwaveform_array1, dwaveform_array1, CUFFT_FORWARD));
	   cudaSyncCheck(cudaDeviceSynchronize());
	   cufftDestroy(plan_array1);


//////////////////////////////////
// start xcorr loop for array2
	     for(iloop_array2=0; iloop_array2<nloop_array2; iloop_array2++){
	       if(iloop_array2 != (nloop_array2-1)){
		 nstat_currentloop_array2 = ARRAY2_CHANNELPERLOOP;
	       }
	       else{
		 nstat_currentloop_array2 = nstat_array2 - ARRAY2_CHANNELPERLOOP * (nloop_array2-1);
	       }

// STEP 1. for array2
// Copy Original data To GPU for Temporal normalization and Spectral whitening
	       if(cudaMalloc((void **)&dwaveform_array2_orig,sizeof(cufftComplex)*npts_whiten*nstat_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for dwaveform_array2_orig\n");
	       if(cudaMemcpy(dwaveform_array2_orig, (waveform_array2_orig+npts_whiten*iloop_array2*ARRAY2_CHANNELPERLOOP), (nstat_currentloop_array2*npts_whiten*sizeof(cufftComplex)), cudaMemcpyHostToDevice) != cudaSuccess) printf("Error in copying waveform_array2_orig\n");
	       grid_num_array2 = floor(1.0 * nstat_currentloop_array2 * npts_whiten / block_num) + 1;
	       checkNanKernel<<<grid_num_array2,block_num>>>(dwaveform_array2_orig, npts_whiten, nstat_currentloop_array2);
//	       for(j=0;j<10;j++) printf("before onebit : %f\n",waveform_array1_orig[j+npts_whiten*ARRAY1_CHANNELPERLOOP*iloop_array1].x);
//	       printf("\n");
//	       for(j=0;j<10;j++) printf("before onebit : %f\n",waveform_array2_orig[j+npts_whiten*ARRAY2_CHANNELPERLOOP*iloop_array2].x);
// 	       printf("\n");


// STEP 2. for array2
// TEMPORAL NORMALIZATION: onebit or running-absolute-mean
	       if(TEMPORAL_NORMALIZATION == 1){
	         Onebit<<<grid_num_array2,block_num>>>(dwaveform_array2_orig, npts_whiten, nstat_currentloop_array2);
	       }
	       else if(TEMPORAL_NORMALIZATION == 2){
	         npts_halfram = floor(1.0 / freql / delta / 2.0);
	         if(cudaMalloc((void **)&weight_array2,sizeof(float)*npts_whiten*nstat_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for weight_array2\n");
	         RunabsmeanKernel<<<grid_num_array2,block_num>>>(dwaveform_array2_orig, weight_array2, npts_whiten, nstat_currentloop_array2, nt, npts_halfram);
	         cudaSyncCheck(cudaDeviceSynchronize());
	         divweight<<<grid_num_array2,block_num>>>(dwaveform_array2_orig, weight_array2, npts_whiten, nstat_currentloop_array2, nt, npts_halfram);
	         cudaSyncCheck(cudaDeviceSynchronize());
	         cudaFree(weight_array2);
	       }
	       else if(TEMPORAL_NORMALIZATION != 0){
	         printf("Unknow TEMPORAL_NORMALIZATION defined in params.h!\n");
	         exit(0);
	       }

//	       cudaMemcpy((waveform_array1_orig+npts_whiten*iloop_array1*ARRAY1_CHANNELPERLOOP), dwaveform_array1_orig, (npts_whiten*nstat_currentloop_array1*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
//	       for(j=0;j<10;j++) printf("onebit : %f\n",waveform_array1_orig[j+npts_whiten*ARRAY1_CHANNELPERLOOP*iloop_array1].x);
//	       printf("\n");
//	       cudaMemcpy(waveform_array2_orig, dwaveform_array2_orig, (npts_whiten*nstat*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
//	       for(j=0;j<10;j++) printf("onebit : %f\n",waveform_array2_orig[j].x);
//	       printf("\n");


// STEP 3. for array2
// SPECTRAL WHITHENING 
	       if(SPECTRAL_WHITENING){
	         df = (1.0 / delta) / npts_whiten;
	         Nlow = floor((double)freql / (double)df);
	         Nhigh = ceil((double)freqh / (double)df);

	         cufftCreateCheck(cufftPlan1d(&plan_whiten_array2, npts_whiten, CUFFT_C2C, nstat_currentloop_array2));
	         cufftExeInverseCheck(cufftExecC2C(plan_whiten_array2, dwaveform_array2_orig, dwaveform_array2_orig, CUFFT_FORWARD));
	         cudaSyncCheck(cudaDeviceSynchronize());
	         SpectralWhitenKernel<<<grid_num_array2,block_num>>>(dwaveform_array2_orig, Nlow, Nhigh, npts_whiten, nstat_currentloop_array2);
	         cufftExeInverseCheck(cufftExecC2C(plan_whiten_array2, dwaveform_array2_orig, dwaveform_array2_orig, CUFFT_INVERSE));
	         cudaSyncCheck(cudaDeviceSynchronize());
	         cufftDestroy(plan_whiten_array2);
	       }
	       waveform_array2_orig_temp = (cufftComplex *)malloc(sizeof(cufftComplex)*npts_whiten*nstat_currentloop_array2);
	       if(waveform_array2_orig_temp == NULL) printf("Error in mallocing waveform_array2_orig_temp\n");
	       cudaMemcpy(waveform_array2_orig_temp, dwaveform_array2_orig, (npts_whiten*nstat_currentloop_array2*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
	       cudaFree(dwaveform_array2_orig);	

//	       for(j=0;j<10;j++) printf("spectral : %f\n",waveform_array2_orig[j+npts_whiten*ARRAY2_CHANNELPERLOOP*iloop_array2].x);
//	       printf("\n");


// STEP 4. for array 2
// RESHAPE WAVEFORM TO SEGMENTED DATA IN CPU FOR LATER FFT AND XCOR
	       Batch_currentloop_array2 = nstat_currentloop_array2 * seg_number;
	       waveform_array2 = (cufftComplex *)malloc(sizeof(cufftComplex) * npts_fft * Batch_currentloop_array2);
	       for(j=0; j<Batch_currentloop_array2; j++){
	         ichan = floor(1.0 * j / seg_number);
	         iseg = j - ichan * seg_number;
	         index = ichan * npts_whiten + iseg * npts_segment;
	         for(k=0; k<npts_segment; k++){
	           waveform_array2[j*npts_fft+k].x = waveform_array2_orig_temp[index+k].x;
	           waveform_array2[j*npts_fft+k].y = 0.0;
	         }
	         for(k=npts_segment; k<npts_fft; k++){
	           waveform_array2[j*npts_fft+k].x = 0.0;
	           waveform_array2[j*npts_fft+k].y = 0.0;
	         } // end for k
	       } // end for j
//	       for(j=0;j<10;j++) printf("%d copy 2: %f\n",iloop_array2,waveform_array2[j].x);
//	       printf("\n");


// STEP 5. for array2
// COPY TO GPU FOR XCORR
	       Batch_currentloop_array2 = nstat_currentloop_array2 * seg_number;
	       if(cudaMalloc((void **)&dwaveform_array2,sizeof(cufftComplex) * npts_fft * Batch_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for dwaveform_array2\n");
	       cudaMemcpy(dwaveform_array2, waveform_array2, (Batch_currentloop_array2 * npts_fft * sizeof(cufftComplex)), cudaMemcpyHostToDevice);
	       grid_num_array2 = floor(1.0 * Batch_currentloop_array2 * npts_fft / block_num) + 1;
//	       for(ichan_array1=0;ichan_array1<10;ichan_array1++) printf("ARRAY1: %f\n",waveform_array1[ichan_array1+(1-1)*npts_fft].x);
//	       printf("\n");
//	       for(ichan_array2=0;ichan_array2<10;ichan_array2++) printf("ARRAY2: %f\n",waveform_array2[ichan_array2+(1-1)*npts_fft].x);
//	       printf("\n");


// STEP 6. for array2
// CALCULATE L2 NORM FOR EACH SEGMENT FOR LATER NORMALIZED XCORR
	       norm_array2 = (float *)malloc(sizeof(float) * Batch_currentloop_array2);
	       if(norm_array2 == NULL) printf("Error in mallocaing norm_array2\n");
	       if(cudaMalloc((void **)&dnorm_array2,sizeof(float)*Batch_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for dnorm_array2\n");
	       cublasCreateCheck(cublasCreate(&handle));
	       for(j=0; j<Batch_currentloop_array2; j++)  cublasSnrmCheck(cublasSnrm2(handle, npts_segment, (cufftReal *)(dwaveform_array2+j*npts_fft), 2, &norm_array2[j]));
	       cublasDestroyCheck(cublasDestroy(handle));
	       cudaMemcpy(dnorm_array2, norm_array2, (Batch_currentloop_array2*sizeof(float)), cudaMemcpyHostToDevice);
//	       for(j=0; j<10; j++) printf("ARRAY2: %f \n",norm_array2[j]);
//	       printf("\n");
	       free(norm_array2);


// STEP 7. for array2
// FAST FOURIER TRANSFORM FOR EACH SEGMENT
	       cufftCreateCheck(cufftPlan1d(&plan_array2, npts_fft, CUFFT_C2C, Batch_currentloop_array2));
	       cufftExeForwardCheck(cufftExecC2C(plan_array2, dwaveform_array2, dwaveform_array2, CUFFT_FORWARD));
	       cudaSyncCheck(cudaDeviceSynchronize());


// STEP 8.
// LOOP FOR NORMALIZED XCORR AND STACKING
// IN EACH LOOP, ONE CHANNEL IN ARRAY1 IS XCORR WITH ALL CHANNELS IN ARRAY2
	       if(cudaMalloc((void **)&dcc,sizeof(float)*npts_cc*nstat_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for dcc\n");
	       if(cudaMalloc((void **)&dxcorr,sizeof(cufftComplex)*npts_fft*Batch_currentloop_array2) != cudaSuccess) printf("Cuda error: Failed to allocate for dxcorr\n");

//	       if(i==0&&iloop_array1==0&&iloop_array2==0){
//	         for(j=0; j<nstat_array1; j++){
//		 for(k=0; k<nstat_array2*npts_cc; k++) cc_stack[j][k] = 0.0;
//	         }
//	       }
//	       fprintf(stderr,"j = %d k = %d\n",j,k);
//	       for(ichan_array2=0;ichan_array2<10;ichan_array2++) printf("%f\n",cc_stack[nstat-1][(nstat-1)*npts_cc+ichan_array2]);

	       for(ichan_array1=0;ichan_array1<nstat_currentloop_array1; ichan_array1++){
//		 printf("ichan_array1: %d iloop_array1: %d iloop_array2: %d\n",ichan_array1,iloop_array1,iloop_array2);
	         XcorKernel<<<grid_num_array2,block_num>>>(dwaveform_array1, dwaveform_array2, dxcorr, ichan_array1, seg_number, npts_fft, Batch_currentloop_array2, dnorm_array1, dnorm_array2);
	         cufftExeInverseCheck(cufftExecC2C(plan_array2, dxcorr, dxcorr, CUFFT_INVERSE));
	         cudaSyncCheck(cudaDeviceSynchronize());
	         cudaMemcpy(dcc, (*(cc_stack+ichan_array1+iloop_array1*ARRAY1_CHANNELPERLOOP)+iloop_array2*ARRAY2_CHANNELPERLOOP*npts_cc), (nstat_currentloop_array2*npts_cc*sizeof(float)), cudaMemcpyHostToDevice);
	         StackKernel<<<grid_num_array2,block_num>>>(dcc, dxcorr, npts_fft, npts_fft2, npts_lag, npts_cc, seg_number, Batch_currentloop_array2, segyfile_number);
//	         fprintf(stderr,"isegy:%d ichan : %d\n",i,ichan_array1);
	         cudaMemcpy((*(cc_stack+ichan_array1+iloop_array1*ARRAY1_CHANNELPERLOOP)+iloop_array2*ARRAY2_CHANNELPERLOOP*npts_cc), dcc, (nstat_currentloop_array2*npts_cc*sizeof(float)), cudaMemcpyDeviceToHost);
//		 printf("ichan_array1: %d iloop_array1: %d iloop_array2: %d\n",ichan_array1,iloop_array1,iloop_array2);
//		 printf("\n");
	       }

// for debug 
//	   cudaMemcpy(cc, dcc, (nstat*npts_cc*sizeof(float)), cudaMemcpyDeviceToHost);
//	   cudaMemcpy(xcorr, dxcorr, (Batch*npts_fft*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);

//	   cufftExeInverseCheck(cufftExecC2C(plan, dwaveform_array1, dwaveform_array1, CUFFT_INVERSE));
//	   cudaSyncCheck(cudaDeviceSynchronize());
//	   cufftExeInverseCheck(cufftExecC2C(plan, dwaveform_array2, dwaveform_array2, CUFFT_INVERSE));
//	   cudaSyncCheck(cudaDeviceSynchronize());

//	   cudaMemcpy(waveform_array1, dwaveform_array1, (Batch*npts_fft*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
//	   cudaMemcpy(waveform_array2, dwaveform_array2, (Batch*npts_fft*sizeof(cufftComplex)), cudaMemcpyDeviceToHost);
//	   for(ichan_array1=0;ichan_array1<10;ichan_array1++) printf("%f %f\n",waveform_array1[ichan_array1].x/npts_fft,waveform_array1[ichan_array1].y/npts_fft);
//	   printf("\n");
//	   for(ichan_array2=0;ichan_array2<10;ichan_array2++) printf("%f %f\n",waveform_array2[ichan_array2].x/npts_fft,waveform_array2[ichan_array2].y/npts_fft);
//	   printf("\n");
//	  if(i==segyfile_number-1){
//	   for(ichan_array2=1249*seg_number;ichan_array2<(1250*seg_number);ichan_array2++){	
//	      for(k=0;k<npts_cc; k++) printf("%f\n",xcorr[ichan_array2*npts_fft+k+npts_fft2-npts_lag].x);
//	   }
//	   for(ichan_array2=0;ichan_array2<npts_cc;ichan_array2++) printf("%f\n",cc_stack[nstat-1][(nstat-1)*npts_cc+ichan_array2]);}
//	   for(ichan_array2=0;ichan_array2<npts_cc;ichan_array2++) printf("%f\n",cc[(nstat-1)*npts_cc+ichan_array2]);

	       cudaFree(dcc);
	       cudaFree(dxcorr);
	       cufftDestroy(plan_array2);
  	       free(waveform_array2);
	       free(waveform_array2_orig_temp);
	       cudaFree(dwaveform_array2);	
	       cudaFree(dnorm_array2);	
	     } // for iloop_array2
	     free(waveform_array1);
	     cudaFree(dwaveform_array1);	
	     cudaFree(dnorm_array1);	
	   } // for iloop_array1
	} // end for ith segy file

	zap((char *) &segy, 240);
	segy.ntbig = nstat_array2 * npts_cc;
	segy.dt = (short) (delta*1.0e6);
	for(i=0; i<nstat_array1; i++){
	   segy.trseql = i;
	   write(fd_output, &segy, SEGY_HEAD_SIZE);
 	   write(fd_output, *(cc_stack+i), (4*npts_cc*nstat_array2));
	}
	close(fd_output);

//	cufftDestroy(plan);
	free(waveform_array1_orig);
	free(waveform_array2_orig);
//	free(xcorr);
//	free(cc);
	for(i=0; i<nstat_array1; i++) free(cc_stack[i]);
	free(cc_stack);
	free(fl);

}
