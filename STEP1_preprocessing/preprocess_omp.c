/*******************************************************************
 *	main
 *	task: gpu based template matching
 *
 *	input: 
 *		directory path for segy files	argv[1]
 *	Usage:	
 *		./preprocess [directory]
*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>

#include "segy_PASSCAL.h"
#include "sub_preprocess.c"
#include "../params.h"

/******************************* main **************************************/
int main(int argc, char *argv[])
{
	clock_t begin, end;
	double time_spent;
	struct filelist *fl;
	char  filename_list[256];
	int   segyfile_number = 0;
	int   i;
	int   threads_number;
	FILE  *fp;

	long int order = 4;
	char type[4] = "BP";
	char proto[4] = "BU";
	begin = clock();
//---------------------------------------------------------------------------
//read input paramters
	if((fp=fopen(SEGY_LIST,"r")) == NULL) perror(SEGY_LIST);
	while(fscanf(fp, "%s", filename_list) != EOF) segyfile_number += 1;

	rewind(fp);
	fl = (struct filelist *)malloc(sizeof(struct filelist) * segyfile_number);
	for(i=0;i<segyfile_number;i++){
	   fscanf(fp, "%s", filename_list);
	   strcpy(fl[i].file, filename_list);
//	   fprintf(stderr,"%s\n",fl[i].file);
	}
	printf("INPUT: %d hour segy files\n", segyfile_number);
	fclose(fp);

	threads_number = segyfile_number >= NUM_THREADS ? NUM_THREADS : segyfile_number;

//---------------------------------------------------------------------------
//count das station numbers

	if(mkdir(ARRAY1_DOWNSAMPLE_DIR,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) != 0) fprintf(stderr,"%s already exist or Creating Wrong !\n",ARRAY1_DOWNSAMPLE_DIR);
	if(mkdir(ARRAY2_DOWNSAMPLE_DIR,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) != 0) fprintf(stderr,"%s already exist or Creating Wrong !\n",ARRAY2_DOWNSAMPLE_DIR);

	#pragma omp parallel for schedule(dynamic,1) num_threads(threads_number)
	for(i=0; i<segyfile_number; i++){
	   struct segy segy_array1, segy_array2;
	   char  filename[256];
	   char  buf[3600];
	   char  outfilename[256];
	   double delta_array1, delta_array2;
	   int   nt_array1, nt_array2;	
	   int   ntnew_array1, ntnew_array2;
	   int   nstat_array1, nstat_array2;
	   int   fdin_array1, fdin_array2;
	   int   fdout_array1, fdout_array2;
	   float *data_array1, *data_array2;
	   float *data_process, *data_new;

	   long int nsects1;
	   float p_sn[30], p_sd[30];

	   // preprocessing on array_1
	   sprintf(filename, "%s/%s",ARRAY1_ORIGIN_DIR,fl[i].file);
	   if((fdin_array1=open(filename,0)) < 0){perror(filename); exit(0);} 
	   sprintf(outfilename, "%s/%s",ARRAY1_DOWNSAMPLE_DIR,fl[i].file);
	   if((fdout_array1 = creat(outfilename,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0){perror(outfilename); exit(0);}  
	   read(fdin_array1, buf, 3600);
	   write(fdout_array1, buf, 3600);
	   nstat_array1 = 0;
	   while(read(fdin_array1, &segy_array1, SEGY_HEAD_SIZE) == SEGY_HEAD_SIZE){
	      nt_array1 = segy_array1.ntbig;
	      ntnew_array1 = floor((nt_array1-1)/ARRAY1_DECIMATE) + 1;
	      delta_array1 = (double)(segy_array1.dt/1e6);
	      data_array1 = (float *)malloc(sizeof(float) * nt_array1);
	      data_process = (float *)malloc(sizeof(float) * nt_array1);
	      data_new = (float *) malloc(sizeof(float) * ntnew_array1);
	      read(fdin_array1,data_array1,sizeof(float) * nt_array1);    // read data

//	      diff(data_array1, data_process, nt_array1, delta_array1);   // differentiate data
	      rtrend1(data_process,nt_array1);			    // remove trend
	      rmean1(data_process,nt_array1);			    // remove mean
//            design(order, type, proto, 1., 1., (double) freql, (double) freqh, (double) delta_array1, p_sn, p_sd, &nsects1);	
//	      apply(data_process, (long int) nt_array1, 1, p_sn, p_sd, nsects1);  // band pass filter
	      decimate_filter(data_process,data_new,nt_array1,ntnew_array1,ARRAY1_DECIMATE); // decimate data
	      segy_array1.dt = (short) (delta_array1 * ARRAY1_DECIMATE * 1e6);
	      segy_array1.ntbig = ntnew_array1;
	      write(fdout_array1, &segy_array1, SEGY_HEAD_SIZE);
	      write(fdout_array1, data_new, 4*ntnew_array1);
	      free(data_array1);
	      free(data_process);
	      free(data_new);
	      nstat_array1++;
	   }
	   printf("%s chan: %d dt=%f %f %d\n",fl[i].file, nstat_array1, delta_array1*ARRAY1_DECIMATE, (segy_array1.dt/1e6), segy_array1.ntbig);
	   close(fdin_array1);
	   close(fdout_array1);

//	   sprintf(filename, "%s/%s",ARRAY2_ORIGIN_DIR,fl[i].file);
//	   if((fdin_array2=open(filename,0)) < 0){perror(filename); exit(0);}
//	   sprintf(outfilename, "%s/%s",ARRAY2_DOWNSAMPLE_DIR,fl[i].file);
//	   if((fdout_array2 = creat(outfilename,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0){perror(outfilename); exit(0);}  
//	   read(fdin_array2, buf, 3600);
//	   write(fdout_array2, buf, 3600);
//	   nstat_array2 = 0;
//	   while(read(fdin_array2,&segy_array2,SEGY_HEAD_SIZE) == SEGY_HEAD_SIZE){
//	      nt_array2 = segy_array2.ntbig;
//	      ntnew_array2 = floor((nt_array2 - 1)/ARRAY2_DECIMATE) + 1;
//	      delta_array2 = (double)(segy_array2.dt/1e6);
//	      data_array2 = (float *)malloc(sizeof(float) * nt_array2);
//	      data_process = (float *)malloc(sizeof(float) * nt_array2);
//	      data_new = (float *) malloc(sizeof(float) * ntnew_array2);
//	      read(fdin_array2,data_array2,sizeof(float) * nt_array2);    // read data

//	      diff(data_array2, data_process, nt_array2, delta_array2);   // differentiate data
//	      rtrend1(data_process, nt_array2);			    // remove trend
//	      rmean1(data_process, nt_array2);			    // remove mean
//	      design(order, type, proto, 1., 1., (double) freql, (double) freqh, (double) delta_array2, p_sn, p_sd, &nsects1);	
//	      apply(data_process, (long int) nt_array2, 1, p_sn, p_sd, nsects1);  // band pass filter
//	      decimate_filter(data_process, data_new, nt_array2, ntnew_array2, ARRAY2_DECIMATE);
//	      segy_array2.dt = (short) (delta_array2 * ARRAY2_DECIMATE * 1e6);
//	      segy_array2.ntbig = ntnew_array2;
//	      write(fdout_array2, &segy_array2, SEGY_HEAD_SIZE);
//	      write(fdout_array2, data_new, 4*ntnew_array2);
//	      free(data_array2);
//	      free(data_process);
//	      free(data_new);
//	      nstat_array2++;
//	   }
//	   fprintf(stderr,"2 chan: %d dt=%f %f %d\n",nstat_array2, delta_array2*ARRAY2_DECIMATE, (segy_array2.dt/1e6), segy_array2.ntbig);
//	   close(fdin_array2);
//	   close(fdout_array2);
//	   if(fabs((delta_array1*ARRAY1_DECIMATE)-(delta_array2*ARRAY2_DECIMATE))>1e-7)
//	      fprintf(stderr,"WARNING: DIFFERENT DELTA BETWEEN TWO DOWNSAMPLED SEGY FILES FOR HOUR: %s! Cannot go to xcorrelation!\n",fl[i].file);
	}
//	if((fp=fopen(SEGY_INFO,"w+")) == NULL) perror(SEGY_INFO);
//	fprintf(fp, "#DOWNSAMPLED SEGY INFO SUMMARY:  nstat  nt  delta\n", nstat_array1, segy_array1.ntbig, delta_array1*ARRAY1_DECIMATE);
//	fprintf(fp, "ARRAY1: %d %d %.3f\n", nstat_array1, segy_array1.ntbig, delta_array1*ARRAY1_DECIMATE);
//	fprintf(fp, "ARRAY2: %d %d %.3f\n", nstat_array2, segy_array2.ntbig, delta_array2*ARRAY2_DECIMATE);
//	fclose(fp);
	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
	fprintf(stderr,"Preprocessing Done! Cost: %f s per threads\n",time_spent/threads_number);
}	
