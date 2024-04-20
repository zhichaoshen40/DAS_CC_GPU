#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "segy_PASSCAL.h"
#include "sac.h"
#include "sacio.c"
#include "../params.h"
#include <unistd.h>

void zap(char *x, int n){
        int i;
        for(i=0; i<n; i++) x[i]= 0;
}

int main(int argc,char *argv[]){
        struct filelist *fl;
	int i,j,k;
	int nfl = 0;
	int fdin;
	int fdout;
	int npts_cc;
	float sum;
	float *weight;
	float **cc_stack;
	int nstat_array1, nstat_array2;
	int nt_array1, nt_array2;
	double delta_array1, delta_array2;
	double dt;
	char filename[256], outname[256];
	char buf[3600];
	struct segy segy;
	float *cc_temp;
	int flag =1;
	FILE *fp;

	if((fp=fopen(argv[1],"r")) == NULL) perror(argv[1]);
        while(fscanf(fp, "%s %f", filename,&sum) != EOF) nfl += 1;
	rewind(fp);
	sum = 0.0;
	printf("nfl %d\n",nfl);
	fl = (struct filelist *) malloc(nfl*sizeof(struct filelist));
	weight = (float *) malloc(nfl*sizeof(float));
	for(i=0;i<nfl;i++){
           fscanf(fp, "%s %f", fl[i].file,&weight[i]);
	   sum += weight[i];
           printf("%d: %s %f %f\n",i, fl[i].file,weight[i],sum);
	}
	fclose(fp);

	printf("  SEGY INFO from %s\n", SEGY_INFO);
        if((fp=fopen(SEGY_INFO,"r")) == NULL) perror(SEGY_INFO);
        while(fgetc(fp) != '\n') i++;
        fscanf(fp, "%s%d%d%lf", buf, &nstat_array1, &nt_array1, &delta_array1);
        fscanf(fp, "%s%d%d%lf", buf, &nstat_array2, &nt_array2, &delta_array2);
        printf("\t        nstat   nt    delta\n");
        printf("\tArray1: %5d %d %.5f\n",nstat_array1, nt_array1, delta_array1);
        printf("\tArray2: %5d %d %.5f\n\n",nstat_array2, nt_array2, delta_array2);
        fclose(fp);

        npts_cc = 2 * floor(MAX_LAG / delta_array1) + 1;
	printf("npts_cc %d\n",npts_cc);
	if(mkdir(STACK_XCORR_OUTPUT_DIR,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) != 0) printf("\t%s already exist or Creating Wrong !\n",STACK_XCORR_OUTPUT_DIR);
	sprintf(outname, "%s/stack.segy", STACK_XCORR_OUTPUT_DIR);
	if((fdout = creat(outname,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0){
           printf("cannot creat %s\n",outname);
           exit(0);}

	cc_temp = (float *)malloc(sizeof(float) * nstat_array2 * npts_cc);	
	cc_stack = (float **) malloc(sizeof(float *) * nstat_array1);
        for(i=0; i<nstat_array1; i++){
           cc_stack[i] = (float *)malloc(sizeof(float) * nstat_array2 * npts_cc);
           if(cc_stack[i] == NULL) printf("Error in ccstack %d\n",i);
           for(j=0; j<(nstat_array2*npts_cc); j++) cc_stack[i][j] = 0.0;
        }

	for(i=0; i<nfl; i++){
	   printf("reading %dth file\n",i);
	   if((fdin = open(fl[i].file,0)) < 0){
              printf("cannot open file: %s\n",fl[i].file);
              exit(0);
           }
           read(fdin, buf, 3600);
           if(i==0) write(fdout, buf, 3600);
           for(j=0; j<nstat_array1; j++){
	     if(read(fdin,&segy,SEGY_HEAD_SIZE) != SEGY_HEAD_SIZE){
                printf("Fail in reading SEGY_HEAD\n");
                exit(0);}
	     if(segy.ntbig != (nstat_array2*npts_cc)){
                printf("error: the npts %d in file %s channel %d is different from nt: %d !\n",segy.ntbig,fl[i].file,j,nstat_array2*npts_cc);
                exit(0);}
	     dt = (double)(segy.dt/1e6);
             if(fabs(dt - delta_array1) > 1e-5){
                printf("error: the delta %f in file %s channel %d is not the same as input delta %f\n",dt,fl[i].file,j,delta_array1);
                exit(0);}

             read(fdin, cc_temp, sizeof(float)*(nstat_array2*npts_cc));
	     printf("adding %dth file %dth chan\n",i,j);
	     for(k=0;k<(nstat_array2*npts_cc);k++) cc_stack[j][k] = cc_stack[j][k] + cc_temp[k] * weight[i] / sum;
	   }
	   close(fdin);
	}

        zap((char *) &segy, 240);
        segy.ntbig = nstat_array2 * npts_cc;
        segy.dt = (short) (delta_array1*1.0e6);
        for(i=0; i<nstat_array1; i++){
	   printf("writing %dth chan: npts: %d delta %f\n",i,segy.ntbig,(float)segy.dt/1e6);
           segy.trseql = i;
           write(fdout, &segy, SEGY_HEAD_SIZE);
           write(fdout, *(cc_stack+i), (4*npts_cc*nstat_array2));
        }

	close(fdout);
	free(fl);
	free(cc_temp);
	for(i=0; i<nstat_array1; i++) free(cc_stack[i]);
        free(cc_stack);
}

