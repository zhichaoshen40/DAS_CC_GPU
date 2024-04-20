#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "segy_PASSCAL.h"
#include "../params.h"

int main(int argc, char *argv[]){
        struct segy segy;
        struct filelist *fl;
	int i, j, nstat;
	int fd;
	int npts=180000;
	double delta;
	double dt=0.02;
	int segyfile_number = 0;
	char filename[365];
	char buf[3600];
	float *data;
	FILE *fp;
	
	if((fp=fopen(argv[1],"r")) == NULL) perror(argv[1]);
	while(fscanf(fp, "%s", filename) != EOF) segyfile_number += 1;
        rewind(fp);
        fl = (struct filelist *)malloc(sizeof(struct filelist) * segyfile_number);
        for(i=0;i<segyfile_number;i++){
           fscanf(fp, "%s", filename);
           strcpy(fl[i].file, filename);
//           convert_datetime2ut(fl[i].file,&fl[i].ut);
        }
        printf("  INPUT: %d HOUR SEGY FILES (SORTED BY TIME)\n", segyfile_number);
        fclose(fp);

	for(i=0; i<segyfile_number; i++){
	   if((fd = open(fl[i].file,0)) < 0){
             printf("cannot open file: %s\n",fl[i].file);
             exit(0);
           }
           read(fd,buf,3600);
	   nstat = 0;
	   data=(float *)malloc(npts*sizeof(float));
           while(read(fd,&segy,240) == 240){
	     if(segy.ntbig != npts) printf("NPTS ERROR IN FILE %s CHAN %d\n",fl[i].file,nstat);
	     delta = (double) (segy.dt/1e6);
	     if(fabs(delta-dt) > 1e-7) printf("DELTA ERROR IN FILE %s CHAN %d\n",fl[i].file,nstat);
	     read(fd,data,sizeof(float)*segy.ntbig);
	     for(j=0;j<npts;j++){
	       if(isnan(data[j])){
		 printf("NAN ERROR IN %d FILE %s CHAN %d\n",i,fl[i].file,nstat);
		 break;
	       }
	     }
	     nstat++;
	   }
	   if(nstat!=1250) printf("NSTA ERROR IN FILE %s\n", fl[i].file);
	   printf("%dth: %s done!\n",i,fl[i].file);
	}
	free(fl);
}
