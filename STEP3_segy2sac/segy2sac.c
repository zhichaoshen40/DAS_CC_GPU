#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "segy_PASSCAL.h"
#include "sac.h"
#include "sacio.c"
#include "../params.h"
#include <unistd.h>

int main(int argc,char *argv[]){
        struct filelist fl[1000];
	DIR *dirp;
	struct dirent *dp;
	int i,j,nfl;
	int fd;
	int nstat_array1, nstat_array2;
	int npts_segment;
	int nt,ichan;
	double dt;
	char filename[256], outname[256];
	char outdir[256];
	char buf[3600];
	struct segy segy;
	float *data, *data_copy;
	int flag =1;
	FILE *fp;

	if((dirp = opendir(XCORR_OUTPUT_DIR)) == NULL){
	   fprintf(stdout,"cannot open dirctory: %s\n",argv[1]);
	   exit(0);
	}
	nfl = 0;
	while((dp = readdir(dirp)) != NULL){
	   if(dp->d_name[0] == '.') continue;
	   if(strncmp((dp->d_name),"stack.segy",10) != 0) continue;
	   strcpy(fl[nfl].file,dp->d_name);
	   nfl++;}
	closedir(dirp);

        for(i=0;i<nfl;i++){
	   sprintf(filename,"%s/%s",XCORR_OUTPUT_DIR,fl[i].file);
	   if((fd = open(filename,0)) < 0){
	      fprintf(stdout,"cannot open file: %s\n",filename);
	      exit(0);
	   }
	   read(fd,buf,3600);
	   while(read(fd,&segy,240) == 240){
	      nt = segy.ntbig;
	      ichan = segy.trseql;
	      dt = (double) (segy.dt/1e6);
	      npts_segment = 2 * floor(MAX_LAG / dt) + 1;
	      nstat_array2 = nt / npts_segment;
	      data=(float *)malloc(nt*sizeof(float));
	      fprintf(stdout,"ichan=%4d nt= %d dt=%7.5f npts_segment=%d nstat_array2=%d\n",ichan,nt,dt,npts_segment,nstat_array2);
	      if((read(fd,data,sizeof(float)*nt)) != (sizeof(float)*nt)){
		fprintf(stdout,"read error in %s nt= %d\n",fl[i].file,nt);
		exit(0);
	      }
	      sprintf(outdir,"%s/%04d",XCORR_OUTPUT_DIR,ichan);
	      if(mkdir(outdir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) != 0) printf("\t%s already exist or Creating Wrong !\n",outdir);
	      for(j=0;j<nstat_array2;j++){
	         sprintf(outname,"%s/ccA%04dB%04d.sac",outdir,ichan,j);
	         wrtsacfromsegy_az(outname, dt, npts_segment, (-1.0*MAX_LAG), (float) (ichan+1), (float) (j+1),(data+j*npts_segment));
	      }
	      free(data);
//	      fprintf(stdout,"%04d %f %f\n",ichan,dt,(float)dt);
//	      printf("%s\n",outname);
//	      if(nstat == 1249){
//		fp=fopen("test","w+");
//		for(j=0;j<nt;j++) fprintf(fp,"%f\n",data[j+index]);
//		fclose(fp);}
	   }
	   close(fd);
	}
}

swap_segy(struct segy *segy)
   {
        char *psegy;
        psegy= (char *)(segy);
        swap4(&psegy[  0],7);
        swap2(&psegy[ 28],4);
        swap4(&psegy[ 36],8);
        swap2(&psegy[ 68],2);
        swap4(&psegy[ 72],4);
        swap2(&psegy[ 88],45);
        swap4(&psegy[180],5);
   }


struct b2 { char b1, b2; };
swap2(struct b2 *x, int n)
   {
        int i;
        char t;
        for(i=0; i<n; i++)
           {
                t= x[i].b1;
                x[i].b1= x[i].b2;
                x[i].b2= t;
           }
   }

struct b4 { char b1, b2, b3, b4; };
swap4(struct b4 *x, int n)
   {
        int i;
        char t;
        for(i=0; i<n; i++)
           {
                t= x[i].b1;
                x[i].b1= x[i].b4;
                x[i].b4= t;
                t= x[i].b2;
                x[i].b2= x[i].b3;
                x[i].b3= t;
           }
   }

