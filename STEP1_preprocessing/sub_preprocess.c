//#include "segy_PASSCAL.h"
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <fcntl.h>
//#include <unistd.h>
#include "decim.c"

int countstation(char *filename, int nt){
	int fd;
	int nstat=0;
	char buf[3600];
	struct segy segy;
	float *data;

	if((fd = open(filename,0)) < 0)
	   fprintf(stderr,"cannot open file: %s\n",filename);
	read(fd,buf,3600);
	data = (float *)malloc(sizeof(float)*nt);
	while(read(fd,&segy,SEGY_HEAD_SIZE) == SEGY_HEAD_SIZE){
	   if(segy.ntbig != nt){
		fprintf(stderr, "error: the length of %dth trace is not the same as the first one !\n",nstat);
	  	return -1;}
	   if((read(fd,data,sizeof(float)*nt)) != (sizeof(float)*nt)){
		fprintf(stderr,"read error in %s nt = %d\n",filename,nt);
	 	return -1;}
	   nstat++;
	}
	close(fd);
	free(data);
//	printf("filename: %s nt:%d nstat: %d\n",filename,nt,nstat);	
	return nstat;
}

int countnt(char *filename){
	int fd;
	int nt;
	char buf[3600];
	struct segy segy;

	if((fd = open(filename,0)) < 0)
           fprintf(stderr,"cannot open file: %s\n",filename);
        read(fd,buf,3600);
	read(fd,&segy,SEGY_HEAD_SIZE);
	nt = segy.ntbig;

	close(fd);
	return nt;
}

double countdelta(char *filename){
        int fd;
        double delta;
        char buf[3600];
        struct segy segy;

        if((fd = open(filename,0)) < 0)
           fprintf(stderr,"cannot open file: %s\n",filename);
        read(fd,buf,3600);
        read(fd,&segy,SEGY_HEAD_SIZE);
        delta = (double)(segy.dt/1e6);

        close(fd);
        return delta;
}


void diff(float *data,float *data_dif, int nt, double delta){
	int i;
	data_dif[0] = (data[1] - data[0]) / delta;
	for(i=1; i<(nt-1); i++) data_dif[i] = (data[i+1]-data[i-1])/(2.0*delta);
	data_dif[i] = (data[i] - data[i-1]) / delta;
}

void    rtrend1(float *y, int n) {
     int i;
     double a, b, a11, a12, a22, y1, y2;
     y1 = y2 = 0.;
     for(i=0;i<n;i++) {
       y1 += i*y[i];
       y2 += y[i];
     }
     a12 = 0.5*n*(n-1);
     a11 = a12*(2*n-1)/3.;
     a22 = n;
     b = a11*a22-a12*a12;
     a = (a22*y1-a12*y2)/b;
     b = (a11*y2-a12*y1)/b;
     for(i=0;i<n;i++) {
       y[i] = y[i] - a*i - b;
     }
}

void rmean1(float *x, int n){
        float ave;
        int i;
        ave = 0.0;
        for(i=0;i<n;i++){
                ave += x[i];
        }
        ave = ave/n;
        for(i=0;i<n;i++){
                x[i] -= ave;
        }
}

void decimate(float *data, float *data_new, int newlength, int decimatefactor){
	int i;

	for(i=0; i<newlength; i++) data_new[i] = data[i*decimatefactor];
}

void decimate_filter(float *data, float *data_new, int length, int newlength, int decimatefactor){
	char filename[5];
	int  nb, nchp1;
	int  i, nloop;
	int  residual;
	float dtemp;	
	float *c, *buf;

	FILE *fp;
	
	sprintf(filename,"dec%d",decimatefactor);
	fp = fopen(filename,"r");
	fscanf(fp,"%f %d",&dtemp,&nchp1);
	nb = nchp1 * 2 + 100;
	buf = (float *) malloc(sizeof(float) * nb);
	c = (float *) malloc(sizeof(float) * nchp1);
	nloop = floor(nchp1 / 5);
	for(i=0; i<nloop; i++){
	   fscanf(fp,"%E%E%E%E%E",(c+i*5),(c+i*5+1),(c+i*5+2),(c+5*i+3),(c+5*i+4));
	}
	residual = nchp1 - nloop * 5;
	if(residual){
	   switch(residual){
	      case 1:
		fscanf(fp,"%E",(c+nloop*5));
		break;
	      case 2:
		fscanf(fp,"%E%E",(c+nloop*5),(c+nloop*5+1));
		break;
	      case 3:
		fscanf(fp,"%E%E%E",(c+nloop*5),(c+nloop*5+1),(c+nloop*5+2));
		break;
	      case 4:
		fscanf(fp,"%E%E%E%E",(c+nloop*5),(c+nloop*5+1),(c+nloop*5+2),(c+nloop*5+3));
		break;
	   }
	}
	fclose(fp);
	decim(data, length, buf, nb, c, nchp1, 1, decimatefactor, data_new, &newlength);
	free(c);
	free(buf);
}
void zap(char *x, int n)
   {
        int i;
        for(i=0; i<n; i++) x[i]= 0;
   }

