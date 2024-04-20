__host__ double make_time(int yr,int mon,int day,int hr,int min,double sec){
        struct tm *tm;
        time_t iut;
        int isec;
        double dsec, ut;

        putenv("TZ=GMT");
        isec= (int)(sec);
        dsec= sec - (double)(isec);
	tm = (struct tm*)malloc(sizeof(struct tm));
        tm->tm_year      =yr - 1900;
        tm->tm_mon       =mon -1;
        tm->tm_mday      =day;
        tm->tm_hour      =hr;
        tm->tm_min       =min;
        tm->tm_sec       =isec;
        iut = mktime(tm);
        ut= (double)(iut) + dsec;
        return(ut);
}

__host__ int get_int(char *str,int n){
        char c;
        int val;
        c= str[n];
        str[n]= '\0';
        val= atol(str);
        str[n]= c;
        return(val);
}

__host__ void convert_datetime2ut(char *timestr,double *ut) {
        int year, mon, day, hour, min;
        float sec= 0.0;

        year =  get_int(&timestr[0],4);
        mon = get_int(&timestr[4],2);
        day = get_int(&timestr[6],2);
        hour = get_int(&timestr[8],2);
        min = 0.0;

        *ut = make_time(year,mon,day,hour,min,sec);
}

__host__ int filecomp (const void * a, const void * b){
        struct filelist *ia = (struct filelist *)a;
        struct filelist *ib = (struct filelist *)b;

        return ia->ut >= ib->ut ? 1 : -1;
}


__host__ void checksegyfiles(char *filename, double delta, int nt, int NSTA){
        int fd;
        int nstat=0;
        char buf[3600];
        struct segy segy;
        float *data;
	double delta_chan;
	
        if((fd = open(filename,0)) < 0){
           printf("cannot open file: %s\n",filename);
	   exit(0);
	}
        read(fd,buf,3600);
	if(read(fd,&segy,SEGY_HEAD_SIZE) != SEGY_HEAD_SIZE){
	   printf("Fail in reading SEGY_HEAD\n");
	   exit(0);
	}

        if(segy.ntbig != nt){
           printf("error: the npts %d in file %s channel %d is different from nt: %d !\n",segy.ntbig,filename,nstat,nt);
           exit(0);
	   }

	delta_chan = (double)(segy.dt/1e6);
	if(fabs(delta_chan - delta) > 1e-5){
           printf("error: the delta %f in file %s channel %d is not the same as input delta %f\n",delta_chan,filename,nstat,delta);
           exit(0);
	}

        data = (float *)malloc(sizeof(float) * nt);
	if((read(fd,data,sizeof(float)*nt)) != (sizeof(float)*nt)){
           printf("read error in %s nt = %d\n",filename,nt);
           exit(0);
	}
	nstat++;
        while(read(fd,&segy,SEGY_HEAD_SIZE) == SEGY_HEAD_SIZE){
           if(segy.ntbig != nt){
              printf("error: the npts %d in file %s channel %d is different from nt: %d !\n",segy.ntbig,filename,nstat,nt);
              exit(0);
	   }
	   delta_chan = (double)(segy.dt/1e6);
	   if(fabs(delta_chan - delta) > 1e-5){
              printf("error: the delta %f in file %s channel %d is not the same as input delta %f\n",delta_chan,filename,nstat,delta);
	      exit(0);
	   }
           if((read(fd,data,sizeof(float)*nt)) != (sizeof(float)*nt)){
              fprintf(stderr,"read error in %s nt = %d\n",filename,nt);
              exit(0);
	      }
           nstat++;
        }
        close(fd);
        free(data);
//      printf("filename: %s nt:%d nstat: %d\n",filename,nt,nstat);     

        if(nstat<NSTA){
	   printf("Error: number of stations %d in file %s is less than NSTA %d\n",nstat,filename,NSTA);
	   exit(0);
	}
	else if(nstat>NSTA){
	   printf("Warning: number of stations %d in file %s is more than NSTA %d, only the first %d channels will be used\n",nstat,filename,NSTA,NSTA);
	}
}

__host__ void cublasCreateCheck(cublasStatus_t status){
	if(status != CUBLAS_STATUS_SUCCESS){
	   printf("CUBLAS error: Initialization failed! error code: %d\n",status);
	   exit(0);}
}

__host__ void cublasSnrmCheck(cublasStatus_t status){
	if(status != CUBLAS_STATUS_SUCCESS){
	   printf("CUBLAS error: NORM2 CALCULATION failed! error code: %d\n",status);
	   exit(0);}
}

__host__ void cublasDestroyCheck(cublasStatus_t status){
	if(status != CUBLAS_STATUS_SUCCESS){
	   printf("CUBLAS error: SHUTDOWN failed! error code: %d\n",status);
	   exit(0);}
}

__host__ void cufftCreateCheck(cufftResult_t Result){
	if(Result != CUFFT_SUCCESS){
           printf("CUFFT ERROR: plan creation failed! error code: %d\n",Result);
	   exit(0);}
}

__host__ void cufftExeForwardCheck(cufftResult_t Result){
	if(Result != CUFFT_SUCCESS){
           printf("CUFFT ERROR: ExecC2C Forward failed! error code: %d\n",Result);
	   exit(0);}
}

__host__ void cufftExeInverseCheck(cufftResult_t Result){
	if(Result != CUFFT_SUCCESS){
           printf("CUFFT ERROR: ExecC2C Inverse failed! error code: %d\n",Result);
	   exit(0);}
}

__host__ void cudaSyncCheck(cudaError_t Result){
	if(Result != cudaSuccess){
           printf("CUDA ERROR: C2C failed to synchronize! error code: %d\n",Result);
	   exit(0);}
}

__host__ void zap(char *x, int n){
	int i;
        for(i=0; i<n; i++) x[i]= 0;
}
