__global__ void checkNanKernel(cufftComplex *dwaveform_orig, int npts_whiten, int nstat){
        int i = blockIdx.x * blockDim.x + threadIdx.x;

	if(i < (nstat*npts_whiten)){
	  if(isnan(dwaveform_orig[i].x)){
	     dwaveform_orig[i].x = 0.0;
	     dwaveform_orig[i].y = 0.0;
	  }
	}
}

__global__ void Onebit(cufftComplex *dwaveform_orig, int npts_whiten, int nstat){
        int i = blockIdx.x * blockDim.x + threadIdx.x;

	if(i < (nstat*npts_whiten)){
	   if(dwaveform_orig[i].x > 0.0){
	      dwaveform_orig[i].x = 1.0;
	   }
	   else if(dwaveform_orig[i].x < 0.0){
	      dwaveform_orig[i].x = -1.0;
	   }
	   else{ 
	      dwaveform_orig[i].x = 0.0;
	   }
	}
        __syncthreads();
}

__global__ void RunabsmeanKernel(cufftComplex *dwaveform_orig, float *weight, int npts_whiten, int nstat, int nt, int npts_halfram){
        int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, index, ichan;

	ichan = floor(1.0 * i / npts_whiten);
	index = i - ichan * npts_whiten;

	if(i<(npts_whiten*nstat)) {
	   weight[i] = 0.0;
	   if(index<npts_halfram && index>=0){
	      for(j=0; j<=(npts_halfram+index); j++){ weight[i] += fabs(dwaveform_orig[i+j-index].x); }
	      weight[i] = weight[i] / (index + npts_halfram + 1);
	   }
	   else if(index>=npts_halfram && index<(nt-npts_halfram)){
	      for(j=0; j<=(2*npts_halfram); j++){ weight[i] += fabs(dwaveform_orig[i+j-npts_halfram].x); }
	      weight[i] = weight[i] / (2 * npts_halfram + 1);
	   }
	   else if(index>=(nt-npts_halfram) && index<nt){
	      for(j=0; j<(nt-index+npts_halfram); j++){ weight[i] += fabs(dwaveform_orig[i+j-npts_halfram].x);}
	      weight[i] = weight[i] / (nt - index + npts_halfram);
	   }
	}
        __syncthreads();
}


__global__ void divweight(cufftComplex *dwaveform_orig, float *weight, int npts_whiten, int nstat, int nt, int npts_halfram){
        int i = blockIdx.x * blockDim.x + threadIdx.x;
	int index, ichan;

	ichan = floor(1.0 * i / npts_whiten);
	index = i - ichan * npts_whiten;

	if(i<(npts_whiten*nstat) && (index<nt)) {
	   if((index<(2*npts_halfram)) || (index>=(nt-2*npts_halfram))){
		dwaveform_orig[i].x = 0.0;
	   }
	   else if(fabs(weight[i]) > 1e-6){
		dwaveform_orig[i].x = dwaveform_orig[i].x / weight[i];
	   }
	   else{
		dwaveform_orig[i].x = 0.0;
	   }
	}
        __syncthreads();
}

__global__ void SpectralWhitenKernel(cufftComplex *dwaveform_orig, int Nlow, int Nhigh, int npts_whiten, int nstat){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int ichan, index;
	double phi;
	
	ichan = floor(1.0 * i / npts_whiten);
	index = i - ichan * npts_whiten;

	if(i<(nstat * npts_whiten)){
	   if(index >= Nlow && index <=Nhigh){
	      phi = atan2(dwaveform_orig[i].y, dwaveform_orig[i].x);
	      dwaveform_orig[i].x = 1.0 * cos(phi) / npts_whiten;
	      dwaveform_orig[i].y = 1.0 * sin(phi) / npts_whiten;
	   }
	   else{
	      dwaveform_orig[i].x = 0.0;
	      dwaveform_orig[i].y = 0.0;
	   }
	}
        __syncthreads();
}


//one segment in array1 cross-correlated with all segments in array2
__global__ void XcorKernel(cufftComplex *dwaveform_array1, cufftComplex *dwaveform_array2, cufftComplex *dxcorr, int ichan_array1, int seg_number, int npts_fft, int Batch, float *dnorm_array1, float *dnorm_array2){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int ichan_array2, sign;
	int index_array1, index_array2;
	int chan_index, iseg_array2;

	float factor;

	ichan_array2 = floor(1.0 * i / (npts_fft * seg_number));
	chan_index = i - ichan_array2 * npts_fft * seg_number;
	iseg_array2 = floor(1.0 * chan_index / npts_fft);
	index_array2 = chan_index - iseg_array2 * npts_fft;
	index_array1 = ichan_array1 * npts_fft * seg_number + chan_index;
	factor = npts_fft * dnorm_array1[ichan_array1*seg_number+iseg_array2] * dnorm_array2[iseg_array2+ichan_array2*seg_number];
	
	if(i < Batch*npts_fft){
	   sign = pow(-1.0, 1.0 * index_array2);
	   if(index_array2==0){
	      dxcorr[i].x = dwaveform_array1[index_array1].x * dwaveform_array2[i].x / factor;
	      dxcorr[i].y = -1.0 * dwaveform_array1[index_array1].y * dwaveform_array2[i].y / factor;
	   }
	   else{
	      dxcorr[i].x = sign * (dwaveform_array1[index_array1].x * dwaveform_array2[i].x + dwaveform_array1[index_array1].y * dwaveform_array2[i].y) / factor;
	      dxcorr[i].y = sign * (dwaveform_array1[index_array1].y * dwaveform_array2[i].x - dwaveform_array1[index_array1].x * dwaveform_array2[i].y) / factor;
	   }
	   if(fabs(factor)<1e-6){
	      dxcorr[i].x = 0.0;
	      dxcorr[i].y = 0.0;
	   }
	}
        __syncthreads();
}


__global__ void StackKernel(float *dstackcc, cufftComplex *dxcorr, int npts_fft, int npts_fft2, int npts_lag, int npts_cc, int seg_number, int Batch, int segyfile_number){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, ichan;	
	int chan_index;
	int cc_lag;

	ichan = floor(1.0 * i / (npts_fft * seg_number));
	chan_index = i - ichan * npts_fft * seg_number;
	cc_lag = npts_fft2 - npts_lag;

	if((i<(npts_fft*Batch)) && (chan_index<(npts_cc+cc_lag)) && (chan_index>=cc_lag)){
//	  dstackcc[ichan * npts_cc + chan_index - cc_lag] = 0.0;
	  for(j=0;j<seg_number;j++) dstackcc[ichan * npts_cc + chan_index - cc_lag] += dxcorr[j * npts_fft + i].x / (seg_number * segyfile_number);
//	  dstackcc[ichan * npts_cc + chan_index - cc_lag] /= seg_number;
	}
        __syncthreads();
}
