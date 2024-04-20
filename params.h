// define data path here; 
// just need to change "ARRAY1/ARRAY2_DOWNSAMPLE_DIR" below if you already have downsampled data
// make sure that both continuous data have same sampling rate

#define ARRAY1_ORIGIN_DIR "/kuafu/zshen/Benchmark/1hoursegy"
#define ARRAY2_ORIGIN_DIR "/net/qin/data3-nobackup/RidgecrestDAS/DASarchive/Ridgecrest"
#define ARRAY1_DOWNSAMPLE_DIR "/kuafu/zshen/Benchmark/downsample"
#define ARRAY2_DOWNSAMPLE_DIR "/kuafu/zshen/Benchmark/downsample"
#define XCORR_OUTPUT_DIR "/net/xia/data2-nobackup/zshen/Benchmark/xcorr"
#define SEGY_LIST "/home/zshen/Codes/cuda_debug/das_cc/segy.list"
#define SEGY_INFO "/home/zshen/Codes/cuda_debug/das_cc/segy.info"
#define STACK_XCORR_OUTPUT_DIR "/kuafu/zshen/cuda_debug/xcorr_Olanchasouth_ram_stack"

//decimate factor; change here to make sure same dt in both downsampled data
#define ARRAY1_DECIMATE 5 // range 2-7; it is possible to add more decimation in the future
#define ARRAY2_DECIMATE 5 

//define frequency range to be filtered
#define freql   0.1  	// lowest frequency
#define freqh   10.0  	// highest frequency

#define TEMPORAL_NORMALIZATION  2	// 0 for no normalization, 1 for one-bit, 2 for running-absolute-mean
#define SPECTRAL_WHITENING  1 		// 0 for no whitening, 1 for spectral whitening
#define MAX_LAG       30   		// [s] how long to be output e.g. +/- 10s
#define XCOR_SEGMENT  40  		// [s] xcor_window should be larger than MAX_LAG, the code will find the 2^n npts larger than XCOR_SEGMENT. But the XCOR_SEGMENT should close to 2^n

#define ARRAY1_CHANNELPERLOOP 1216	// number of channels used per xcorr loop for array1 (debug)
#define ARRAY2_CHANNELPERLOOP 1216	// number of channels used per xcorr loop for array2 (debug)

// Number of threads for preprocessing
#define NUM_THREADS 24

struct filelist{
        char file[256];
        double ut;};
