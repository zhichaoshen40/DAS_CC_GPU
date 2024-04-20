#define SEGY_HEAD_SIZE	240	/* segy header size in bytes */
struct segy		/* standard SEG-Y trace */
	{
	int trseql;	/*   0: trace sequence no. in line */
	int trseqr;	/*   4: trace sequence no. in reel */
	int frecno;	/*   8: original field record no. */
	int chan;	/*  12: trace no. within original field rec */
	int sp;		/*  16: energy source point number */
	int cdpensno;	/*  20: CDP ensemble no. */
	int enstrno;	/*  24: trace number within CDP ensemble */
	short trid;	/*  28: trace id code */
	short vsumno;	/*  30: no. of vertically summed traces in trace */
	short hsumno;	/*  32: no. of horizontally summed traces in trace */
	short use;	/*  34: data use */
	int offset;	/*  36: distance from source point to receiver */
	int gz;		/*  40: elevation of receiver group */
	int szsurf;	/*  44: surface elevation at source */
	int sdepth;	/*  48: source depth below surface */
	int gdatz;	/*  52: datum elevation at receiver group */
	int sdatz;	/*  56: datum elevation at source */
	int swatdep;	/*  60: water depth at source */
	int gwatdep;	/*  64: water depth at receiver group */
	short lscaler;	/*  68: elevation and depth scaler */
	short coscaler;	/*  70: coordinate scaler */
	int sx;		/*  72: x source coordinate */
	int sy;		/*  76: y source coordinate */
	int gx;		/*  80: x receiver group coordinate */
	int gy;		/*  84: y receiver group coordinate */
	short counit;	/*  88: coordinate units */
	short vweath;	/*  90: weathering velocity */
	short vsbweath;	/*  92: subweathering velocity */
	short tus;	/*  94: uphole time at source */
	short tug;	/*  96: uphole time at receiver group */
	short dts;	/*  98: source static correction */
	short dtg;	/* 100: receiver group static correction */
	short dttotal;	/* 102: total static applied */
	short tlA;	/* 104: lag time A */
	short tlB;	/* 106: lag time B */
	short delay;	/* 108: recording delay time */
	short tmute0;	/* 110: mute time start */
	short tmutef;	/* 112: mute time end */
	unsigned short nt;	/* 114: no. of samples in trace */
	short dt;	/* 116: sample interval in trace */
	short gaintype;	/* 118: gain type of field instruments */
	short kgain;	/* 120: instrument gain constant */
	short gain0;	/* 122: instrument early or initial gain */
	short correl;	/* 124: correlation switch */
	short swf0;	/* 126: sweep frequency at start */
	short swff;	/* 128: sweep frequency at end */
	short swt;	/* 130: sweep length */
	short swtype;	/* 132: sweep type */
	short swtrtpt0;	/* 134: sweep trace taper length at start */
	short swtrtptf;	/* 136: sweep trace taper length at end */
	short swtaptype;/* 138: sweep taper type */
	short alff;	/* 140: alias filter frequency */
	short alfsl;	/* 142: alias filter slope */
	short nchff;	/* 144: notch filter frequency */
	short nchfsl;	/* 146: notch filter slope */
	short lcf;	/* 148: low cut frequency */
	short hcf;	/* 150: high cut frequency */
	short lcsl;	/* 152: low cut slope */
	short hcsl;	/* 154: high cut slope */
	short year;	/* 156: year data recorded */
	short jday;	/* 158: day of year */
	short hour;	/* 160: hour of day */
	short minute;	/* 162: minute of hour */
	short second;	/* 164: second of minute */
	short tbasis;	/* 166: time basis code */
	short trwt;	/* 168: trace weighting factor */
	short rsw1gno;	/* 170: geophone group no. of roll switch pos. 1 */
	short frtr1gno;	/* 172: grp. no. of trace 1 w/i orig. field rec. */
	short frtrfgno;	/* 174: grp. no. of last trace w/i field rec. */
	short gapsize;	/* 176: gap size */
	short overtrav;	/* 178: overtravel assoc. w/ taper at line ends */
	/* PASSCAL mods */
	char  name[6];	/* 180: station name; */
	char  serial[8];/* 186: sensor serial # */
	char  cha[4];   /* 194: channel name */
	short tstatic;  /* 198: static shift im millisec */
	int dtbig;      /* 200: dt in microsecs */
	short dformat;	/* 204: data format 0=16-bit int, 1=32-bit int, */
	short msec;	/* 206: millisec for 1st sample */
	short year2;	/* 208: year */
	short jday2;	/* 210: day of year */
	short hour2;	/* 212: hour of day */
	short minute2;	/* 214: minute of hour */
	short second2;	/* 216: second of minute */
	short msec2;	/* 218: milli sec. */
	float scale;	/* 220: scale factor */
	short serial2;	/* 224: serial number */
	short unused2;	/* 226: unused */
	int  ntbig;	/* 228: nt big */
	int   max;	/* 232: max value of trace */
	int   min;	/* 236: max value of trace */
	};

