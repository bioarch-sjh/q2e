#ifndef CONSTANTS_H
#define CONSTANTS_H


#define STEP 0.1

#define MAX_STRING 1000
#define MAXLINE 1000000

typedef struct alldata{
  char name[1024];
  double *vars;
  double *peaks;
  int *num;
}SAMPLE;

typedef struct pep{
    char name[1024];
	int found;
	double noise;
	double signal;
    double diff;
    double *peaks;
    double *betas;
}SAMPLEPEP;

typedef struct pepinfo{
  double pepmass;
  int found;
  int length;
  int zs;	
  int pg;
  int cbm;
  int *numseq;
  int nq;
  SAMPLEPEP *sample;
}PEPTIDE;

typedef struct table {
	double avmass;
	double *isotope;
	double *abundance;
	char name[1024];
} ISOTAB;

typedef struct distribution {
	double *mass;
	double *prob;
} ISODIST;

typedef struct mass {
	int *element;
	/* 0 = C; 1 = H; 2 = N; 3 = O; 4 = S */
} AMINO;

#endif
