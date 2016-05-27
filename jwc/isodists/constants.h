#ifndef CONSTANTS_H
#define CONSTANTS_H

#define MAX_STRING 1000
#define MAXLINE 1000000

typedef struct pepinfo{
  double pepmass;
  int length;
  int zs;	
  int *numseq;
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
