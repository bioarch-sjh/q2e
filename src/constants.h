/************************************************************************/
/* Copyright (C) 2009-2012 Julie Wilson                                 */
/* When you use this, send an email to: julie.wilson@york.ac.uk         */
/* with an appropriate reference to your work.                          */

/* This file is part of Q2E						*/

/* Q2E is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* Q2E is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

/************************************************************************/




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
	int found;
	double noise;
	double signal;
    double diff;
    double *peaks;
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
