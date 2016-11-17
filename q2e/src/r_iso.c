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
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <math.h>

//#include <R.h>

#include "constants.h"

#include "getline.h"

#include "isodists.h"
#include "r_iso.h"

//Use this to get printouts when needed (may have overdone it with the printfs!)
//#define DEBUG_R

/* We need these definitions in R_iso */
#define NUMELEMENTS 5
#define NUMLETTERS 26
#define NUMISOTOPES 5





/*C Procedures*/
float iso(PEPTIDE *pep, int numpep, ISODIST  *dist, int check);
void R_iso_seq(char **seq,  double *resultmass, double *resultprob, int *failed);
void R_iso_main(char **argv1,  double *resultmass, double *resultprob, int *failed);

/* forward declarations */
PEPTIDE * allocPepArray(const int numpep);
void freePepArray(PEPTIDE *pep, const int numpep);

ISODIST * allocDistArray(const int numpep);

void pepFromSequence(PEPTIDE *pep,char *sequence);



/*********************************************************************************
Procedure:     R iso seq, using data instead of a file reference

This will only work on one peptide at a time, and should be called repeatedly
if you want to get the peak distributions of >1 sequence.

 *********************************************************************************/
void R_iso_seq(char **seq,  double *resultmass, double *resultprob, int *failed){

	const int numpep = 1;
    const int numiso = 5;

	int k,m,tmp;
	float ftmp;
    FILE *fin = NULL;

#ifdef DEBUG_R
    printf("Inside R_iso_seq, sequence is %s, length is %d\n",*seq,strlen(*seq));
#endif

    /* allocate and initialise the PEPTIDE array */
    PEPTIDE *pep = NULL;
    pep = allocPepArray(numpep);


	/* allocate and initialise the ISODIST array */
    ISODIST  *dist = NULL;
    dist = allocDistArray(numpep);


	/* Now calculate the masses and probs */
    for (k = 0; k < numpep; k++)
    {
      ///* read peptide's m/z (this is monoisotopic mass) */
      //tmp = fscanf(fin, "%f ", &ftmp);
      //
      ///* peptide should be observed at monoisotopic mass + 1 (for charge) */
      //pep[k].pepmass = ftmp + 1.0;
      //
      ///* Get the sequence and its length from the remainder of the line*/
      //pep[k].length = isogetLine(fin, pep[k].numseq, &zs);
	  //
      //pep[k].zs = zs;
      pepFromSequence( &(pep[k]), *seq);

#ifdef DEBUG_R
      printf("calculating theoretical isotope distribution:\n");
#endif 
      /* Get the theoretical mass (adding one for the charge) and the theoretical distribution */
      pep[k].pepmass = 1.0 + R_iso(pep, numpep, dist, 0);



#ifdef DEBUG_R
      printf("theoretical isotope distribution (including charge):\n");
      printf("%0.2f: %0.9f\n", pep[k].pepmass, dist[k].prob[0]);
      printf("%0.2f: %0.9f \n" ,pep[k].pepmass + 1.0, dist[k].prob[1]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 2.0, dist[k].prob[2]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 3.0, dist[k].prob[3]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 4.0, dist[k].prob[4]);

      printf("Populating R structures\n");
#endif      
      
      
      for(m=0;m<numiso;m++){
#ifdef DEBUG_R      
      	  printf("Isotope %d: mass is %f, prob is %f\n",m,pep[k].pepmass,dist[k].prob[m]);
#endif

    	  resultmass[m] = pep[k].pepmass + (1.0 * m);
    	  resultprob[m] = dist[k].prob[m];
      }

      for (k = 0; k < numpep; k++)
      {
        free (dist[k].mass);
        free (dist[k].prob);
      }
      free (dist);
    }


    freePepArray(pep,numpep);
}







/****************************************
Procedure:     R iso main
 *****************************************/
void R_iso_main (char **argv1, double *resultmass, double *resultprob, int *failed)
{

  /* Check that we have the correct number of files */
  //if (argc != 2 ) printf("Correct program useage: ./iso [peptide filelist] \n");
  //else
  //{
    FILE *fin = NULL;
    int k,m, zs;
    int numpep;
    float ftmp;
    int tmp, dtmp;

    const int numiso=5;

    *failed = 0;

#ifdef DEBUG_R
    printf("Loading %s\n",argv1[0]);fflush(stdout);
#endif


    /**************
    read the number of peptides that isotope distribution is to be calculated for
     **************/

    /* Open the file */
    if((fin = fopen(argv1[0], "r"))==NULL){
        printf("Failed to open file %s, exiting isodists\n",argv1[0]);
	    *failed=1;
        return;
    }

	/* Read the number of peptides, store it in numpep */
    tmp = fscanf(fin, "%d\n", &dtmp);
    numpep = dtmp;


    /* allocate and initialise the PEPTIDE array */
    PEPTIDE *pep = NULL;
    //pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
    //for (k = 0; k < numpep; k++)
    //{
    //  pep[k].numseq = (int*) malloc (MAXLINE * sizeof(int));
    //  /* found is not used here but is necessary to share isodists.c with QtoE */
    //  pep[k].found = 1;
    //}
    pep = allocPepArray(numpep);

	/* allocate and initialise the ISODIST array */
    ISODIST  *dist = NULL;
    //dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));
    //for (k = 0; k < numpep; k++)
    //{
    //  dist[k].mass = (double*) malloc (5 * sizeof(double));
    //  dist[k].prob = (double*) malloc (5 * sizeof(double));
    //}
    dist = allocDistArray(numpep);


	/* Now calculate the masses and probs */
    for (k = 0; k < numpep; k++)
    {
      /* read peptide's m/z (this is monoisotopic mass) */
      tmp = fscanf(fin, "%f ", &ftmp);

      /* peptide should be observed at monoisotopic mass + 1 (for charge) */
      pep[k].pepmass = ftmp + 1.0;

      /* Get the sequence and its length from the remainder of the line*/
      pep[k].length = isogetLine(fin, pep[k].numseq, &zs);

      pep[k].zs = zs;

      iso(pep, numpep, dist, 1);
      
#ifdef DEBUG_R
      printf("observed isotope distribution:\n");
      printf("%0.1f: %0.9f\n", pep[k].pepmass, dist[k].prob[0]);
      printf("%0.1f: %0.9f \n" ,pep[k].pepmass + 1.0, dist[k].prob[1]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 2.0, dist[k].prob[2]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 3.0, dist[k].prob[3]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 4.0, dist[k].prob[4]);

      printf("Populating R structures\n");
#endif

      for(m=0;m<numiso;m++){
#ifdef DEBUG_R
      	  printf("Isotope %d: mass is %f, prob is %f\n",m,pep[k].pepmass,dist[k].prob[m]);
#endif 
    	  resultmass[m] = pep[k].pepmass + (1.0 * m);
    	  resultprob[m] = dist[k].prob[m];
      }

      for (k = 0; k < numpep; k++)
      {
        free (dist[k].mass);
        free (dist[k].prob);
      }
      free (dist);
    }

    freePepArray(pep,numpep);

  //}
  //return 0;
}






/******************************************************************
 Utility declarations
 ******************************************************************/
PEPTIDE * allocPepArray(const int numpep){

	int k;

    /* allocate and initialise the PEPTIDE array */
    PEPTIDE *pep = NULL;
    pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
    for (k = 0; k < numpep; k++)
    {
      pep[k].numseq = (int*) malloc (MAXLINE * sizeof(int));
      /* found is not used here but is necessary to share isodists.c with QtoE */
      pep[k].found = 1;
    }

    return pep;
}


void freePepArray(PEPTIDE *pep, const int numpep){

	int k;

    for (k = 0; k < numpep; k++)
    {
      free(pep[k].numseq);
    }

    free(pep);
    pep = NULL;
}


/* allocate and initialise the ISODIST array */
ISODIST * allocDistArray(const int numpep){

	int k;

	ISODIST *dist = NULL;
    dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));
    for (k = 0; k < numpep; k++){
      dist[k].mass = (double*) malloc (5 * sizeof(double));
      dist[k].prob = (double*) malloc (5 * sizeof(double));
    }

    return dist;
}


/**/
void pepFromSequence(PEPTIDE *pep,char *sequence){

	int i=0;
	pep->zs=0;

	while(sequence[i]){
		pep->numseq[i] = sequence[i]-65;
		if(pep->numseq[i]==25){
			pep->zs++;
		}
		i++;
	}
	pep->length=strlen(sequence);
}



/**************************************************
int isogetLine(FILE *fp, int v[], int *zs)
{
    int i, ch;
    int z = 0;

    i = 0;
    ch = fgetc(fp);

    v[i] = ch-65;

	while (ch != 10)
	{
	    i++;
	    ch = fgetc(fp);
	    v[i] = ch-65;
	    if (v[i] == 25) z++;
	}

    (*zs) = z;

 	return i;
}
***************************************************/
/****************************************************
 *
 * FROM HERE ON, WE ARE COPYING FROM R/src/isodists.C
 *
 */


/****************************************
Procedure:      R_iso
*****************************************/
float R_iso (PEPTIDE *pep, int numpep, ISODIST *dist, int check)
{
  //FILE *fp = NULL;
  /* number of elements in isotope table */

#ifdef DEBUG_R
  printf("Loading amino table\n");fflush(stdout);
#endif


  /* THIS STRUCTURE SHOULD BE ITENTICAL TO THE FILE "aminomasses" */
  //const int AMINODATA[numletters][numelements] = {

  const int TAMINODATA[3][2] = {{0,1},{2,3},{4,5}};

  int x = TAMINODATA[0][0];

#ifdef DEBUG_R  
  printf("x = %d\n",x);
#endif

  const int AMINODATA[NUMLETTERS][NUMELEMENTS] = {
  //  ELEMENT:
  //  C   H   N   O   S   //Amino Acid
    { 3,  7,  1,  2,  0}, //A
	{ 0, -1,  0, -1,  0}, //B
	{ 3,  7,  1,  2,  1}, //C
	{ 4,  7,  1,  4,  0}, //D
	{ 5,  9,  1,  4,  0}, //E
	{ 9, 11,  1,  2,  0}, //F
	{ 2,  5,  1,  2,  0}, //G
	{ 6,  9,  3,  2,  0}, //H
	{ 6, 13,  1,  2,  0}, //I
	{ 0,  0,  0,  0,  0}, //J
	{ 6, 14,  2,  2,  0}, //K
	{ 6, 13,  1,  2,  0}, //L
	{ 5, 11,  1,  2,  1}, //M
	{ 4,  8,  2,  3,  0}, //N
	{ 0,  0,  0,  0,  0}, //O
	{ 5,  9,  1,  2,  0}, //P
	{ 5, 10,  2,  3,  0}, //Q
	{ 6, 14,  4,  2,  0}, //R
	{ 3,  7,  1,  3,  0}, //S
	{ 4,  9,  1,  3,  0}, //T
	{ 0,  0,  0,  0,  0}, //U
	{ 5, 11,  1,  2,  0}, //V
	{11, 12,  2,  2,  0}, //W
	{ 2,  2,  0,  2,  0}, //X
	{ 9, 11,  1,  3,  0}, //Y
	{ 0,  0,  0,  1,  0}  //Z
  };


  /* number of different isotopes in isotope table */
  int numiso = 5;//TODO: This should use the constant NUMISOTOPES

  double imass;

  int i, j, k;

  double *peptide = NULL;
  peptide = (double*) malloc (NUMELEMENTS *  sizeof(double));

  ISOTAB  *element = NULL;
  element = (ISOTAB*) malloc (NUMELEMENTS *  sizeof(ISOTAB));
  for (i = 0; i < NUMELEMENTS; i++)
  {
    element[i].isotope = (double*) malloc (numiso * sizeof(double));
    element[i].abundance = (double*) malloc (numiso * sizeof(double));
  }

	AMINO  *letter = NULL;
	letter = (AMINO*) malloc (NUMLETTERS *  sizeof(AMINO));
	for (i = 0; i < NUMLETTERS; i++){
		letter[i].element = (int*) malloc (NUMELEMENTS * sizeof(int));
	}

  //TODO: Although we've replaced the file aminomasses with the AMINODATA array, we need a better way to load other amino masses if we need to.
  /* read in elements for each amino acid */
  //fp = fopen(AMINO_FILE, "r");
  for (i = 0; i < NUMLETTERS; i++)
  {
    for (k = 0; k < NUMELEMENTS; k++)
    {
      //tmp = fscanf(fp,"%d ", &dtmp);
      //letter[i].element[k] = dtmp;
      letter[i].element[k] = AMINODATA[i][k];
    }
  }
  //fclose(fp);

  R_loadIsotopeTable(/*ISO_TABLE,*/ element, NUMELEMENTS, numiso);
#ifdef DEBUG_R
  printf("success reading iso table\n");
#endif

#ifdef DEBUG_R
	//Print the table for debugging
#endif


	for (i = 0; i < numpep; i++)
	{
		/*find the peptides that isotope distributions are required for */
		if (pep[i].found == 1)
		{

			for (k = 0; k < NUMELEMENTS; k++)
			{

#ifdef DEBUG_R
				printf("loading peptide element %d\n",k);fflush(stdout);
#endif
				peptide[k] = 0.0;
				for (j = 0; j < pep[i].length; j++)
				{

#ifdef DEBUG_R
					printf("loading peptide letter %d (=%c)\n",j,  pep[i].numseq[j]+65);fflush(stdout);
#endif

					peptide[k] = peptide[k] + (float)letter[ pep[i].numseq[j] ].element[k];
				}
			}

			/* remove 2H10 for each link between amino acids */
			peptide[1] = peptide[1] - 2.0*(float)(pep[i].length - pep[i].zs - 1);
			peptide[3] = peptide[3] - (float)(pep[i].length - pep[i].zs - 1);


#ifdef DEBUG_R
			printf("getting iso dist..\n");fflush(stdout);
#endif

			imass = getIsoDist(i, element, NUMELEMENTS, peptide, dist);

			printf("calculated mass %f\n", imass);
			if(check)
				if (fabs(imass-pep[i].pepmass) > 1.5) printf("ERROR: SOMETHING IS WRONG HERE, THE DIFERENCE BETWEEN THE GIVEN AND CALCULATED MASSES IS %d\n", (int)(fabs(imass-pep[i].pepmass)+0.5));
			}
		}

	for (i = 0; i < NUMELEMENTS; i++)
	{
		free(element[i].isotope);
		free(element[i].abundance);
	}
	free(element);

	for (i = 0; i < NUMLETTERS; i++){
		free(letter[i].element);
	}
	free(letter);

	free(peptide);

	return imass;
}


/*************************************************************************************************************
Procedure: readIsotopeTable
inputs the isotope info
************************************************************************************************************/
void R_loadIsotopeTable(/*char *filename,*/ ISOTAB *element, const int num, const int numiso)
{
//  FILE  *ft = NULL;
//ft = fopen(filename, "r");

#ifdef DEBUG_R
  printf("Loading iso table\n");fflush(stdout);
#endif


  if(num != 5)
    printf("WARNING! Number of elements is not 5 - defaults won't work!\n");
  if(numiso != 5)
    printf("WARNING! Number of isotopes is not 5 - defaults won't work!\n");

  //char name[1024];
  //double ftmp, ftmp1;
  int i, j;//, tmp;

//TODO: Make these global variables if they are consts
//TODO: Need to read these values into the table:
/*
C 12.01070 12.00000 98.89220 13.00335 1.10780  0.00000 0.00000  0.00000 0.00000  0.00000 0.00000
H  1.00794  1.00783 99.98443  2.01410 0.01557  0.00000 0.00000  0.00000 0.00000  0.00000 0.00000
N 14.00670 14.00307 99.63370 15.00011 0.36630  0.00000 0.00000  0.00000 0.00000  0.00000 0.00000
O 15.99940 15.99491 99.76280 16.99913 0.03720 17.99916 0.20004  0.00000 0.00000  0.00000 0.00000
S 32.06500 31.97207 95.01800 32.97146 0.75000 33.96787 4.21500  0.00000 0.00000 35.96708 0.01700
*/

#ifdef DEBUG_R
  printf("loading name array\n");fflush(stdout);
#endif

  const char ISOname[NUMELEMENTS][3] = {{"C\0"},{"H\0"},{"N\0"},{"O\0"},{"S\0"}};


#ifdef DEBUG_R
  printf("loading float arrays\n");fflush(stdout);
#endif


  /* avmass comes from the first column in the above */
  const float ISOavmass[NUMELEMENTS] = {12.01070,
                                 1.00794,
                                14.00670,
                                15.99940,
                                32.06500};


  const float ISOmass[NUMELEMENTS][NUMISOTOPES] = {
                              {12.00000, 13.00335,  0.00000, 0.00000,  0.00000},
                              { 1.00783,  2.01410,  0.00000, 0.00000,  0.00000},
                              {14.00307, 15.00011,  0.00000, 0.00000,  0.00000},
                              {15.99491, 16.99913, 17.99916, 0.00000,  0.00000},
                              {31.97207, 32.97146, 33.96787, 0.00000, 35.96708}};



  const float ISOabundance[NUMELEMENTS][NUMISOTOPES] = {
                              {98.89220,  1.10780,  0.00000, 0.00000,  0.00000},
                              {99.98443,  0.01557,  0.00000, 0.00000,  0.00000},
                              {99.63370,  0.36630,  0.00000, 0.00000,  0.00000},
                              {99.76280,  0.03720,  0.20004, 0.00000,  0.00000},
                              {95.01800,  0.75000,  4.21500, 0.00000,  0.01700}};


#ifdef DEBUG_R
  printf("Basic arrays set\n");fflush(stdout);
#endif

  for (i = 0; i < NUMELEMENTS; i++)
  {
    //tmp = fscanf(ft, "%s %lf",  name, &ftmp);
    //strcpy(element[i].name, name);
    //element[i].avmass = ftmp;

#ifdef DEBUG_R
    printf("loading element %d\n",i);fflush(stdout);
#endif


    strcpy(element[i].name, ISOname[i]);
    element[i].avmass = ISOavmass[i];


    for (j = 0; j < NUMISOTOPES; j++)
    {
      //tmp = fscanf(ft, "%lf %lf",  &ftmp, &ftmp1);
      //element[i].isotope[j] = ftmp;
      //element[i].abundance[j] = ftmp1/100.0;

      element[i].isotope[j] = ISOmass[i][j];
      element[i].abundance[j] = ISOabundance[i][j]/100.0;
    }
  }
}


