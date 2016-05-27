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

//#include <R.h>

#include "constants.h"

#include "getline.h"


#define DEBUG_R

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
    
    printf("Inside R_iso_seq, sequence is %s, length is %d\n",*seq,strlen(*seq));
	    
	    
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

      printf("calculating theoretical isotope distribution:\n");
      /* Get the theoretical mass (adding one for the charge) and the theoretical distribution */
      pep[k].pepmass = 1.0 + iso(pep, numpep, dist, 0);
      printf("theoretical isotope distribution (including charge):\n");
      printf("%0.2f: %0.9f\n", pep[k].pepmass, dist[k].prob[0]);
      printf("%0.2f: %0.9f \n" ,pep[k].pepmass + 1.0, dist[k].prob[1]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 2.0, dist[k].prob[2]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 3.0, dist[k].prob[3]);
      printf("%0.2f: %0.9f \n", pep[k].pepmass + 4.0, dist[k].prob[4]);

      printf("Populating R structures\n");
      for(m=0;m<numiso;m++){
      	  printf("Isotope %d: mass is %f, prob is %f\n",m,pep[k].pepmass,dist[k].prob[m]);
      	
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

    free (pep);
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
      printf("observed isotope distribution:\n");
      printf("%0.1f: %0.9f\n", pep[k].pepmass, dist[k].prob[0]);
      printf("%0.1f: %0.9f \n" ,pep[k].pepmass + 1.0, dist[k].prob[1]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 2.0, dist[k].prob[2]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 3.0, dist[k].prob[3]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 4.0, dist[k].prob[4]);

      printf("Populating R structures\n");
      for(m=0;m<numiso;m++){
      	  printf("Isotope %d: mass is %f, prob is %f\n",m,pep[k].pepmass,dist[k].prob[m]);
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

    free (pep);

  //}
  //return 0;
}






/****************************************
Procedure:      iso_main - this was the original
 *****************************************/
int iso_main (int argc, char *argv[])
{  

  /* Check that we have the correct number of files */
  if (argc != 2 ) printf("Correct program useage: ./iso [peptide filelist] \n");
  else
  {
    FILE *fin = NULL;
    int k, zs;
    int numpep;
    float ftmp;
    int tmp, dtmp; 
 
    /* read number of peptides that isotope distribution is to be calculated for */
    fin = fopen(argv[1], "r");
    tmp = fscanf(fin, "%d\n", &dtmp);
    numpep = dtmp; 
    PEPTIDE *pep = NULL;
    pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
    for (k = 0; k < numpep; k++)
    {
      pep[k].numseq = (int*) malloc (MAXLINE * sizeof(int));
      /* found is not used here but is necessary to share isodists.c with QtoE */
      pep[k].found = 1;
    }

    ISODIST  *dist = NULL;
    dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));  
    for (k = 0; k < numpep; k++)
    {
      dist[k].mass = (double*) malloc (5 * sizeof(double));
      dist[k].prob = (double*) malloc (5 * sizeof(double));
    }

    for (k = 0; k < numpep; k++)
    {
      /* read peptide's m/z (this is monoisotopic mass) */
      tmp = fscanf(fin, "%f ", &ftmp);
      /* peptide should be observed at monoisotopic mass + 1 (for charge) */
      pep[k].pepmass = ftmp + 1.0;  
      pep[k].length = isogetLine(fin, pep[k].numseq, &zs);

      pep[k].zs = zs;
    
      printf("calculating theoretical isotope distribution:\n");
      iso(pep, numpep, dist, 1);
      printf("observed isotope distribution:\n");
      printf("%0.1f: %0.9f\n", pep[k].pepmass, dist[k].prob[0]);
      printf("%0.1f: %0.9f \n" ,pep[k].pepmass + 1.0, dist[k].prob[1]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 2.0, dist[k].prob[2]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 3.0, dist[k].prob[3]);
      printf("%0.1f: %0.9f \n", pep[k].pepmass + 4.0, dist[k].prob[4]);

      for (k = 0; k < numpep; k++)
      {
        free (dist[k].mass);
        free (dist[k].prob);
      }
      free (dist);
    }
        
    free (pep);
    
  }
  return 0;
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







































