#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "constants.h"
#include "r_iso.h"

/*C Procedures*/
/*TODO: These need to go in their own header files!*/
float iso(PEPTIDE *pep, int numpep, ISODIST  *dist); 
int R_iso_getLine(FILE *fp, int v[], int *zs);

/******************************************
 Tests the R_iso_seq function that is called
 in the R q2e package
 ******************************************/
int test_R_iso_seq(char * seq){


	const int nres = 5;
	//char seq[] = "IGQPGAVGPAGIR";
	double resmass[nres];
	double resprob[nres];
	int failed = 0;

	//void R_iso_seq(char **seq,  double *resultmass, double *resultprob, int *failed);
	R_iso_seq(&seq, resmass, resprob, &failed);

	return failed;
}




/****************************************
Procedure:      main

This is based on Julie Wilson's original main() for iso,
but it tests
 *****************************************/
int main (int argc, char *argv[])
{  
	/******************************
	 * test the sequence call from R - see file peptideList in tests/iso/
	 * you can get the masses from here: http://education.expasy.org/student_projects/isotopident/htdocs/
	 */

	int test;

	test = test_R_iso_seq("IGQPGAVGPAGIR");
	test = test_R_iso_seq("GPPGPQGAR");
	test = test_R_iso_seq("ACDEFGHIKLMNPQRSTVWY");
	test = test_R_iso_seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ");//<- this sequence contains illegal characters

	printf("================================================================\n");


	/* Check that we have the correct number of files */
	if (argc != 2 ) printf("Correct program usage: ./iso [peptide filelist] \n");
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
		  printf("monoisotopic mass: %0.1f\n", ftmp);
		  /* peptide should be observed at monoisotopic mass + 1 (for charge) */
		  pep[k].pepmass = ftmp + 1.0;
		  pep[k].length = R_iso_getLine(fin, pep[k].numseq, &zs);

		  pep[k].zs = zs;
		}
		fclose(fin);

		  iso(pep, numpep, dist);
		  for (k = 0; k < numpep; k++)
		  {
			free (dist[k].mass);
			free (dist[k].prob);
		  }
		  free (dist);

		  for (k = 0; k < numpep; k++)
			{
			free(pep[k].numseq);
			}
		free (pep);
    
	}
	return 0;
}



/*****************************************************
Procedure: getLine- OBSOLETE - see getline.h and getline.c
Description: reads in whole line as a string
*****************************************************/
int R_iso_getLine(FILE *fp, int v[], int *zs)
{
	int i, ch;
    int z = 0;

    i = 0;
    ch = fgetc(fp);
    printf("%d\n", ch);

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
/*******************************************************/


