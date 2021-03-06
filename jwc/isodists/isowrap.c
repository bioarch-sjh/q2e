#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "constants.h"

/*C Procedures*/
/*TODO: These need to go in their own header files!*/
float iso(PEPTIDE *pep, int numpep, ISODIST  *dist); 
int getLine(FILE *fp, int v[], int *zs);

/****************************************
Procedure:      main
 *****************************************/
int main (int argc, char *argv[])
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
      pep[k].length = getLine(fin, pep[k].numseq, &zs); 

      pep[k].zs = zs;
    
      iso(pep, numpep, dist);
 
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

/*****************************************************
Procedure: getLine
Description: reads in whole line as a string
******************************************************/
int getLine(FILE *fp, int v[], int *zs)
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


