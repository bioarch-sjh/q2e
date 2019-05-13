#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "constants.h"

#define DEBUG_R

/*C Procedures*/
void iso(PEPTIDE *pep, int numpep, ISODIST  *dist); 
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


#ifdef DEBUG_R
  printf("opening file\n");fflush(stdout);
#endif

 
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
      pep[k].length = getLine(fin, pep[k].numseq, &zs); 

      pep[k].zs = zs;
    
      iso(pep, numpep, dist);
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


