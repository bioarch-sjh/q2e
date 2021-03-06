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
#include "constants.h"

#include "getline.h"



/*C Procedures*/
int pepfiles (char *name, SAMPLE *data, int nsamples, int nreps, int nmasses, PEPTIDE *pep, int numpep, int FITPEAKS, double FIRSTMASS, double SNRLIM);              
void iso(PEPTIDE *pep, int numpep, ISODIST  *dist); 
//int getLine(FILE *fp, int v[], int *zs, int *pg, int *cbm);
void runga (char *name, PEPTIDE *pep, int numpep, ISODIST  *dist, int msamples, int FITPEAKS, double GALIM);




/****************************************
Procedure:	Rmain - main function for calling from R
 *****************************************/



void Rmain(char **argv1, char **argv2, char **argv3){


	  /* Check that we have the correct number of files */
	  if(0)// (argc != 4 )
		  printf("ERROR: Correct program useage: ./QtoEcalc [data filelist] [peptide filelist] [parameter file] \n");
	  else
	  {
	    FILE *fp = NULL;
	    FILE *fin = NULL;
	    int i, j, k, ii, nsamples, nreps, zs, cbm, pg;
	    int numpep, newnumpep, missing;
	    float ftmp, ftmp1;
	    int tmp, dtmp, dtmp1;
	    char name[MAX_STRING];

	    int FITPEAKS;
	    double GALIM;
	    double SNRLIM;
	    double FIRSTMASS;
	    double LASTMASS;
	    int CSV;

	    printf("reading parameters from %s\n", argv3[0]);
	    fp = fopen(argv3[0], "r");
	    if (fp == NULL) printf("ERROR: couldn't open parameter file %s\n", argv3[0]);
	    else
	    {
	      tmp = fscanf(fp, "%s %d\n", name, &dtmp);
	      FITPEAKS = dtmp;
	      printf("%s = %d\n",name, FITPEAKS);
	      tmp = fscanf(fp, "%s %f\n", name, &ftmp);
	      GALIM = ftmp;
	      printf("%s = %0.2f\n",name, GALIM);
	      tmp = fscanf(fp, "%s %f\n", name, &ftmp);
	      SNRLIM = ftmp;
	      printf("%s = %0.1f\n",name, SNRLIM);
	      tmp = fscanf(fp, "%s %f\n", name, &ftmp);
	      FIRSTMASS = ftmp;
	      printf("%s = %0.1f\n",name, FIRSTMASS);
	      tmp = fscanf(fp, "%s %f\n", name, &ftmp);
	      LASTMASS = ftmp;
	      printf("%s = %0.1f\n",name, LASTMASS);
	      tmp = fscanf(fp, "%s %d\n", name, &dtmp);
	      CSV = dtmp;
	      printf("%s = %d\n",name, CSV);
	      fclose(fp);
	    }


	    printf("Parameters have loaded correctly\n");


	    /* keepname for output files */
	    char outfilename[1024];
	    strcpy(outfilename, argv1[0]);

	    printf("reading data files from %s\n", argv1[0]);
	    fp = fopen(outfilename, "r");
	    if (fp == NULL)
	    {
	      printf("ERROR: couldn't open file %s\n", argv1[0]);
	    }
	    else
	    {
	      missing = 0;
	      tmp = fscanf(fp, "%d %d\n", &dtmp, &dtmp1);
	      nsamples = dtmp;
	      nreps = dtmp1;
	      int msamples = nsamples/nreps;
	      printf("%d individual samples with %d replicates each of %d samples\n", nsamples, nreps, msamples);

	      /* number of masses to be kept: nmasses = lastmass - firstmass in step m/z units */
	      int nmasses = (int)((LASTMASS - FIRSTMASS)/STEP);

	      SAMPLE *data = NULL;
	      data = (SAMPLE *) malloc (nsamples * sizeof(SAMPLE));
	      for (i = 0; i < nsamples; i++)
	      {
	        data[i].vars = (double *) malloc (nmasses * sizeof(double));
	        data[i].num = (int *) malloc (nmasses * sizeof(int));
	        for (j = 0; j < nmasses; j++)
	        {
	          data[i].vars[j] = 0.0;
	          data[i].num[j] = 0;
	        }
	      }

	      double mult = 1.0/STEP;
	      /* read in filenames */
	      for (i = 0; i < nsamples; i++)
	      {
	        tmp = fscanf(fp, "%s\n", name);
	        fin = fopen(name, "r");
	        if (fin == NULL)
	        {
	          printf("ERROR: couldn't find file %s\n", name);
	          missing ++;
	        }
	        else
	        {
	          printf("%s\n", name);
	          strcpy(data[i].name, name);
	          /* read in data */
	          tmp = 0;
	          while (tmp != EOF)
	          {
	            /* read in mass */
	            if (CSV == 1) tmp = fscanf(fin, "%f,%f\n", &ftmp, &ftmp1);
	            else
	            {
	              tmp = fscanf(fin, "%f %d\n", &ftmp, &dtmp);
	              ftmp1 = (float)dtmp;
	            }
	            if (tmp != EOF)
	            {
	              ii = (int)((ftmp - FIRSTMASS)*mult + 0.5);
	              if ((ii > 0) && (ii < nmasses))
	              {
	                data[i].vars[ii] += ftmp1;
	                data[i].num[ii]++;
		          }
		        }
		      }
	   	      fclose(fin);
		    }
	      }
	      fclose(fp);

	      if (missing > 0) printf("You have %d samples missing\n", missing);
	      if (missing == nsamples) printf("ERROR: couldn't find any of your files!\n");
	      else
	      {
	        for (i = 0; i < nsamples; i++)
	        {
	          for (k = 0; k < nmasses; k++)
	          {
	            if (data[i].num[k] > 0) data[i].vars[k] /= (float) data[i].num[k];
	          }
	        }

	        /* read number of peptides */
	        fin = fopen(argv2[0], "r");
	        if (fin == NULL) printf("ERROR: couldn't open peptide file %s\n", argv2[0]);
	        else
	        {
	          tmp = fscanf(fin, "%d\n", &dtmp);
	          numpep = dtmp;
	          PEPTIDE *pep = NULL;
	          pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
	          for (i = 0; i < numpep; i++)
	          {
	            pep[i].sample = (SAMPLEPEP*) malloc (msamples * sizeof(SAMPLEPEP));
	            pep[i].numseq = (int*) malloc (MAXLINE * sizeof(int));
	            for (j = 0; j < msamples; j++)
	            {
	              pep[i].sample[j].peaks = (double*) malloc ((FITPEAKS + 2) * sizeof(double));
	              pep[i].sample[j].found = -1;
	            }
	          }

	          printf("checking %d m/zs:\n", numpep);

	          for (k = 0; k < numpep; k++)
	          {

		    printf("%d m/zs...\n", k);
	            /* read peptide's m/z (this is monoisotopic mass) and number of Qs */
	            tmp = fscanf(fin, "%f %d ", &ftmp, &dtmp);
	            pep[k].nq = dtmp;
	            /* peptide should be observed at monoisotopic mass + 1 (for charge) */
	            pep[k].pepmass = ftmp + 1.0;
	            pep[k].length = getLine(fin, pep[k].numseq, &zs, &pg, &cbm);
	            pep[k].zs = zs;
	            pep[k].pg = pg;
	            pep[k].cbm = cbm;
	            /*printf("%f %d\n", pep[k].pepmass, pep[k].nq);*/
	          }
	          newnumpep = pepfiles(outfilename, data, nsamples, nreps, nmasses, pep, numpep, FITPEAKS, FIRSTMASS, SNRLIM);
	          if (newnumpep > 0)
	          {
	            /* note: this will not be filled for those peptides with pep[i].found = 0 */
	            ISODIST  *dist = NULL;
	            dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));
	            for (k = 0; k < numpep; k++)
	            {
	              dist[k].mass = (double*) malloc (5 * sizeof(double));
	              dist[k].prob = (double*) malloc (5 * sizeof(double));
	            }

	            iso(pep, numpep, dist);

	            /* number of peaks in isotope distribution to be fitted */
	            runga (outfilename, pep, numpep, dist, msamples, FITPEAKS, GALIM);

	            for (k = 0; k < newnumpep; k++)
	            {
	              free (dist[k].mass);
	              free (dist[k].prob);
	            }
	            free (dist);

	          }
	          else printf("none of these peptide found!\n");

	          for (i = 0; i < nsamples; i++)
	          {
	            free (data[i].vars);
	          }
	          free (data);

	          for (i = 0; i < numpep; i++)
	          {
	            for (j = 0; j < msamples; j++)
	            {
	              free (pep[i].sample[j].peaks);
	            }
	            free (pep[i].sample);
	          }
	          free (pep);
	        }
	      }
	    }
	  }
	  //return 0;

	printf("FINISED Q2E\n");



}



/****************************************
Procedure:      main
 *****************************************/

int main(int argc, char *argv[]){
	if (argc != 4 ) printf("ERROR: Correct program useage: ./QtoEcalc [data filelist] [peptide filelist] [parameter file] \n");
	else{
		Rmain(&(argv[1]),&(argv[2]),&(argv[3]));
	}
	return 0;
}


int oldmain (int argc, char *argv[])
{  

  /* Check that we have the correct number of files */
  if (argc != 4 ) printf("ERROR: Correct program useage: ./QtoEcalc [data filelist] [peptide filelist] [parameter file] \n");
  else
  {
    FILE *fp = NULL;
    FILE *fin = NULL;
    int i, j, k, ii, nsamples, nreps, zs, cbm, pg;
    int numpep, newnumpep, missing;
    float ftmp, ftmp1;
    int tmp, dtmp, dtmp1; 
    char name[MAX_STRING];
    
    int FITPEAKS;
    double GALIM;
    double SNRLIM;
    double FIRSTMASS;
    double LASTMASS;
    int CSV;
    
    printf("reading parameters from %s\n", argv[3]); 
    fp = fopen(argv[3], "r");
    if (fp == NULL) printf("ERROR: couldn't open parameter file %s\n", argv[3]);
    else
    {
      tmp = fscanf(fp, "%s %d\n", name, &dtmp); 
      FITPEAKS = dtmp;
      printf("%s = %d\n",name, FITPEAKS);
      tmp = fscanf(fp, "%s %f\n", name, &ftmp); 
      GALIM = ftmp;    
      printf("%s = %0.2f\n",name, GALIM);
      tmp = fscanf(fp, "%s %f\n", name, &ftmp); 
      SNRLIM = ftmp;
      printf("%s = %0.1f\n",name, SNRLIM);
      tmp = fscanf(fp, "%s %f\n", name, &ftmp); 
      FIRSTMASS = ftmp;
      printf("%s = %0.1f\n",name, FIRSTMASS);
      tmp = fscanf(fp, "%s %f\n", name, &ftmp); 
      LASTMASS = ftmp;
      printf("%s = %0.1f\n",name, LASTMASS);
      tmp = fscanf(fp, "%s %d\n", name, &dtmp); 
      CSV = dtmp;
      printf("%s = %d\n",name, CSV);
      fclose(fp);         
    }
 

    printf("Parameters have loaded correctly\n");


    /* keepname for output files */
    char outfilename[1024];
    strcpy(outfilename, argv[1]);

    printf("reading data files from %s\n", argv[1]); 
    fp = fopen(outfilename, "r");
    if (fp == NULL) 
    {
      printf("ERROR: couldn't open file %s\n", argv[1]);
    }
    else
    {
      missing = 0;
      tmp = fscanf(fp, "%d %d\n", &dtmp, &dtmp1); 
      nsamples = dtmp;
      nreps = dtmp1;
      int msamples = nsamples/nreps;
      printf("%d individual samples with %d replicates each of %d samples\n", nsamples, nreps, msamples);

      /* number of masses to be kept: nmasses = lastmass - firstmass in step m/z units */
      int nmasses = (int)((LASTMASS - FIRSTMASS)/STEP);
    
      SAMPLE *data = NULL;
      data = (SAMPLE *) malloc (nsamples * sizeof(SAMPLE));
      for (i = 0; i < nsamples; i++)
      {
        data[i].vars = (double *) malloc (nmasses * sizeof(double));
        data[i].num = (int *) malloc (nmasses * sizeof(int));
        for (j = 0; j < nmasses; j++)
        {
          data[i].vars[j] = 0.0;
          data[i].num[j] = 0;
        }
      }  
   
      double mult = 1.0/STEP;
      /* read in filenames */
      for (i = 0; i < nsamples; i++)
      {
        tmp = fscanf(fp, "%s\n", name); 
        fin = fopen(name, "r");
        if (fin == NULL)
        {
          printf("ERROR: couldn't find file %s\n", name);
          missing ++;
        }
        else
        {
          printf("%s\n", name);
          strcpy(data[i].name, name);
          /* read in data */
          tmp = 0;
          while (tmp != EOF)
          {
            /* read in mass */
            if (CSV == 1) tmp = fscanf(fin, "%f,%f\n", &ftmp, &ftmp1);
            else
            {
              tmp = fscanf(fin, "%f %d\n", &ftmp, &dtmp);
              ftmp1 = (float)dtmp;
            }
            if (tmp != EOF)
            {
              ii = (int)((ftmp - FIRSTMASS)*mult + 0.5);
              if ((ii > 0) && (ii < nmasses))
              {
                data[i].vars[ii] += ftmp1;
                data[i].num[ii]++;
	          }
	        }
	      } 	
   	      fclose(fin);
	    }
      }
      fclose(fp);         
  
      if (missing > 0) printf("You have %d samples missing\n", missing);
      if (missing == nsamples) printf("ERROR: couldn't find any of your files!\n");
      else
      {
        for (i = 0; i < nsamples; i++)
        {
          for (k = 0; k < nmasses; k++) 
          {
            if (data[i].num[k] > 0) data[i].vars[k] /= (float) data[i].num[k];
          }
        }
    
        /* read number of peptides */
        fin = fopen(argv[2], "r");
        if (fin == NULL) printf("ERROR: couldn't open peptide file %s\n", argv[2]);
        else
        {
          tmp = fscanf(fin, "%d\n", &dtmp);
          numpep = dtmp; 
          PEPTIDE *pep = NULL;
          pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
          for (i = 0; i < numpep; i++)
          {
            pep[i].sample = (SAMPLEPEP*) malloc (msamples * sizeof(SAMPLEPEP));
            pep[i].numseq = (int*) malloc (MAXLINE * sizeof(int));
            for (j = 0; j < msamples; j++)
            {
              pep[i].sample[j].peaks = (double*) malloc ((FITPEAKS + 2) * sizeof(double));
              pep[i].sample[j].found = -1;
            }
          }  
  
          printf("checking %d m/zs:\n", numpep);
    
          for (k = 0; k < numpep; k++)
          {

	    printf("%d m/zs...\n", k);
            /* read peptide's m/z (this is monoisotopic mass) and number of Qs */
            tmp = fscanf(fin, "%f %d ", &ftmp, &dtmp);
            pep[k].nq = dtmp;
            /* peptide should be observed at monoisotopic mass + 1 (for charge) */
            pep[k].pepmass = ftmp + 1.0;  
            pep[k].length = getLine(fin, pep[k].numseq, &zs, &pg, &cbm); 
            pep[k].zs = zs;
            pep[k].pg = pg;
            pep[k].cbm = cbm;
            /*printf("%f %d\n", pep[k].pepmass, pep[k].nq);*/
          }
          newnumpep = pepfiles(outfilename, data, nsamples, nreps, nmasses, pep, numpep, FITPEAKS, FIRSTMASS, SNRLIM);
          if (newnumpep > 0)
          {
            /* note: this will not be filled for those peptides with pep[i].found = 0 */
            ISODIST  *dist = NULL;
            dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));  
            for (k = 0; k < numpep; k++)
            {
              dist[k].mass = (double*) malloc (5 * sizeof(double));
              dist[k].prob = (double*) malloc (5 * sizeof(double));
            }
          
            iso(pep, numpep, dist);
          
            /* number of peaks in isotope distribution to be fitted */
            runga (outfilename, pep, numpep, dist, msamples, FITPEAKS, GALIM);
    
            for (k = 0; k < newnumpep; k++)
            {
              free (dist[k].mass);
              free (dist[k].prob);
            }
            free (dist);
    
          }
          else printf("none of these peptide found!\n");
      
          for (i = 0; i < nsamples; i++)
          {
            free (data[i].vars);
          }
          free (data);
      
          for (i = 0; i < numpep; i++)
          {
            for (j = 0; j < msamples; j++)
            {
              free (pep[i].sample[j].peaks);
            }
            free (pep[i].sample);
          }
          free (pep);
        }
      }
    }
  }
  return 0;
} 

