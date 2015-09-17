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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "constants.h"

void outputData(char *name, double *data, int nvars);
void outputSampleData(char *name, SAMPLE *data, int nsamples, int nvars);
void outputMasses(double *mzs, int nvars);
void getPeakIntensities(int jj, int mult, int npeaks, SAMPLE *data, int nsamples);
void fSort (double *array, int *index, int nRef);
void rmbase(int i, SAMPLE *data, int nmsz, double *fitted, int width);
void debug               (char *text, ...);
void fatalError          (char *text, ...);

int pepfiles (char *outfilename, SAMPLE *data, int nsamples, int nreps, int nmasses, PEPTIDE *pep, int numpep, int FITPEAKS, double FIRSTMASS, double SNRLIM)
{    
  int i, j, k, l;
  double diff, noise, max, totsnr, signal;
  int keepj;
           
  /* npeaks includes 2 noise peaks before distribution */
  int npeaks = FITPEAKS + 2;
  
  int msamples = nsamples/nreps;
  double snr[nreps];

  double peaks[npeaks];

  double *mzs = NULL;
  mzs = (double *) malloc (nmasses *sizeof(double));
  for (i = 0; i < nmasses; i++)
  {
    mzs[i] = FIRSTMASS + STEP*(float)i;
  }
    
  /* write out full data */
  char keepname[MAX_STRING];
  strcpy(keepname, outfilename);
  strcat(keepname, "Data\0");
  /* outputSampleData(keepname, data, nsamples, nmasses); */
  /* and masses 
  outputMasses(mzs, nmasses);*/

  /* remove baseline */
  double *fitted   = NULL;
  fitted   = (double *) malloc (nmasses * sizeof(double));
  for (i = 0; i < nsamples; i++)
  {
    rmbase(i, data, nmasses, fitted, 20);
    for (j = 0; j < nmasses; j++)
    {
      data[i].vars[j] = data[i].vars[j] - fitted[j];
      if (data[i].vars[j] < 0.0) data[i].vars[j] = 0.0;
    }
  }
  free (fitted);
  
  /* strcpy(keepname, outfilename);
  strcat(keepname, "Base\0");
  outputSampleData(keepname, data, nsamples, nmasses);*/

  for (i = 0; i < nsamples; i++)
  {
    data[i].peaks = (double *) malloc (npeaks * sizeof(double));
  }     
  int newnumpep = 0;  

  double mult = 1.0/STEP;
  int nmult = (int)mult;

  for (k = 0; k < numpep; k++)
  {
    /* use all peptides - can use this to make sure not missing for most samples */
    pep[k].found = 1;
    keepj = (int)(pep[k].pepmass*mult - (FIRSTMASS*mult-0.5));
    
    getPeakIntensities(keepj, nmult, npeaks, data, nsamples);
   
    for (i = 0; i < msamples; i++)
    {
      /* these will be zero if no replicates to be combined */
      pep[k].sample[i].noise = 0.0;
      pep[k].sample[i].signal = 0.0;
      /* first find out how to combine the replicates */
      for (l = 0; l < nreps; l++)
      {
        noise = (data[l + i*nreps].peaks[0] +data[l + i*nreps].peaks[1])/2.0;
        diff = fabs(data[l + i*nreps].peaks[0] - data[l + i*nreps].peaks[1]);
        max = 0.0;
        totsnr = 0.0;
        /* find maximum signal */
        for (j = 0; j < npeaks; j++)
        {          
          if (data[l + i*nreps].peaks[j] > max) max = data[l + i*nreps].peaks[j];
        }
        /* calculate signal to noise if noise difference < smallest noise peak or < 10* signal, else snr for rep is zero */
        snr[l] = 0.0;
        if ((diff + 1.0 < data[l + i*nreps].peaks[0] + 1.0) && (diff + 1.0 < data[l + i*nreps].peaks[1] + 1.0)) snr[l] = max - noise;
        else if (diff + 1.0 < 10.0*max + 1.0) snr[l] = max - noise;
        totsnr = totsnr + snr[l];
      }
      /* now combine replicates according to snr */
      for (j = 0; j < npeaks; j++)
      {          
        peaks[j] = 0.0;
        for (l = 0; l < nreps; l++)
        {
          peaks[j] = peaks[j] + snr[l]*data[l + i*nreps].peaks[j]/totsnr;
       }
      }
      /* only write out if combined signal to noise is good enough */  
      noise = (peaks[0] + peaks[1])/2.0;
      signal = peaks[2];
      /*just check first 3 peaks for signal intensity */
      for (j = 2; j < 5; j++)
      {          
        if (peaks[j] > signal) signal = peaks[j];
      }
      pep[k].sample[i].signal = signal;
      pep[k].sample[i].noise = noise;
      pep[k].sample[i].found = 0;
      for (j = 0; j < npeaks-2; j++)
      {     
        pep[k].sample[i].peaks[j] = 0.0;
      }
      if (signal + 1.0 > SNRLIM*(noise+1.0))
      {
        pep[k].sample[i].found = 1;
        pep[k].sample[i].diff = fabs(peaks[0] - peaks[1]);
        /* only keep the FITPEAKS peaks from monoisotopic mass */
        for (j = 0; j < npeaks-2; j++)
        {     
          pep[k].sample[i].peaks[j] = peaks[j+2];
        }
      }
    }
    newnumpep++;
  }
  free (mzs);
  
  return (newnumpep);
}


/*************************************************************************************************************
Procedure:  outputData
outs a single spectra to "name"
************************************************************************************************************/
void outputData(char *name, double *data, int nvars)
{
  int j;
  FILE *fileHandle = NULL;

  fileHandle = NULL;
  fileHandle = fopen(name, "w");   

  for (j = 0; j < nvars; j++)
  {
    fprintf(fileHandle, "%f ", data[j] ); 
  }
  fprintf(fileHandle, "\n");     
  fclose(fileHandle);      
}


/*************************************************************************************************************
Procedure:  getPeakIntensities
************************************************************************************************************/
void getPeakIntensities(int jj, int diff, int npeaks, SAMPLE *data, int nsamples)
{
   int i, j, k;
   int start, keepj;
   double max;
   
   for (i = 0; i < nsamples; i++)
   {
     /* find peak maxima */
     max = data[i].vars[jj];
     keepj = jj;
     for (k = jj-4; k < jj+5; k++)
     {
       if (data[i].vars[k] > max) 
       {
         max = data[i].vars[k];
         keepj = k;
       }
     }  
     
     /* start two peaks before distribution */
     start = keepj - 2*diff;
     for (j = 0; j < npeaks; j++)
     {
       data[i].peaks[j] = 0.0;
       /* use 3 data points to get peak intensity */
       for (k = start-1; k < start+2; k++)
       {
         data[i].peaks[j] += data[i].vars[k];
       }
       start = start + diff;
     }
   }

}


/****************************************************************
Procedure:  fSort
****************************************************************/
void fSort (double *array, int *index, int nRef)
{
  int nInt, ii, i, j, nf;
  double tf;
  nInt = nRef;
  
  do
  {
    nInt = nInt / 2;
    if ((nInt % 2) == 0)
    {
      nInt = nInt-1;
    }
    for (ii= 0; ii < nRef - nInt; ii++)
    {
      i = ii;
      j = i + nInt;
      if (array[i] > array[j])
      {
        tf = array[j];
        nf = index[j];
        do
        {
          array[j] = array[i];
          index[j] = index[i];
          j = i;
          i = i - nInt;
        }
        while ((i >= 0) && (array[i] > tf));
        array[j] = tf;
        index[j] = nf;
      }
    }
  }
  while (nInt > 1);
}

/*************************************************************************************************************
Procedure:  outputSampleData
************************************************************************************************************/
void outputSampleData(char *name, SAMPLE *data, int nsamples, int nvars)
{
  int i,j;   
  FILE *fileHandle = NULL;
  fileHandle = fopen(name, "w");   

  for (i = 0; i < nsamples ; i++)
  {
    for (j = 0; j < nvars; j++)
    {
      fprintf(fileHandle, "%f  ", data[i].vars[j]);    
    }
    fprintf(fileHandle, "\n");     
  }
  fclose(fileHandle);
}  
    
/*************************************************************************************************************
Procedure:  outputMasses
************************************************************************************************************/
void outputMasses(double *mzs, int nvars)
{
  int j;   

  FILE *fp = NULL;
  fp = fopen("mzs", "w");   

  for (j = 0; j < nvars; j++)
  {
     fprintf(fp, "%f ", mzs[j]);
  }
  fprintf(fp, "\n");
  fclose(fp);
}      
 
/*************************************************************************************************************
Procedure: smooth
performs moving average smoothing over w points
************************************************************************************************************/
void smooth(double *indata, int nmzs, double *outdata, int w) 
{
  int j, k, n;
  double av;
  
  for (j = 0; j < w; j++)
  {
    av = 0.0;
    n = 0;
    for (k = 0; k < j+w; k++)
    {
      av += indata[k];
      n++;
    }
    outdata[j] = av/(float)n;
  }

  for (j = w; j < nmzs-w-1; j++)
  {
    av = 0.0;
    n = 0;
    for (k = j-w; k < j+w; k++)
    {
      av += indata[k];
      n++;
    }
    outdata[j] = av/(float)n;
  }

  for (j = nmzs-w-1; j < nmzs; j++)
  {
    av = 0.0;
    n = 0;
    for (k = j-w; k < nmzs; k++)
    {
      av += indata[k];
      n++;
    }
    outdata[j] = av/(float)n;
  }
}

/*************************************************************************************************************
Procedure: newbase
finds the baseline estimate from the local minima
************************************************************************************************************/
void rmbase(int i, SAMPLE *data, int nmzs, double *fitted, int w) 
{    
  int j, k;
  double min;
  
  double temp[nmzs];
  
  for (j = 10; j < w; j++)
  {
    min = data[i].vars[j];
    temp[j] = min;
    for (k = 0; k < j+w; k++)
    {
      if (data[i].vars[k] < min) 
      {
        temp[j] = data[i].vars[k];
        min = temp[j];
      }
    }
  }

  for (j = w; j < nmzs-w-1; j++)
  {
    min = data[i].vars[j];
    temp[j] = min;
    for (k = j-w; k < j+w; k++)
    {
      if (data[i].vars[k] < min) 
      {
        temp[j] = data[i].vars[k];
        min = temp[j];
      }
    }
  }

  for (j = nmzs-w-1; j < nmzs; j++)
  {
    min = data[i].vars[j];
    temp[j] = min;
    for (k = j-w; k < nmzs; k++)
    {
      if (data[i].vars[k] < min) 
      {
        temp[j] = data[i].vars[k];
        min = temp[j];
      }
    }
  }
  
 /* now smooth to remove odd jumps */
  smooth(temp, nmzs, fitted, 2*w); 
    
}



