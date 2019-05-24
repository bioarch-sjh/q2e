#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "constants.h"

typedef struct allData {
  double *peak;
} DATA;

typedef struct allAlphas {
  double *pep;
} ALFAS;

double ga(double *peaks, ISODIST  *dist, int kk, int np, int ng, double *alpha, int nc, int nal, int FITPEAKS);
void testfit1(double *peaks, ISODIST  *dist, int kk, int npeaks, double alpha1);
void testbeta2(double *peaks, ISODIST  *dist, int kk, int npeaks, double beta1, double beta2);
void testbeta3(double *peaks, ISODIST  *dist, int kk, int npeaks, double beta1, double beta2, double beta3);
void fSort (double *array, int *index, int nRef);

/****************************************
Procedure:      runga
*****************************************/
void runga (char *outfilename, PEPTIDE *pep, int numpep, ISODIST  *dist, int msamples, int FITPEAKS, double GALIM)
{  
  FILE *fout = NULL;
  FILE *fout2 = NULL;

  /* Number of generations  */
  int ngen = 100;
  /* Number in the population  */
  int npop = 2000;    
  /* number of ga runs */
  int nruns = 10; 
  
  int nc = 10;
      
  int i, j, k, l;
  double sum, fit, bestfit;
  double noise, diff, signal;
  	
  printf("fitting %d peaks\n", FITPEAKS);
  double *peaks = NULL;
  peaks = (double *) malloc (FITPEAKS * sizeof(double));

  char keepname[MAX_STRING];
  strcpy(keepname, outfilename);
  strcat(keepname, "Results\0");
  fout = fopen(keepname, "w");

  strcpy(keepname, outfilename);
  strcat(keepname, "Betas.csv\0");
  fout2 = fopen(keepname, "w");

  int countminus1; 
  double last;

  for (k = 0; k < numpep; k++)
  {
    countminus1 = 0;
    double alpha[pep[k].nq];
    double bestalpha[pep[k].nq];
    double meanalpha[pep[k].nq];
    double stderralpha[pep[k].nq];
    for (l = 0; l < pep[k].nq; l++)
    {
      meanalpha[l] = 0.0;
      stderralpha[l] = 0.0;
    }
    if (pep[k].found == 1)
    {
      printf("%0.1f \n", pep[k].pepmass);
      fprintf(fout, "%0.1f\n", pep[k].pepmass);

      /* rescale the isotope distribution */
      sum = 0.0;
      for (i = 0; i < FITPEAKS; i++)
      {
        sum = sum + dist[k].prob[i];
      }
      /* rescale isotope values to give a percentage */
      for (i = 0; i < FITPEAKS; i++)
      {
        dist[k].prob[i] = dist[k].prob[i]/sum;
      }
      
      for (i = 0; i < msamples; i++)
      {
        for (l = 0; l < pep[k].nq; l++)
        {
          pep[k].sample[i].betas[l] = -1.0;
        }
        signal = 0.0; 
        /* check whether this peptide was seen in this sample */
        if (pep[k].sample[i].found == 1)
        {
          noise = pep[k].sample[i].noise;   
          diff = pep[k].sample[i].diff;   
          for (j = 0; j < FITPEAKS; j++)
          {     
             peaks[j] = pep[k].sample[i].peaks[j];
          }    
          sum = 0.0;
          /* get signal intensity from first 3 peaks and subtract the noise from all peaks */
          for (j = 0; j < FITPEAKS; j++)
          {     
            if ((j < 4) && (peaks[j] > signal)) signal = peaks[j];
            peaks[j] = peaks[j] - noise;
            sum = sum + peaks[j];          
          }
          /* rescale experimental peaks to give a percentage */
          for (j = 0; j < FITPEAKS; j++)
          {     
            peaks[j] = peaks[j]/sum;
          }              
            
          for (j = 0; j < nruns; j++)
          {
            fit = ga(peaks, dist, k, npop, ngen, alpha, nc, pep[k].nq, FITPEAKS);
            if ((j == 0) || (fit < bestfit)) 
            {
              bestfit = fit;
              for (l = 0; l < pep[k].nq; l++)
              {
                bestalpha[l] = alpha[l];
              }
            } 
          }
          if (bestfit < GALIM) 
          {
            printf("sample %d: best fit is : %f with signal to noise level:  %0.2f    noise diff:  %0.2f  ",  i+1, bestfit, signal/noise, diff/noise); 
            fprintf(fout,"sample %d: best fit is : %f with signal to noise level:  %0.2f    noise diff:  %0.2f  ",  i+1, bestfit, signal/noise, diff/noise);            
            last = 1.0;
            for (l = 0; l < pep[k].nq; l++)
            {
              meanalpha[l] += bestalpha[l];
              stderralpha[l] += bestalpha[l]*bestalpha[l];          
              printf("beta%d = %0.3f  ", l+1, bestalpha[l]);
              fprintf(fout, "beta%d = %0.3f  ", l+1, bestalpha[l]);
              pep[k].sample[i].betas[l] = bestalpha[l];         
              last = last - bestalpha[l];
            }
            printf("beta%d %0.3f\n", pep[k].nq + 1, last);
            fprintf(fout,"\n");            
          }
          else
          {
            printf("sample %d: rejected due to poor fit in GA\n", i+1);
            fprintf(fout,"sample %d: rejected due to poor fit in GA\n", i+1);
            countminus1++;
          }
        }
        else
        {
          printf("sample %d: rejected due to poor SNR: signal %0.2f; noise %0.2f\n", i+1, pep[k].sample[i].signal, pep[k].sample[i].noise);
          fprintf(fout, "sample %d: rejected due to poor SNR: signal %0.2f; noise %0.2f\n", i+1, pep[k].sample[i].signal, pep[k].sample[i].noise);
          countminus1++;
        }
      }
    }
    else
    {
      for (i = 0; i < msamples; i++)
      {
        countminus1++;
      }
    }
    printf("peptide %d: %0.1f number of samples that beta values are calculated for is %d.\n", k, pep[k].pepmass, msamples - countminus1);
   }
  fclose(fout);

  /* output betas for each peptide (columns) for each sample (rows) 
  For multiple Qs, betas are separated by semi-colons */
  fprintf(fout2, "sample name,");
  for (k = 0; k < numpep; k++)
  {
    fprintf(fout2, "%0.2f,", pep[k].pepmass);
  }
  fprintf(fout2, "\n");
  for (j = 0; j < msamples; j++)
  {
    fprintf(fout2, "%s,", pep[0].sample[j].name);
    for (k = 0; k < numpep; k++)
    {
      for (l = 0; l < pep[k].nq; l++)
      {
        if (pep[k].sample[j].betas[l] > -0.5) fprintf(fout2, "%0.2f", pep[k].sample[j].betas[l]);
        if ((pep[k].nq > 1) && (l < pep[k].nq-1)) fprintf(fout2, "; ");
      }
      fprintf(fout2, ",");
    }
    fprintf(fout2, "\n");
  }
  fclose(fout2);

  free(peaks);
}

/*************************************************************************************************************
Procedure:  evaluatebetaFitness
function to evaluate the fitness ... small is best!
*************************************************************************************************************/
double evaluatebetaFitness(double *peaks, ISODIST  *dist, int kk, int npeaks, int *c, int nc)
{
  int i, digit;
 
  double power, beta;
  double eq, error;
  
  beta = 0.0;
  for (i = 0; i < nc; i++)
  { 
     digit = c[i];
     power = i + 1;
     beta = beta + (float)(digit)/(float)pow(10, power); 
  }
  
  error = (peaks[0] - beta*dist[kk].prob[0])*(peaks[0] - beta*dist[kk].prob[0]);
  for (i = 1; i < npeaks; i++)
  {
     eq = beta * dist[kk].prob[i] + (1.0 - beta) * dist[kk].prob[i-1]; 
     error = error + (peaks[i] - eq)*(peaks[i] - eq);
  }
  
  return (error);	
} 

/*************************************************************************************************************
Procedure:  evaluatebetaFitness2
function to evaluate the fitness ... small is best!
*************************************************************************************************************/
double evaluatebetaFitness2(double *peaks, ISODIST  *dist, int kk, int npeaks, int *c, int nc)
{
  int i, j, digit;
  int n = nc/2;
 
  double power;
  double eq, error;
  
  double beta[2];
  
  for (j = 0; j < 2; j++)
  {
    beta[j] = 0.0;
    for (i = 0; i < n; i++)
    { 
       digit = c[j*n + i];
       power = i + 1;   
       beta[j] = beta[j] + (float)(digit)/(float)pow(10, power); 
    }
  }
  
  double gamma1, gamma2;
  gamma1 = beta[0];
  gamma2 = beta[1]*(1.0 - gamma1);

  /* first peak */
  eq = gamma1*dist[kk].prob[0];
  error = (peaks[0] - eq) * (peaks[0] - eq);
  /* second peak */
  eq = gamma1*dist[kk].prob[1] + gamma2*dist[kk].prob[0];
  error += (peaks[1] - eq) * (peaks[1] - eq);
  /* other peaks */
  for (i = 2; i < npeaks; i++)
  {
    eq =  gamma1*dist[kk].prob[i] + gamma2*dist[kk].prob[i-1] + (1.0 - gamma1 - gamma2)*dist[kk].prob[i-2];
    error += (peaks[i] - eq)*(peaks[i] - eq);
  }
  
  return (error);	
} 

/*************************************************************************************************************
Procedure:  evaluatebetaFitness3
function to evaluate the fitness ... small is best!
*************************************************************************************************************/
double evaluatebetaFitness3(double *peaks, ISODIST  *dist, int kk, int npeaks, int *c, int nc)
{
  int i, j, digit;
  int n = nc/3;
 
  double power;
  double eq, error;
  
  double beta[3];
  
  for (j = 0; j < 3; j++)
  {
    beta[j] = 0.0;
    for (i = 0; i < n; i++)
    { 
       digit = c[j*n + i];
       power = i + 1;   
       beta[j] = beta[j] + (float)(digit)/(float)pow(10, power); 
    }
  }
  
  double gamma1, gamma2, gamma3;
  gamma1 = beta[0];
  gamma2 = beta[1]*(1.0 - gamma1);
  gamma3 = beta[2]*(1.0 - gamma1 - gamma2);

  /* first peak */
  eq = gamma1*dist[kk].prob[0];
  error = (peaks[0] - eq) * (peaks[0] - eq);
  /* second peak */
  eq =  gamma1*dist[kk].prob[1] + gamma2*dist[kk].prob[0];
  error += (peaks[1] - eq) * (peaks[1] - eq);
  /* third peak */
  eq =  gamma1*dist[kk].prob[2] + gamma2*dist[kk].prob[1] + gamma3*dist[kk].prob[0];
  error += (peaks[2] - eq) * (peaks[2] - eq);
  /* other peaks */
  for (i = 3; i < npeaks; i++)
  {
    eq =  gamma1*dist[kk].prob[i] + gamma2*dist[kk].prob[i-1] + gamma3*dist[kk].prob[i-2] + (1.0 - gamma1 - gamma2 - gamma3)*dist[kk].prob[i-3];
    error += (peaks[i] - eq)*(peaks[i] - eq);
  }

 return (error);	
} 

/*************************************************************************************************************
   Procedure:  ga
   Number in the population  np 
   length of chromosomes  nc  
   Number of generations (if huge leaves running until satisfactory) ng
************************************************************************************************************/
double ga(double *peaks, ISODIST  *dist, int kk, int np, int ng, double *alpha, int nc, int nal, int FITPEAKS)
{  
   nc = nc*nal; 
  /* Percentage crossover ie probablility of performing crossover */
  int pcrss = 90;  
  /* Percentage mutation ie probablility of performing mutation */
  int mut = 40;    
  /* Number of offspring to reproduce at each population-must be even number
  Not reproducing entire population each time ie steady state algorithm 
  must be less than np */
  int numoffspring = 50*np/100; 
  /* fitness for each member of population */
  double fitness[np+numoffspring];
  /* arrays for entire populations */
  int *a = NULL;
  a = (int *) malloc (nc * (np + numoffspring) * sizeof(int));
  int *b = NULL;
  b = (int *) malloc (nc * numoffspring * sizeof(int));
  int *c = NULL;
  c = (int *) malloc (nc * sizeof(int));

  int ind, ind1, ind2, jnd, n, i;  
  int cc, cp;  /* chromosome counters */
  int cg;  /* generation counter */
  int gap, si, sj;  /* used by sort routine */
  float temp;     
  int random1, random2;
  int parent1, parent2;
  int tmp, crss, swap;            
  int range, power, digit;
  int npeaks = FITPEAKS;

  /* seed random number generator */
  time_t tt;
  tt = time(NULL);
  srand(tt);
    
  /* produce inital population */
  for (cp = 0; cp < np; cp++)
  {
    for (cc = 0; cc < nc; cc++)
    {
      /* random integer is between 0 and range */
      range = 10;
      ind = cc + nc*cp;
      a[ind] = rand()%range;
      c[cc] = a[ind]; 
    }    
    if (nal == 1) fitness[cp] = evaluatebetaFitness(peaks, dist, kk, npeaks, c, nc);
    else if (nal == 2) fitness[cp] = evaluatebetaFitness2(peaks, dist, kk, npeaks, c, nc);
    else fitness[cp] = evaluatebetaFitness3(peaks, dist, kk, npeaks, c, nc);
  }

  /* This main loop going through each generation */
  cg = 1;
  while (cg <= ng)
  {  
    /* sorts population  in order of fitness */
    for (gap = (np)/2; gap > 0; gap /= 2)
    {
      for (si = gap; si < np; si++)
      {  
         for (sj = si-gap; sj>=0 && fitness[sj] > fitness[sj+gap]; sj -= gap)
         {
           temp = fitness[sj];
           fitness[sj] = fitness[sj+gap];
           fitness[sj+gap] = temp;
           for (cc = 0; cc < nc; cc++)
           {
             ind = cc + nc*sj;
             ind1 = cc + nc*(sj+gap);
             tmp = a[ind];
             a[ind] = a[ind1];
             a[ind1] = tmp;
           }
         }
       }
     }     

     if (cg == ng)
     {
       double beta[3]; 
       for (i = 0; i < nal; i++)
       {
          beta[i] = 0.0;
       }
       n = nc/nal;
       for (i = 0; i < nal; i++)
       {
         for (cc = 0; cc < n; cc++)
         { 
           digit = a[n*i+cc];
           power = cc + 1;
           beta[i] += (float)(digit)/(float)pow(10, power);
         }
       }

       /* make sure all betas are in range 0 to 1 */
       alpha[0] = beta[0];
       if (nal > 1) 
       {
         alpha[1] = beta[1]*(1.0 - alpha[0]);
         if (nal > 2) 
         {
           alpha[2] = beta[2]*(1.0 - alpha[0] - alpha[1]);
         }
       }
     }
     else
     {
       /*  produce offspring */
       for(cp = 0; cp < numoffspring; cp += 2)
       {
         /* chooses one parent */
         random1 = rand()%np;
         random2 = rand()%np;
         if (fitness[random1] <= fitness[random2]) {parent1 = random1;}
         else  {parent1 = random2;}
       
         /* chooses another parent  */
         random1 = rand()%np;
         random2 = rand()%np;
         if (fitness[random1] <= fitness[random2]) {parent2 = random1;}
          else {parent2=random2;}
       
        /* crossover if random number less than pcrss of producing offsprings */     
         swap = rand()%100+1;
         if (swap < pcrss)
         {         
           for (cc = 0; cc < nc; cc++)
           {
             swap = rand()%(nc-1) + 1;
             ind = cc + nc*cp;
             jnd = cc + nc*(cp+1);
             ind1 = cc + nc*parent1;
             ind2 = cc + nc*parent2;
             if (cc < swap)
             {
               b[ind] = a[ind1];
               b[jnd] = a[ind2];
             }
             else
             {
               b[ind] = a[ind2];
               b[jnd] = a[ind1];
             }
           } 
         } 
         else
         {
           for (cc = 0; cc < nc; cc++)
           {
             ind = cc + nc*cp;
             jnd = cc + nc*(cp+1);
             ind1 = cc + nc*parent1;
             ind2 = cc + nc*parent2;
             b[ind] = a[ind1];
             b[jnd] = a[ind2];
           }
         }
       }
     
       /* add offspring to extended population */
       for (cp = 0; cp < numoffspring; cp++)
       {
         for (cc = 0; cc < nc; cc++)
         {
           ind = cc + nc*(cp + np);
           jnd = cc + nc*cp;
           a[ind] = b[jnd];
         }
       }
  
       /**** mutate new offspring ****/
       for (cp = np; cp < np+numoffspring; cp++)
       { 
         /* only mutate if if random number less than mut */
         swap = rand()%100+1;
         if (swap < mut)
         {
           /* choose randon gene in chromosome to mutate */
           swap = rand()%(nc-1);
           ind = swap + nc*cp;
           /* new random integer is between 0 and range */
           range = 10;
           crss = rand()%range;
           a[ind] = crss;
         }
       }    
   
       /* evaluate the new offspring */
       for (cp = np; cp < np+numoffspring; cp++)
       {
         for (cc = 0; cc < nc; cc++)
         {
           ind = cc + nc*cp;
           c[cc] = a[ind];
         }
         if (nal == 1) fitness[cp] = evaluatebetaFitness(peaks, dist, kk, npeaks, c, nc);
         else if (nal == 2) fitness[cp] = evaluatebetaFitness2(peaks, dist, kk, npeaks, c, nc);
         else fitness[cp] = evaluatebetaFitness3(peaks, dist, kk, npeaks, c, nc);
       }
      
       /* add the offspring to the population */
       for (cp = np-numoffspring; cp < np; cp++)
       {
         for (cc = 0; cc < nc; cc++)
         {
           ind = cc + nc*cp;
           ind1 = cc + nc*(cp+numoffspring);
           a[ind] = a[ind1];
         }
         fitness[cp] = fitness[cp+numoffspring];
       }    
     }
     cg++;
  } /*end of loop through generations */
   
  free(a);
  free(b);
  free(c);
     
  return (fitness[0]);  
} 
