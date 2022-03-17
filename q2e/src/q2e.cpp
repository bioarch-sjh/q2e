#include <Rcpp.h>
using namespace Rcpp;


#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <chrono>

#include "constants.h"

#include "makepeptidefiles.h"
#include "isodists.h"
#include "ga.h"


//ERROR CODES RETURNED BY RQ2E
#define Q2E_BAD_ARGS                     1
#define Q2E_CANT_OPEN_PARAMFILE          2
#define Q2E_CANT_OPEN_PEPTIDEFILE        3
#define Q2E_CANT_OPEN_DATAFILE           4
#define Q2E_BAD_NFILES_NREPLICATES_COMBO 5
#define Q2E_BAD_LINE_IN_DATAFILE         6

//#define DBG_PRINT

//' Multiply a number by three
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesThree(int x) {
  return x * 3;
}



/*C Procedures*/
//int pepfiles (char *name, SAMPLE *data, int nsamples, int nreps, int nmasses, PEPTIDE *pep, int numpep, int FITPEAKS, double FIRSTMASS, double SNRLIM);
//void iso(PEPTIDE *pep, int numpep, ISODIST  *dist);
int getLine(FILE *fp, int v[], int *zs, int *pg, int *cbm);
//void runga (char *name, PEPTIDE *pep, int numpep, ISODIST  *dist, int msamples, int FITPEAKS, double GALIM);




/* Debugging the pep structure via printfs..

 pep[k].pepmass = ftmp + 1.0;
 pep[k].length = getLine(fin, pep[k].numseq, &zs, &pg, &cbm);
 pep[k].zs = zs;
 pep[k].pg = pg;
 pep[k].cbm = cbm;*/
void printpep(PEPTIDE *pep, int numpep){

  Rprintf("%d peptides in pep array\n",numpep);
  Rprintf("#\tnq\tmass\t\tlen\tzs\tpg\tcbm\n");
  for(int pp=0;pp<numpep;pp++){
    Rprintf("%d\t%d\t%4.3f\t%d\t%d\t%d\t%d\n",
            pp,
            pep[pp].nq,
            pep[pp].pepmass,
            pep[pp].length,
            pep[pp].zs,
            pep[pp].pg,
            pep[pp].cbm);
  }
}


/****************************************
 Procedure:      rq2e, access point for the
 *****************************************/
//' Run q2e via datafiles as in the C++ code
//'
//' @param paramFile the parameter file
// [[Rcpp::export]]
int rq2e ( Rcpp::StringVector argVector, Rcpp::StringVector fnVector, int Rnreplicates)
{
  FILE *fp_param = NULL;
  //FILE *fp = NULL;
  FILE *fin = NULL;
  int i, j, k, ii, nsamples, nreps, zs, cbm, pg;
  int numpep, newnumpep;
  float ftmp, ftmp1;
  int tmp, dtmp;//, dtmp1;
  char name[MAX_STRING];

  int FITPEAKS;
  double GALIM;
  double SNRLIM;
  double FIRSTMASS;
  double LASTMASS;
  int CSV;

  int nargs = argVector.size();

  //* Check that we have the correct number of files */
  if(nargs != 4){
    Rprintf("ERROR Q2E_BAD_ARGS: Correct program useage: ./QtoEcalc [data filelist] [peptide filelist] [parameter file] \n");
    return Q2E_BAD_ARGS;
  }

  std::string fileListFile = Rcpp::as<std::string>(argVector[1]);
  std::string peptideFile = Rcpp::as<std::string>(argVector[2]);
  std::string paramFile = Rcpp::as<std::string>(argVector[3]);

  Rprintf("reading parameters from %s\n", paramFile.c_str());
  fp_param = fopen(paramFile.c_str(), "r");
  if (fp_param == NULL)
  {
    Rprintf("ERROR Q2E_CANT_OPEN_PARAMFILE: couldn't open parameter file %s\n", paramFile.c_str());
    return Q2E_CANT_OPEN_PARAMFILE;
  }
  else
  {
    tmp = fscanf(fp_param, "%s %d\n", name, &dtmp);
    FITPEAKS = dtmp;
    Rprintf("%s = %d\n",name, FITPEAKS);
    tmp = fscanf(fp_param, "%s %f\n", name, &ftmp);
    GALIM = ftmp;
    Rprintf("%s = %0.2f\n",name, GALIM);
    tmp = fscanf(fp_param, "%s %f\n", name, &ftmp);
    SNRLIM = ftmp;
    Rprintf("%s = %0.1f\n",name, SNRLIM);
    tmp = fscanf(fp_param, "%s %f\n", name, &ftmp);
    FIRSTMASS = ftmp;
    Rprintf("%s = %0.1f\n",name, FIRSTMASS);
    tmp = fscanf(fp_param, "%s %f\n", name, &ftmp);
    LASTMASS = ftmp;
    Rprintf("%s = %0.1f\n",name, LASTMASS);
    tmp = fscanf(fp_param, "%s %d\n", name, &dtmp);
    CSV = dtmp;
    Rprintf("%s = %d\n",name, CSV);
    fclose(fp_param);
  }



  /* keepname for output files */
  char outfilename[1024];
  strcpy(outfilename, fileListFile.c_str());

  /*
   Rprintf("reading data files from %s\n", fileListFile.c_str());
   fp = fopen(outfilename, "r");
   if (fp == NULL)
   {
   Rprintf("ERROR Q2E_CANT_OPEN_FILELIST: something is wrong with the filelist, %s\n", fileListFile.c_str());
   return Q2E_CANT_OPEN_FILELIST;
   }*/

  //Check that the number of files is a multiple of the number of replicates:
  if(fnVector.size()%Rnreplicates != 0){

    Rprintf("ERROR Q2E_BAD_NFILES_NREPLICATES_COMBO: nfiles = %d, nreplicates = %d\n", fnVector.size(), Rnreplicates);
    return Q2E_BAD_NFILES_NREPLICATES_COMBO;
  }
  else
  {
    /*
     tmp = fscanf(fp, "%d %d\n", &dtmp, &dtmp1);
     nsamples = dtmp;
     nreps = dtmp1;
     int msamples = nsamples/nreps;
     */

    nsamples = fnVector.size();
    nreps = Rnreplicates;
    int msamples = nsamples/nreps;

    Rprintf("%d individual samples with %d replicates each of %d samples\n", nsamples, nreps, msamples);

    // number of masses to be kept: nmasses = lastmass - firstmass in step m/z units
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

    // read in filenames
    for (i = 0; i < nsamples; i++)
    {
      //tmp = fscanf(fp, "%s\n", name);
      std::string fileName = Rcpp::as<std::string>(fnVector[i]);
      Rprintf("Loading data from file %s (entry %d in the file list)\n",fileName.c_str(),i);
      fin = fopen(fileName.c_str(), "r");
      if (fin == NULL)
      {
        Rprintf("ERROR Q2E_CANT_OPEN_DATAFILE: couldn't find data file, %s (entry %d in the file list).\n", fileName.c_str(),i);
        return(Q2E_CANT_OPEN_DATAFILE );
      }
      else
      {
        Rprintf("%s, CSV = %d; ftmp = %f, ftmp1 = %f\n", fileName.c_str(), CSV, ftmp, ftmp1);
        strcpy(data[i].name, fileName.c_str());
        // read in data
        tmp = 0;
        int line = 0;
        while (tmp != EOF)
        {
          line++;

          // read in mass
          if (CSV == 1){
            tmp = fscanf(fin, "%f,%f\n", &ftmp, &ftmp1);
          }
          else          {
            tmp = fscanf(fin, "%f %d\n", &ftmp, &dtmp);
            ftmp1 = (float)dtmp;
          }

          if(tmp != 2 && tmp != EOF){
            Rprintf("ERROR Q2E_BAD_LINE_IN_DATAFILE: Line = %d, Read %d arguments: ftmp = %f\tftpmp1 = %f\n",line,tmp,ftmp,ftmp1);
            return(Q2E_BAD_LINE_IN_DATAFILE);
          }

          if (tmp != EOF){
            ii = (int)((ftmp - FIRSTMASS)*mult + 0.5);
            if ((ii > 0) && (ii < nmasses))
            {
              data[i].vars[ii] += ftmp1;
              data[i].num[ii]++;
            }
          }
        }
        fclose(fin);
        Rprintf("Successfully loaded data from file %s (entry %d in the file list)\n",fileName.c_str(),i);
      }
    }
    //fclose(fp);


    for (i = 0; i < nsamples; i++)
    {
      for (k = 0; k < nmasses; k++)
      {
        if (data[i].num[k] > 0) data[i].vars[k] /= (float) data[i].num[k];
      }
    }

    // read number of peptides
    fin = fopen(peptideFile.c_str(), "r");
    if (fin == NULL)
    {
      Rprintf("ERROR %d: couldn't open file with list of peptides, %s\n",Q2E_CANT_OPEN_PEPTIDEFILE, peptideFile.c_str());
      return(Q2E_CANT_OPEN_PEPTIDEFILE);
    }
    else
    {
      tmp = fscanf(fin, "%d\n", &dtmp);
      numpep = dtmp;
      PEPTIDE *pep = NULL;
      pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));

      for (k = 0; k < numpep; k++)
      {
        pep[k].sample = (SAMPLEPEP*) malloc (msamples * sizeof(SAMPLEPEP));
        pep[k].numseq = (int*) malloc (MAXLINE * sizeof(int));
      }

#ifdef DBG_PRINT
      Rprintf("checking %d m/zs:\n", numpep);
#endif

      for (k = 0; k < numpep; k++)
      {
        // read peptide's m/z (this is monoisotopic mass) and number of Qs
        tmp = fscanf(fin, "%f %d ", &ftmp, &dtmp);
        pep[k].nq = dtmp;
        // peptide should be observed at monoisotopic mass + 1 (for charge)
        pep[k].pepmass = ftmp + 1.0;
        pep[k].length = getLine(fin, pep[k].numseq, &zs, &pg, &cbm);
        pep[k].zs = zs;
        pep[k].pg = pg;
        pep[k].cbm = cbm;
#ifdef DBG_PRINT
        Rprintf("\t%f nq:%d  zs:%d\n", pep[k].pepmass, pep[k].nq, pep[k].zs);
#endif
        for (j = 0; j < msamples; j++)
        {
          pep[k].sample[j].peaks = (double*) malloc ((FITPEAKS + 2) * sizeof(double));
          pep[k].sample[j].betas = (double*) malloc (pep[k].nq * sizeof(double));
          pep[k].sample[j].found = -1;
          strcpy(pep[k].sample[j].name, data[j].name);
        }

      }
#ifdef DBG_PRINT
      Rprintf("Finished checking peptideList\n");
      printpep(pep,numpep);
#endif

      newnumpep = pepfiles(outfilename, data, nsamples, nreps, nmasses, pep, numpep, FITPEAKS, FIRSTMASS, SNRLIM);


      if (newnumpep > 0)
      {
        // note: this will not be filled for those peptides with pep[i].found = 0
        ISODIST  *dist = NULL;
        dist = (ISODIST*) malloc (numpep * sizeof(ISODIST));
        for (k = 0; k < numpep; k++)
        {
          dist[k].mass = (double*) malloc (5 * sizeof(double));
          dist[k].prob = (double*) malloc (5 * sizeof(double));
        }

#ifdef DBG_PRINT
        Rprintf("Running iso()\n");
#endif
        iso(pep, numpep, dist);

        // number of peaks in isotope distribution to be fitted
#ifdef DBG_PRINT
        Rprintf("Running runga()\n");
#endif
        runga (outfilename, pep, numpep, dist, msamples, FITPEAKS, GALIM);

        for (k = 0; k < newnumpep; k++)
        {
          free (dist[k].mass);
          free (dist[k].prob);
        }
        free (dist);

      }
      else Rprintf("none of these peptide found!\n");

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
          free (pep[i].sample[j].betas);
        }
        free (pep[i].sample);
      }
      free (pep);

    }

  }

  return 0;

}










/****************************************
 Procedure:      main
 *****************************************/
int q2emain (int argc, char *argv[])
{

  /* Check that we have the correct number of files */
  if (argc != 4 ) printf("ERROR: Correct program useage: ./QtoEcalc [data filelist] [peptide filelist] [parameter file] \n");
  else
  {
    FILE *fp = NULL;
    FILE *fin = NULL;
    int i, j, k, ii, nsamples, nreps, zs, cbm, pg;
    int numpep, newnumpep;
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
    if (fp == NULL)
    {
      printf("ERROR: couldn't open parameter file %s\n", argv[3]);
      exit(-1);
    }
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

    /* keepname for output files */
    char outfilename[1024];
    strcpy(outfilename, "RQ2E_OUTFILE.txt");//argv[1]);

    printf("reading data files from %s\n", argv[1]);
    fp = fopen(outfilename, "r");
    if (fp == NULL)
    {
      printf("ERROR: something is wrong with the filelist, %s\n", argv[1]);
      exit(-1);
    }
    else
    {
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
          printf("ERROR: couldn't find file, %s. Names may not match those in the filelist \n", name);
          exit(-1);
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

      for (i = 0; i < nsamples; i++)
      {
        for (k = 0; k < nmasses; k++)
        {
          if (data[i].num[k] > 0) data[i].vars[k] /= (float) data[i].num[k];
        }
      }

      /* read number of peptides */
      fin = fopen(argv[2], "r");
      if (fin == NULL)
      {
        printf("ERROR: couldn't open file with list of peptides, %s\n", argv[2]);
        exit(-1);
      }
      else
      {
        tmp = fscanf(fin, "%d\n", &dtmp);
        numpep = dtmp;
        PEPTIDE *pep = NULL;
        pep = (PEPTIDE *) malloc (numpep * sizeof(PEPTIDE));
        for (k = 0; k < numpep; k++)
        {
          pep[k].sample = (SAMPLEPEP*) malloc (msamples * sizeof(SAMPLEPEP));
          pep[k].numseq = (int*) malloc (MAXLINE * sizeof(int));
        }
        printf("checking %d m/zs:\n", numpep);

        for (k = 0; k < numpep; k++)
        {
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
          for (j = 0; j < msamples; j++)
          {
            pep[k].sample[j].peaks = (double*) malloc ((FITPEAKS + 2) * sizeof(double));
            pep[k].sample[j].betas = (double*) malloc (pep[k].nq * sizeof(double));
            pep[k].sample[j].found = -1;
            strcpy(pep[k].sample[j].name, data[j].name);
          }

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
            free (pep[i].sample[j].betas);
          }
          free (pep[i].sample);
        }
        free (pep);
      }
    }
  }
  return 0;
}

/*****************************************************
 Procedure: getLine
 Description: reads in whole line as a string
 ******************************************************/
int getLine(FILE *fp, int v[], int *zs, int *pg, int *cbm)
{
  int i, ch;
  int z = 0;
  int x = 0;
  int b = 0;

  i = 0;
  ch = fgetc(fp);

  v[i] = ch-65;

  //JW had it like this but this can be flaky for different systems that might use CR LF etc.
  //while (ch != 10)
  while (ch >= 32)
  {
    i++;
    ch = fgetc(fp);
    v[i] = ch-65;
    if (v[i] == 25) z++;
    /* X means a modified cysteine with carboxymethylation (which adds 58 Da) */
    if (v[i] == 23) x++;
    /* B means modified with Gln->pyro-Glu (N-term Q) which subtracts 17 Da */
    if (v[i] == 1) b++;
  }

  (*zs) = z;
  (*pg) = b;
  (*cbm) = x;
  return i;
}


