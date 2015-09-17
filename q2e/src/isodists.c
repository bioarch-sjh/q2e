#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "constants.h"
#include "files.h"

#define DEBUG_R

#define NUMELEMENTS 5
#define NUMLETTERS 26
#define NUMISOTOPES 5

void readIsotopeTable(char *filename, ISOTAB *element, int num, int numiso);
double getIsoDist(int ii, ISOTAB *element, int num, double *peptide, ISODIST *dist);

/****************************************
Procedure:      iso
*****************************************/
void iso (PEPTIDE *pep, int numpep, ISODIST *dist)
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
  printf("x = %d\n",x);

  const int AMINODATA[NUMLETTERS][NUMELEMENTS] = {
     {3,  7, 1,  2, 0},
	 {0, -1, 0, -1, 0},
	 {3,  7, 1,  2, 1},
	 {4,  7, 1,  4, 0},
	 {5,  9, 1,  4, 0},
	 {9, 11, 1,  2, 0},
	 {2,  5, 1,  2, 0},
	 {6,  9, 3,  2, 0},
	 {6, 13, 1,  2, 0},
	 {0,  0, 0,  0, 0},
	 {6, 14, 2,  2, 0},
	 {6, 13, 1,  2, 0},
	 {5, 11, 1,  2, 1},
	 {4,  8, 2,  3, 0},
	 {0,  0, 0,  0, 0},
	 {5,  9, 1,  2, 0},
	 {5, 10, 2,  3, 0},
	 {6, 14, 4,  2, 0},
	 {3,  7, 1,  3, 0},
	 {4,  9, 1,  3, 0},
	 {0,  0, 0,  0, 0},
	 {5, 11, 1,  2, 0},
	{11, 12, 2,  2, 0},
	 {2,  2, 0,  2, 0},
	 {9, 11, 1,  3, 0},
	 {0,  0, 0,  1, 0}
	};


  /* number of different isotopes in isotope table */
  int numiso = 5;

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
  letter = (AMINO*) malloc (26 *  sizeof(AMINO));
  for (i = 0; i < 26; i++)
  {
    letter[i].element = (int*) malloc (NUMELEMENTS * sizeof(int));
  }

  //TODO: Although we've replaced the file aminomasses with the AMINODATA array, we need a better way to load other amino masses if we need to.
  /* read in elements for each amino acid */
  //fp = fopen(AMINO_FILE, "r");
  for (i = 0; i < 26; i++)
  {
    for (k = 0; k < NUMELEMENTS; k++)
    {
      //tmp = fscanf(fp,"%d ", &dtmp);
      //letter[i].element[k] = dtmp;
      letter[i].element[k] = AMINODATA[i][k];
    }
  }
  //fclose(fp);

  readIsotopeTable(ISO_TABLE, element, NUMELEMENTS, numiso);

  printf("success reading iso table\n");


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
      if (fabs(imass-pep[i].pepmass) > 1.5) printf("ERROR: SOMETHING IS WRONG HERE, THE DIFERENCE BETWEEN THE GIVEN AND CALCULATED MASSES IS %d\n", (int)(fabs(imass-pep[i].pepmass)+0.5));
    }
  }

  for (i = 0; i < NUMELEMENTS; i++)
  {
    free(element[i].isotope);
    free(element[i].abundance);
    free(letter[i].element);
  }
  free(element);
  free(letter);
  free(peptide);

}

/*************************************************************************************************************
Procedure: readIsotopeTable
inputs the isotope info
************************************************************************************************************/
void readIsotopeTable(char *filename, ISOTAB *element, const int num, const int numiso)
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

/*************************************************************************************************************
Procedure: gammln
From Numerical Recipes: use to calculate binomial coefficients
************************************************************************************************************/
float gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return (-tmp+log(2.5066282746310005*ser/x));
}


/*************************************************************************************************************
Procedure: bico
From Numerical Recipes: calculates binomial coefficients
************************************************************************************************************/
float bico(int n, int k)
{
  float lnfactn, lnfactk, lnfactnk, bin;

  if (k > n) bin = 0.0;
  else if (k < 0) bin = 0.0;
  else if (k == 0) bin = 1.0;
  else
  {  
    lnfactn = gammln((float)(n+1));
    lnfactk = gammln((float)(k+1));
    lnfactnk = gammln((float)(n-k+1));
    bin = floor(0.5+exp(lnfactn - lnfactk - lnfactnk));
  }
  
  return (bin);
}

/*************************************************************************************************************
Procedure: getIsoDist

************************************************************************************************************/
double getIsoDist(int ii, ISOTAB *element, int num, double *peptide, ISODIST *dist)
{
  double mult, power, probj, probk, probl, probm;
  double probjk, probkl, probjkl, probklm, probjklm;
  int i, j, k, l, m;
  
  /* get monoisotopic mass and it's probabilty */
  dist[ii].mass[0] = 0.0;
  dist[ii].prob[0] = 100.0;

  for (i = 0; i < num; i++)
  {
    power = peptide[i];
    dist[ii].mass[0] = dist[ii].mass[0] + peptide[i]*element[i].isotope[0];
    dist[ii].prob[0] = dist[ii].prob[0] * pow(element[i].abundance[0], power);
  }
  
  /* get probabilty of monoisotopic mass plus 1 */
  dist[ii].prob[1] = 0.0;
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[1] =   dist[ii].prob[1] + probj;
  }
  dist[ii].prob[1] =  dist[ii].prob[1]*100.0; 

  /* get probabilty of monoisotopic mass plus 2 */
  dist[ii].prob[2] = 0.0;
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[2] = dist[ii].prob[2] + probj;
  }  
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = bico(peptide[j], 2);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[1];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[2] = dist[ii].prob[2] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[2] = dist[ii].prob[2] + probjk;
    }
  }
  dist[ii].prob[2] =  dist[ii].prob[2]*100.0; 
 
  /* get probabilty of monoisotopic mass plus 3 */
  dist[ii].prob[3] = 0.0;
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[3];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[3] = dist[ii].prob[3] + probj;
  } 
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = peptide[j]*(peptide[j]-1);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[2];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[3] = dist[ii].prob[3] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 3.0;
    mult = bico(peptide[j], 3);
    probj = mult * pow(element[j].abundance[0], power)*pow(element[j].abundance[1],3);
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[3] = dist[ii].prob[3] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[3] = dist[ii].prob[3] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = bico(peptide[j], 2);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[3] = dist[ii].prob[3] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (l = 0; l < num; l++)
      {
        probl = 0.0;
        if (l > k)
        {
          power = peptide[l] - 1.0;
          mult = peptide[l];
          probl = mult * pow(element[l].abundance[0], power)*element[l].abundance[1];
        }
        probkl = probk*probl;
        for (i = 0; i < num; i++)
        {
          if ((i != j) && (i != k) && (i != l))
          {
            power = peptide[i];
            probkl = probkl * pow(element[i].abundance[0], power);
          }
        }
        probjkl = probj*probkl;
        dist[ii].prob[3] = dist[ii].prob[3] + probjkl;
      }
    }
  }
  dist[ii].prob[3] =  dist[ii].prob[3]*100.0; 
 
  /* get probabilty of monoisotopic mass plus 4 */
  dist[ii].prob[4] = 0.0;
  /* note there is no abundance for iso 4 or more */
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = peptide[j]*(peptide[j]-1);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[3];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[4] = dist[ii].prob[4] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = bico(peptide[j], 2);
    probj = mult * pow(element[j].abundance[0], power)*pow(element[j].abundance[2],2);
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[4] = dist[ii].prob[4] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 3.0;
    mult = peptide[j]*bico(peptide[j]-1, 2);
    probj = mult * pow(element[j].abundance[0], power)*pow(element[j].abundance[1],2)*element[j].abundance[2];
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[4] = dist[ii].prob[4] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 4;
    mult = bico(peptide[j], 4);
    probj = mult * pow(element[j].abundance[0], power)*pow(element[j].abundance[1],4);
    for (i = 0; i < num; i++)
    {
      if (i != j)
      {
        power = peptide[i];
        probj = probj * pow(element[i].abundance[0], power);
      }
    }
    dist[ii].prob[4] = dist[ii].prob[4] + probj;
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[3];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[2];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = peptide[j]*(peptide[j]-1);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2]*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (l = 0; l < num; l++)
      {
        probl = 0.0;
        if (l > k)
        {
          power = peptide[l] - 1.0;
          mult = peptide[l];
          probl = mult * pow(element[l].abundance[0], power)*element[l].abundance[1];
        }
        probkl = probk*probl;
        for (i = 0; i < num; i++)
        {
          if ((i != j) && (i != k) && (i != l))
          {
            power = peptide[i];
            probkl = probkl * pow(element[i].abundance[0], power);
          }
        }
        probjkl = probj*probkl;
        dist[ii].prob[4] = dist[ii].prob[4] + probjkl;
      }
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[2];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 2.0;
        mult = bico(peptide[k], 2);
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1]*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 3.0;
    mult = bico(peptide[j],3);
    probj = mult * pow(element[j].abundance[0], power)*pow(element[j].abundance[1],3);
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = bico(peptide[j], 2);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 2.0;
        mult = bico(peptide[k], 2);
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1]*element[k].abundance[1];
      }
      for (i = 0; i < num; i++)
      {
        if ((i != j) && (i != k))
        {
          power = peptide[i];
          probk = probk * pow(element[i].abundance[0], power);
        }
      }
      probjk = probj*probk;
      dist[ii].prob[4] = dist[ii].prob[4] + probjk;
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 2.0;
    mult = bico(peptide[j], 2);
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1]*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k != j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (l = 0; l < num; l++)
      {
        probl = 0.0;
        if ((l != j) && (l > k))
        {
          power = peptide[l] - 1.0;
          mult = peptide[l];
          probl = mult * pow(element[l].abundance[0], power)*element[l].abundance[1];
        }
        probkl = probk*probl;
        for (i = 0; i < num; i++)
        {
          if ((i != j) && (i != k) && (i != l))
          {
            power = peptide[i];
            probkl = probkl * pow(element[i].abundance[0], power);
          }
        }
        probjkl = probj*probkl;
        dist[ii].prob[4] = dist[ii].prob[4] + probjkl;
      }
    }
  }
  for (j = 0; j < num; j++)
  {
    power = peptide[j] - 1.0;
    mult = peptide[j];
    probj = mult * pow(element[j].abundance[0], power)*element[j].abundance[1];
    for (k = 0; k < num; k++)
    {
      probk = 0.0;
      if (k > j)
      {
        power = peptide[k] - 1.0;
        mult = peptide[k];
        probk = mult * pow(element[k].abundance[0], power)*element[k].abundance[1];
      }
      for (l = 0; l < num; l++)
      {
        probl = 0.0;
        if (l > k)
        {
          power = peptide[l] - 1.0;
          mult = peptide[l];
          probl = mult * pow(element[l].abundance[0], power)*element[l].abundance[1];
        }
        probkl = probk*probl;
        for (m = 0; m < num; m++)
        {
          probm = 0.0;
          if (m > l)
          {
            power = peptide[m] - 1.0;
            mult = peptide[m];
            probm = mult * pow(element[m].abundance[0], power)*element[m].abundance[1];
          }
          probklm = probkl*probm;
          for (i = 0; i < num; i++)
          {
            if ((i != j) && (i != k) && (i != l) && (i != m))
            {
              power = peptide[i];
              probklm = probklm * pow(element[i].abundance[0], power);
            }
          }
          probjklm = probj*probklm;
          dist[ii].prob[4] = dist[ii].prob[4] + probjklm;
        }
      }
    }
  }
  dist[ii].prob[4] =  dist[ii].prob[4]*100.0; 
 
 double max = 0.0;
  for (i = 0; i < 5; i++)
  {
    if (dist[ii].prob[i] > max) max =   dist[ii].prob[i]; 
  }
  for (i = 0; i < 5; i++)
  {
    dist[ii].prob[i] = dist[ii].prob[i]/max; 
  }

  return(dist[ii].mass[0]);  
}














