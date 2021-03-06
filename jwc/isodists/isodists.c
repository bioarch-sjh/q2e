/* edited this file to include zeros for S35. */

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

void readIsotopeTable(char *filename, ISOTAB *element, int num, int numiso);
double getIsoDist(int ii, ISOTAB *element, int num, double *peptide, ISODIST *dist);

/****************************************
Procedure:      iso
*****************************************/
void iso (PEPTIDE *pep, int numpep, ISODIST *dist)
{  
  FILE *fp = NULL;
  /* number of elements in isotope table */
  int numelements = 5; 
  /* number of different isotopes in isotope table */
  int numiso = 5; 
  int tmp, dtmp;
  double imass;
  
  int i, j, k;
  
  double *peptide = NULL;
  peptide = (double*) malloc (numelements *  sizeof(double));

  ISOTAB  *element = NULL;
  element = (ISOTAB*) malloc (numelements *  sizeof(ISOTAB));
  for (i = 0; i < numelements; i++)
  {
    element[i].isotope = (double*) malloc (numiso * sizeof(double));
    element[i].abundance = (double*) malloc (numiso * sizeof(double));
  }

  AMINO  *letter = NULL;
  letter = (AMINO*) malloc (26 *  sizeof(AMINO));
  for (i = 0; i < 26; i++)
  {
    letter[i].element = (int*) malloc (numelements * sizeof(int));
  }

  /* read in elements for each amino acid */
  fp = fopen(AMINO_FILE, "r");
  for (i = 0; i < 26; i++)
  {
    for (k = 0; k < numelements; k++)
    {
      tmp = fscanf(fp,"%d ", &dtmp);
      letter[i].element[k] = dtmp;
    }
  }
  fclose(fp);

  readIsotopeTable(ISO_TABLE, element, numelements, numiso);
  
  for (i = 0; i < numpep; i++)
  {
    /*find the peptides that isotope distributions are required for */
    printf("observed mass: %0.1f\n", pep[i].pepmass);
		
    for (k = 0; k < numelements; k++)
    {
      peptide[k] = 0.0;
      for (j = 0; j < pep[i].length; j++)
      {
		peptide[k] = peptide[k] + (float)letter[pep[i].numseq[j]].element[k];
      }
    }

    /* remove 2H10 for each link between amino acids */
    peptide[1] = peptide[1] - 2.0*(float)(pep[i].length - pep[i].zs - 1);
    peptide[3] = peptide[3] - (float)(pep[i].length - pep[i].zs - 1);
    imass = getIsoDist(i, element, numelements, peptide, dist);
      
    printf("calculated %f\n", imass);

    printf("monoisotopic mass            : %0.9f\n", dist[i].prob[0]);
    printf("monoisotopic mass plus one   : %0.9f \n", dist[i].prob[1]);
    printf("monoisotopic mass plus two   : %0.9f \n", dist[i].prob[2]);
    printf("monoisotopic mass plus three : %0.9f \n", dist[i].prob[3]);
    printf("monoisotopic mass plus four  : %0.9f \n", dist[i].prob[4]);
    if (fabs(imass-pep[i].pepmass) > 1.5) printf("SOMETHING IS WRIONG WITH THIS MASS: DIFF = %d\n", (int)(fabs(imass-pep[i].pepmass)+0.5));
  }

  for (i = 0; i < numelements; i++)
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
void readIsotopeTable(char *filename, ISOTAB *element, int num, int numiso)
{
  FILE  *ft = NULL;
  ft = fopen(filename, "r"); 
  
  char name[1024];
  double ftmp, ftmp1;
  int i, j, tmp;
  
  for (i = 0; i < num; i++)
  {
    tmp = fscanf(ft, "%s %lf",  name, &ftmp);
    strcpy(element[i].name, name);
    element[i].avmass = ftmp;
    for (j = 0; j < numiso; j++)
    {
      tmp = fscanf(ft, "%lf %lf",  &ftmp, &ftmp1);
      element[i].isotope[j] = ftmp;
      element[i].abundance[j] = ftmp1/100.0;
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















