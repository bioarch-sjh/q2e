/* edited this file to include zeros for S35. */
#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "constants.h"

void readIsotopeTable(char *filename, ISOTAB *element, int num, int numiso);
double getIsoDist(int ii, ISOTAB *element, int num, double *peptide, ISODIST *dist);

/****************************************
Procedure:      iso
*****************************************/
void iso (PEPTIDE *pep, int numpep, ISODIST *dist)
{
  /* number of elements in isotope table */
  const int numelements = 5;
  /* number of different isotopes in isotope table */
  const int numiso = 5;
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

  /* elements for each amino acid */
  for (i = 0; i < 26; i++)
  {
    for (k = 0; k < numelements; k++)
    {
	  letter[i].element[k] = 0;
    }
  }

	/* A */
	letter[0].element[0] = 3;
	letter[0].element[1] = 7;
	letter[0].element[2] = 1;
	letter[0].element[3] = 2;
	/* B */
	letter[1].element[1] = -1;
	letter[1].element[3] = -1;
	/* C */
	letter[2].element[0] = 3;
	letter[2].element[1] = 7;
	letter[2].element[2] = 1;
	letter[2].element[3] = 2;
	letter[2].element[4] = 1;
	/* D */
	letter[3].element[0] = 4;
	letter[3].element[1] = 7;
	letter[3].element[2] = 1;
	letter[3].element[3] = 4;
	/* E */
	letter[4].element[0] = 5;
	letter[4].element[1] = 9;
	letter[4].element[2] = 1;
	letter[4].element[3] = 4;
	/* F */
	letter[5].element[0] = 9;
	letter[5].element[1] = 11;
	letter[5].element[2] = 1;
	letter[5].element[3] = 2;
	/* G */
	letter[6].element[0] = 2;
	letter[6].element[1] = 5;
	letter[6].element[2] = 1;
	letter[6].element[3] = 2;
	/* H */
	letter[7].element[0] = 6;
	letter[7].element[1] = 9;
	letter[7].element[2] = 3;
	letter[7].element[3] = 2;
	/* I */
	letter[8].element[0] = 6;
	letter[8].element[1] = 13;
	letter[8].element[2] = 1;
	letter[8].element[3] = 2;
	/* K */
	letter[10].element[0] = 6;
	letter[10].element[1] = 14;
	letter[10].element[2] = 2;
	letter[10].element[3] = 2;
	/* L */
	letter[11].element[0] = 6;
	letter[11].element[1] = 13;
	letter[11].element[2] = 1;
	letter[11].element[3] = 2;
	/* M */
	letter[12].element[0] = 5;
	letter[12].element[1] = 11;
	letter[12].element[2] = 1;
	letter[12].element[3] = 2;
	letter[12].element[4] = 1;
	/* N */
	letter[13].element[0] = 4;
	letter[13].element[1] = 8;
	letter[13].element[2] = 2;
	letter[13].element[3] = 3;
	/* P */
	letter[15].element[0] = 5;
	letter[15].element[1] = 9;
	letter[15].element[2] = 1;
	letter[15].element[3] = 2;
	/* Q */
	letter[16].element[0] = 5;
	letter[16].element[1] = 10;
	letter[16].element[2] = 2;
	letter[16].element[3] = 3;
	/* R */
	letter[17].element[0] = 6;
	letter[17].element[1] = 14;
	letter[17].element[2] = 4;
	letter[17].element[3] = 2;
	/* S */
	letter[18].element[0] = 3;
	letter[18].element[1] = 7;
	letter[18].element[2] = 1;
	letter[18].element[3] = 3;
	/* T */
	letter[19].element[0] = 4;
	letter[19].element[1] = 9;
	letter[19].element[2] = 1;
	letter[19].element[3] = 3;
	/* V */
	letter[21].element[0] = 5;
	letter[21].element[1] = 11;
	letter[21].element[2] = 1;
	letter[21].element[3] = 2;
	/* W */
	letter[22].element[0] = 11;
	letter[22].element[1] = 12;
	letter[22].element[2] = 2;
	letter[22].element[3] = 2;
	/* X THIS IS FOR CAROLINES MODIFICATION */
	letter[23].element[0] = 2;
	letter[23].element[1] = 2;
	letter[23].element[2] = 0;
	letter[23].element[3] = 2;
	/* Y */
	letter[24].element[0] = 9;
	letter[24].element[1] = 11;
	letter[24].element[2] = 1;
	letter[24].element[3] = 3;
	/* Z FOR ADDING HYDROXYLATIONS */
	letter[25].element[3] = 1;


	/* 	isotope info */
	for (i = 0; i < numelements; i++)
	{
	  for (j = 0; j < numiso; j++)
	  {
		element[i].isotope[j] = 0.0;
		element[i].abundance[j] = 0.0;
	  }
	}

    /* C */
	element[0].avmass = 12.0107;
	element[0].isotope[0] = 12.0;
	element[0].abundance[0] = 0.988922;
	element[0].isotope[1] = 13.00335;
	element[0].abundance[1] = 0.0110780;

    /* H */
	element[1].avmass = 1.00794;
	element[1].isotope[0] = 1.00783;
	element[1].abundance[0] = 0.9998443;
	element[1].isotope[1] = 2.0141;
	element[1].abundance[1] = 0.0001557;

    /* N */
	element[2].avmass = 14.0067;
	element[2].isotope[0] = 14.00307;
	element[2].abundance[0] = 0.996337;
	element[2].isotope[1] = 15.00011;
	element[2].abundance[1] = 0.0036630;

    /* O */
	element[3].avmass = 15.9994;
	element[3].isotope[0] = 15.99491;
	element[3].abundance[0] = 0.997628;
	element[3].isotope[1] = 16.99913;
	element[3].abundance[1] = 0.000372;
	element[3].isotope[2] = 17.99916;
	element[3].abundance[2] = 0.0020004;

    /* S */
	element[4].avmass = 32.065;
	element[4].isotope[0] = 31.97207;
	element[4].abundance[0] = 0.95018;
	element[4].isotope[1] = 32.97146;
	element[4].abundance[1] = 0.0075;
	element[4].isotope[2] = 33.96787;
	element[4].abundance[2] = 0.04215;
	element[4].isotope[4] = 35.96708;
	element[4].abundance[4] = 0.00017;

  for (i = 0; i < numpep; i++)
  {
    /*find the peptides that isotope distributions are required for */
    if (pep[i].found == 1)
    {
#ifdef DBG_PRINT
      Rprintf("%0.1f\n", pep[i].pepmass);
#endif
      for (k = 0; k < numelements; k++)
      {
        peptide[k] = 0.0;
        for (j = 0; j < pep[i].length; j++)
        {
           peptide[k] = peptide[k] + (float)letter[pep[i].numseq[j]].element[k];
        }
      }


#ifdef DBG_PRINT
      Rprintf("peptide[1] = %0.2f; peptide[3] = %0.2f\n",peptide[1],peptide[3]);
#endif

      /* remove 2H10 for each link between amino acids */
      peptide[1] = peptide[1] - 2.0*(float)(pep[i].length - pep[i].zs - pep[i].pg - pep[i].cbm - 1);
      peptide[3] = peptide[3] - (float)(pep[i].length - pep[i].zs - pep[i].pg - pep[i].cbm - 1);

#ifdef DBG_PRINT
      Rprintf("peptide[1] = %0.2f; peptide[3] = %0.2f\n",peptide[1],peptide[3]);
#endif

      imass = getIsoDist(i, element, numelements, peptide, dist);

#ifdef DBG_PRINT
      Rprintf("calculated %f\n", imass);
      Rprintf("monoisotopic mass            : %0.9f\n", dist[i].prob[0]);
      Rprintf("monoisotopic mass plus one   : %0.9f \n", dist[i].prob[1]);
      Rprintf("monoisotopic mass plus two   : %0.9f \n", dist[i].prob[2]);
      Rprintf("monoisotopic mass plus three : %0.9f \n", dist[i].prob[3]);
      Rprintf("monoisotopic mass plus four  : %0.9f \n", dist[i].prob[4]);
#endif
      if (fabs(imass-pep[i].pepmass) > 1.5)
	  {
		Rprintf("ERROR: SOMETHING IS WRONG HERE, THE DIFERENCE BETWEEN THE GIVEN (%f) AND CALCULATED (%f) MASSES IS %f\n",(float)pep[i].pepmass, (float)imass, (float)(fabs(imass-pep[i].pepmass)+0.5));
		//exit(-1);
		stop("IMASS calculation doesn't match value in peptideList");
	  }
    }
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


#ifdef DBG_PRINT
  Rprintf("getIsoDist Loop 1:\n");
#endif
  for (i = 0; i < num; i++){
    power = peptide[i];
    dist[ii].mass[0] = dist[ii].mass[0] + peptide[i]*element[i].isotope[0];
    dist[ii].prob[0] = dist[ii].prob[0] * pow(element[i].abundance[0], power);

#ifdef DBG_PRINT
    Rprintf("peptide[%d] = %0.2f, mass = %0.2f, prob = %0.2f\n",i,peptide[i], dist[ii].mass[0],dist[ii].prob[0]);
#endif

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















