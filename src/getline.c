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


/*****************************************************
Procedure: getLine
Description: reads in whole line as a string
******************************************************/
int isogetLine(FILE *fp, int v[], int *zs)
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
	    	
	while (ch != 10)
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

