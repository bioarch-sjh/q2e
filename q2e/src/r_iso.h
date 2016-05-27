/*
 * r_iso.h
 *
 *  Created on: 29 Apr 2016
 *      Author: sjh
 */

#ifndef R_ISO_H_
#define R_ISO_H_



float R_iso (PEPTIDE *pep, int numpep, ISODIST *dist, int check);

void R_loadIsotopeTable(/*char *filename,*/ ISOTAB *element, const int num, const int numiso);

void R_iso_seq(char **seq,  double *resultmass, double *resultprob, int *failed);


#endif /* R_ISO_H_ */
