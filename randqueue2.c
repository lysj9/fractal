#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "type.h"
#include "func.h"

void randqueue2(int n, int randnum, struct vector_s *star, double eps)
{
	int i=0,j=0,k=0;
	int itemp=0;
	i=n-1;
	int randi;
	int si;
	double min_sep;
	double *sep;
	sep = (double*) (malloc(randnum*sizeof(double)));
	for (j=0;j<randnum;++j){
		randi = (int) ((i+1)*randomz());
		min_sep = DBL_MAX;
		for (k=0;k<j;++k){
			sep[k] = (star[randidx[k]].x - star[randi].x) * (star[randidx[k]].x - star[randi].x) +
				(star[randidx[k]].y - star[randi].y) * (star[randidx[k]].y - star[randi].y) +
				(star[randidx[k]].z - star[randi].z) * (star[randidx[k]].z - star[randi].z);
			if (sep[k]<min_sep) min_sep=sep[k];
			if (min_sep<eps) break;
		}
		itemp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = itemp;
		randidx[j] = idx[i];
		--i;

void randqueue( int        n,
				int     *idx,
				int  randnum,
				int *randidx )
{
	if ( randnum>n ){
		fprintf(stderr,"random queue length longer than \
			original array length\n");
		exit(-1);
	}
	int i=0,temp=0,j=0;
	i = n-1;
	int randi;
//	randomz_seed(1);
	for ( j=0;j<randnum;++j ){
//	for ( j=0;j<n;++j ){
		randi = (int) ( (i+1)*randomz() );
		temp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = temp;
		randidx[j] = idx[i];
		i--;
	}
}
