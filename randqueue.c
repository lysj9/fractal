#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "type.h"
#include "func.h"

void randqueue2(int n, int *idx, struct vector_s *star, int randnum, int *randidx, double eps)
{
	int i=0,j=0,k=0;
	int itemp=0;
	int randi;
//	int si;
	int accept;
	double eps2=eps*eps;
	double min_sep;
	double *sep;
	sep = (double*) (malloc(randnum*sizeof(double)));

	for (j=0, i=n-1;j<randnum; ){
		if (i<0){
			fprintf(stderr,"not enough stars!\n");
			exit(0);
		}
		randi = (int) ((i+1)*randomz());
		itemp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = itemp;

		accept=1;
		min_sep = DBL_MAX;
		for (k=0;k<j;++k){
			sep[k] = (star[randidx[k]].x - star[randi].x) * (star[randidx[k]].x - star[randi].x) +
				(star[randidx[k]].y - star[randi].y) * (star[randidx[k]].y - star[randi].y) +
				(star[randidx[k]].z - star[randi].z) * (star[randidx[k]].z - star[randi].z);
			if (sep[k]<min_sep) min_sep=sep[k];
			if (min_sep<eps2){
				accept=0;
				break;
			}
		}
		if (accept){
			randidx[j] = idx[i];
			j++;
//			fprintf(stdout,"%e\n",sqrt(min_sep));
		}
		i--;
	}
}

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
