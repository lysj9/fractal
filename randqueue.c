#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "type.h"
#include "func.h"

void randqueue2(int n, int *idx, int randnum, int *randidx, struct vector_s *star, double eps)
{
	int i=0,j=0,k=0;
	int itemp=0;
	int randi;
	int accept;
	double eps2=eps*eps;
	//double min_sep;
	double *sep;
	sep = (double*) (malloc(randnum*sizeof(double)));
	for (j=0, i=n-1;j<randnum; ) {
		if (i+1 <= randnum-j) {
			fprintf(stderr,"only select %d stars satisfy separation > %le\n",j,eps);
			for (k=j;k<randnum;++k) randidx[k] = idx[k];
			i -= randnum-j-1;
			break;
		}
		randi = (int) ((i+1)*randomz());
		itemp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = itemp;

		accept=1;
		//min_sep = DBL_MAX;
//#pragma omp parallel for private(sep)
		for (k=0;k<j;++k) {
			/*
			double dx,dy,dz;
			double dxy;
			dx = star[randidx[k]].x - star[randi].x;
			if (dx>eps) continue;
			dy = star[randidx[k]].y - star[randi].y;
			if (dy>eps) continue;
			dxy = dx*dx + dy*dy;
			if (dxy>eps2) continue;
			dz = star[randidx[k]].z - star[randi].z;
			if (dz>eps) continue;
			sep[k] = dxy + dz*dz;
//			sep[k] = dx*dx + dy*dy + dz*dz;
			if (sep[k] < eps2) {
				accept = 0;
				break;
			}
			*/
			sep[k] = (star[randidx[k]].x - star[randi].x) * (star[randidx[k]].x - star[randi].x) +
				(star[randidx[k]].y - star[randi].y) * (star[randidx[k]].y - star[randi].y) +
				(star[randidx[k]].z - star[randi].z) * (star[randidx[k]].z - star[randi].z);
			if (sep[k]<eps2) {
				accept=0;
				break;
			}
		}
		if (accept) {
			randidx[j] = idx[i];
			j++;
//			fprintf(stdout,"%e\n",sqrt(min_sep));
		}
		i--;
	}
	free(sep);
	fprintf(stderr,"n = %d, randnum = %d, %d stars tested...\n",n,randnum,n-i);
	exit(0);
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
