#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "type.h"

#include "randomz.h"

#define SQ(x) ((x) * (x))

void randqueue2(int n, int *idx, int randnum, int *randidx, struct star *star_x, double eps)
{
	int i=0,j=0,k=0;
	int itemp=0;
	int randi;
	int accept;
	double eps2=eps*eps;
	double sep2;
	for (j=0, i=n-1;j<randnum; ) {
		if (i+1 <= randnum-j) {
			fprintf(stderr,"only select %d stars satisfy separation > %le\n",j,eps);
			for (k=j;k<randnum;++k) randidx[k] = idx[k];
			i -= randnum-j;
			break;
		}
		randi = (int) ((i+1)*randomz());
		itemp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = itemp;

		accept=1;
#ifdef _OPENMP
#pragma omp parallel private(sep2) shared(accept)
	{
		int omp_start,omp_stop;
		omp_start = omp_get_thread_num()*j/omp_get_num_threads();
		omp_stop = (omp_get_thread_num()+1)*j/omp_get_num_threads();

		while (omp_start < omp_stop && accept) {
			sep2 = SQ(star_x[randidx[omp_start]].x[0] - star_x[randi].x[0]) +
				   SQ(star_x[randidx[omp_start]].x[1] - star_x[randi].x[1]) +
				   SQ(star_x[randidx[omp_start]].x[2] - star_x[randi].x[2]);
			if (sep2<eps2) {
#pragma omp critical
				accept=0;
			}
			omp_start++;
		}
	}
#else
		for (k=0;k<j;++k) {
			sep2 = SQ(star_x[randidx[k]].x[0] - star_x[randi].x[0]) +
				   SQ(star_x[randidx[k]].x[1] - star_x[randi].x[1]) +
				   SQ(star_x[randidx[k]].x[2] - star_x[randi].x[2]);
			if (sep2<eps2) {
				accept=0;
				break;
			}
		}
#endif

		if (accept) {
			randidx[j] = idx[i];
			j++;
		}
		i--;
	}
	fprintf(stderr,"n = %d, randnum = %d, %d stars tested...\n",n,randnum,n-i-1);
}

void randqueue(int n, int *idx, int randnum, int *randidx)
{
	if (randnum>n) {
		fprintf(stderr,"random queue length longer than \
			original array length\n");
		exit(-1);
	}
	int i=0,temp=0,j=0;
	i = n-1;
	int randi;
	for (j=0;j<randnum;++j) {
		randi = (int) ( (i+1)*randomz() );
		temp = idx[randi];
		idx[randi] = idx[i];
		idx[i] = temp;
		randidx[j] = idx[i];
		i--;
	}
}
