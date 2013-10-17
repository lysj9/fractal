#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "type.h"

#include "quick_select.h"
#include "quick_sort.h"

#define SQ(x) ((x) * (x))

double nbody_scale(int N_node, double q, struct star *s, int *nnbmax_out, double *rs0_out)
{
	int NNBMAX;
	NNBMAX = 1.25*sqrt(N_node);
	if (NNBMAX<30) NNBMAX=30;
	if (N_node<=NNBMAX) NNBMAX=0.5*N_node;
	*nnbmax_out = NNBMAX;
	double RS0;
	int i,j;
	double Ek=0,Ep=0;
	// scale to totle energy = -0.25
	double Ep_nbody = -0.25/(1-q);
	double Ek_nbody = 0.25*q/(1-q);
	double r2,rscale,vscale;
	double *rij;

	struct star centre={0,{0,0,0,0,0,0}};
	// centre of mass correction
	for (i=0;i<N_node;++i) {
		centre.m += s[i].m;
		centre.x[0] += s[i].m * s[i].x[0];
		centre.x[1] += s[i].m * s[i].x[1];
		centre.x[2] += s[i].m * s[i].x[2];
		centre.x[3] += s[i].m * s[i].x[3];
		centre.x[4] += s[i].m * s[i].x[4];
		centre.x[5] += s[i].m * s[i].x[5];
	}
	for (i=0;i<6;++i) centre.x[i] /= centre.m;

	// all stars including binaries and com-binaries
	for (i=0;i<N_node;++i) {
		s[i].m /= centre.m;
		s[i].x[0] -= centre.x[0];
		s[i].x[1] -= centre.x[1];
		s[i].x[2] -= centre.x[2];
		s[i].x[3] -= centre.x[3];
		s[i].x[4] -= centre.x[4];
		s[i].x[5] -= centre.x[5];
	}

	RS0=0;

#ifdef _OPENMP
#pragma omp parallel private(rij)
	{
	rij = (double*) malloc (N_node*sizeof(double));
#pragma omp for private(i,j,r2) reduction(-:Ep) reduction(+:RS0) reduction(+:Ek)
	for (i=0;i<N_node;++i) {
		rij[i] = DBL_MAX;
		for (j=0;j<N_node;++j) {
			if (i!=j) {
				r2 = (s[i].x[0] - s[j].x[0])*(s[i].x[0] - s[j].x[0]) +
					 (s[i].x[1] - s[j].x[1])*(s[i].x[1] - s[j].x[1]) +
					 (s[i].x[2] - s[j].x[2])*(s[i].x[2] - s[j].x[2]);
				rij[j] = r2;
			}
			if (j>i) Ep -= s[i].m*s[j].m / sqrt(r2);
		}
		RS0 += sqrt( quick_select(rij,NNBMAX/5,N_node) );
		Ek  += s[i].m * (s[i].x[3]*s[i].x[3] + s[i].x[4]*s[i].x[4] + s[i].x[5]*s[i].x[5]);
	}
	free(rij);
	}
	Ek *= 0.5;
	RS0 /= N_node;

#else

	rij = (double*) malloc (N_node*sizeof(double));
	for (i=0;i<N_node;++i) {
		rij[i] = DBL_MAX;
		for (j=0;j<N_node;++j) {
			if (i != j) {
				r2 = (s[i].x[0] - s[j].x[0])*(s[i].x[0] - s[j].x[0]) +
					 (s[i].x[1] - s[j].x[1])*(s[i].x[1] - s[j].x[1]) +
					 (s[i].x[2] - s[j].x[2])*(s[i].x[2] - s[j].x[2]);
				rij[j] = r2;
			}
			if (j>i) Ep -= s[i].m*s[j].m / sqrt(r2);
		}
		RS0 += sqrt( quick_select(rij,NNBMAX/5,N_node) );
		Ek  += s[i].m * (s[i].x[3]*s[i].x[3] + s[i].x[4]*s[i].x[4] + s[i].x[5]*s[i].x[5]);
	}
	free(rij);

	Ek *= 0.5;
	RS0 /= N_node;

#endif

	rscale = Ep/Ep_nbody;
	if (Ek == 0) {
		vscale = 1;
		q = 0;
		fprintf(stderr,"scale: kinetic energy = 0! no scale for velocity!\n");
	} else {
		vscale = sqrt(Ek_nbody/Ek);
	}
	fprintf(stderr,"scale: rscale, vscale: %lf, %lf\n",rscale,vscale);

	*rs0_out = RS0*rscale;
	fprintf(stderr,"NNBMAX=%d, RS0=%lf\n",NNBMAX,*rs0_out);
	for (i=0;i<N_node;++i) {
		s[i].x[0] *= rscale;
		s[i].x[1] *= rscale;
		s[i].x[2] *= rscale;
		s[i].x[3] *= vscale;
		s[i].x[4] *= vscale;
		s[i].x[5] *= vscale;
	}

	return centre.m;
}

void sort_radius(int N_cm, struct star *s, double *r2_sort, int *idx)
{
	int i;
	for (i=0;i<N_cm;++i) {
		r2_sort[i] = (s[i].x[0]*s[i].x[0]) + 
					 (s[i].x[1]*s[i].x[1]) + 
					 (s[i].x[2]*s[i].x[2]);
		idx[i] = i;
	}
	quick_sort2(r2_sort,idx,N_cm);
	return;
}

double get_radius(int N_cm, struct star *s, double truncate, double *r2_sort, int *idx)
{
	int i;
	double r_tr;
	double m_tr=0;
	double *m;
	double rl,rh,ml,mh;
	if ( NULL == (m = (double*) malloc(N_cm*sizeof(double))) ) {
		fprintf(stderr,"malloc m failed...\n");
		exit(0);
	}
	for (i=0;i<N_cm;++i) {
		m[i] = s[i].m;
	}
	m_tr = m[idx[0]];
	if (m_tr>truncate) {
		r_tr = truncate*sqrt(r2_sort[0]);
	} else {
		for (i=1;i<N_cm;++i) {
			m_tr += m[idx[i]];
			if (m_tr>truncate) break;
		}
		rl = sqrt(r2_sort[i-1]);
		rh = sqrt(r2_sort[i]);
		mh = m_tr;
		ml = m_tr - m[idx[i-1]];
		r_tr = rl + (truncate-ml)/(mh-ml)*(rh-rl);
	}
	free(m);
	return r_tr;
}
