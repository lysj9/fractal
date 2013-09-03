#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "type.h"

#include "quick_sort.h"

double nbody_scale(int N_node, double q, struct vector_s *star, int *nnbmax_out, double *rs0_out)
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
//	double rvir_nbody = 0.5/(1-q); // virial radius in nbody-unit
//	double vvir_nbody = 0.5*q/(1-q); // mean square velocity in nbody-unit
	double r2,rscale,vscale;
	double *rij;
//	rij = (double*) malloc (N_node*sizeof(double));

	struct vector_s centre={0,0,0,0,0,0,0};
	// centre of mass correction
	for (i=0;i<N_node;++i){
		centre.m  += star[i].m;
		centre.x  += star[i].m * star[i].x;
		centre.y  += star[i].m * star[i].y;
		centre.z  += star[i].m * star[i].z;
		centre.vx += star[i].m * star[i].vx;
		centre.vy += star[i].m * star[i].vy;
		centre.vz += star[i].m * star[i].vz;
	}
	centre.x  /= centre.m;
	centre.y  /= centre.m;
	centre.z  /= centre.m;
	centre.vx /= centre.m;
	centre.vy /= centre.m;
	centre.vz /= centre.m;

	// all stars including binaries and com-binaries
	for (i=0;i<N_node;++i){
		star[i].m  /= centre.m;
		star[i].x  -= centre.x;
		star[i].y  -= centre.y;
		star[i].z  -= centre.z;
		star[i].vx -= centre.vx;
		star[i].vy -= centre.vy;
		star[i].vz -= centre.vz;
	}

	RS0=0;
//	for (i=0;i<N_node;++i) rij[i] = DBL_MAX;
//#pragma omp parallel for private(i,j,r2,rij) reduction(-:Ep) reduction(+:RS0) reduction(+:Ek)
#pragma omp parallel private(rij)
	{
	rij = (double*) malloc (N_node*sizeof(double));
#pragma omp for private(i,j,r2) reduction(-:Ep) reduction(+:RS0) reduction(+:Ek)
	for (i=0;i<N_node;++i){
//		rij[i] = DBL_MAX;
		for (j=0;j<N_node;++j){
			if (i!=j) {
				r2 = ( star[i].x - star[j].x )*( star[i].x - star[j].x ) +\
					 ( star[i].y - star[j].y )*( star[i].y - star[j].y ) +\
					 ( star[i].z - star[j].z )*( star[i].z - star[j].z );
				rij[j] = r2;
			}
			if (j>i) Ep -= star[i].m*star[j].m / sqrt(r2);
		}
		RS0 += sqrt( quick_select(rij,NNBMAX/5,N_node) );
//		quick_sort1(rij,N_node);
//		RS0 += sqrt( rij[NNBMAX/5] );

//		v2 = star[i].vx*star[i].vx + star[i].vy*star[i].vy + star[i].vz*star[i].vz;
		Ek += star[i].m * ( star[i].vx*star[i].vx + star[i].vy*star[i].vy + star[i].vz*star[i].vz );
	}
	}
	Ek *= 0.5;
	RS0 /= N_node;
//	*rs0_out = RS0/N_node;

/*
	for (i=0;i<N_node;++i){
		for (j=i+1;j<N_node;++j){
			r2 = ( star[i].x - star[j].x )*( star[i].x - star[j].x ) +\
				 ( star[i].y - star[j].y )*( star[i].y - star[j].y ) +\
				 ( star[i].z - star[j].z )*( star[i].z - star[j].z );
			Ep -= star[i].m*star[j].m / sqrt(r2);
		}
//		v2 = star[i].vx*star[i].vx + star[i].vy*star[i].vy + star[i].vz*star[i].vz;
		Ek += star[i].m * ( star[i].vx*star[i].vx + star[i].vy*star[i].vy + star[i].vz*star[i].vz );
	}
	Ek *= 0.5;
*/

	rscale = Ep/Ep_nbody;
	vscale = sqrt(Ek_nbody/Ek);
	fprintf(stderr,"scale: rscale, vscale: %lf, %lf\n",rscale,vscale);

	*rs0_out = RS0*rscale;
	fprintf(stderr,"NNBMAX=%d, RS0=%lf\n",NNBMAX,*rs0_out);
	for (i=0;i<N_node;++i){
		star[i].x  *= rscale;
		star[i].y  *= rscale;
		star[i].z  *= rscale;
		star[i].vx *= vscale;
		star[i].vy *= vscale;
		star[i].vz *= vscale;
	}

	return centre.m;
}

/*
void binary_scale()
{
}

void binary_energy(int nbin, struct vector_s *star)
{
	int i,j;
	double Ek,Ep;
}
*/

void sort_radius(int N_cm, struct vector_s *star, double *r2_sort, int *idx)
{
	int i;
	for (i=0;i<N_cm;++i) {
		r2_sort[i] = (star[i].x*star[i].x) + 
					 (star[i].y*star[i].y) + 
					 (star[i].z*star[i].z);
		idx[i] = i;
	}
	quick_sort2(r2_sort,idx,N_cm);
	return;
}

double get_radius(int N_cm, struct vector_s *star, double truncate, double *r2_sort, int *idx)
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
		m[i] = star[i].m;
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
