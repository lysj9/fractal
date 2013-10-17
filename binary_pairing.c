/* Kouwenhoven et al 2009 A&A: Exploring the consequences of pairing algorithms for binary stars */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "type.h"

#include "generate_mass.h"
#include "make_mass.h"
#include "randomz.h"

#define SWAP(a,b,type) do {type temp; temp=a; a=b; b=temp;} while(0)

static
double hq(double q0, double q1)
{
	return q0 + (q1-q0)*randomz();
}

static
void RP(struct star *x, int n, int *nbin, double ml, double mh)
{
	int i,j;
	int nb,ns;
	double *m2;
	nb = *nbin;
	ns = n - nb;
	m2 = (double*) malloc(nb*sizeof(double));
	generate_mass(nb,ml,mh,m2);
	for (i=0;i<nb;++i) {
		j  = (ns - i)*randomz();
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		x[2*i+1].m = m2[i];
	}
	free(m2);
}

static
void PCRP(struct star *x, int n, int *nbin, double ml)
{
	int i,j;
	int nb,ns;
	double m1,m2;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		j  = (ns - i)*randomz();
		m1 = x[j].m;
		m2 = kroupa_IMF(ml,m1);
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		x[2*i+1].m = m2;
	}
}

static
void PCP_I(struct star *x, int n, int *nbin)
{
	int i,j;
	int nb,ns;
	double m1,m2;
	double q;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		j  = (ns - i)*randomz();
		m1 = x[j].m;
		q  = randomz();
		m2 = m1*q;
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		x[2*i+1].m = m2;
	}
}

static
void PCP_II(struct star *x, int n, int *nbin, double ml)
{
	int i,j;
	int nb,ns;
	double m1,m2;
	double q;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		j  = (ns - i)*randomz();
		m1 = x[j].m;
		q  = randomz();
		m2 = m1*q;
		if (m2 < ml) {
			nb--;
			continue;
		}
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		x[2*i+1].m = m2;
	}
	*nbin = nb;
}

static
void PCP_III(struct star *x, int n, int *nbin, double ml)
{
	int i,j;
	int nb,ns;
	double m1,m2;
	double q,q0;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		j  = (ns - i)*randomz();
		m1 = x[j].m;
		q0 = ml/m1;
		q  = hq(q0,1);
		m2 = m1*q;
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		x[2*i+1].m = m2;
	}
}

static
void SCP_I(struct star *x, int n, int *nbin, double ml, double mh)
{
	int i,j;
	int nb,ns;
	double Mc,m1,m2;
	double q;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		Mc = kroupa_IMF(2*ml,mh);
		q  = randomz();
		j  = (ns - i)*randomz();
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		m1 = Mc/(1+q);
		m2 = Mc/(1+1./q);
		x[2*i  ].m = m1;
		x[2*i+1].m = m2;
	}
}

static
void SCP_II(struct star *x, int n, int *nbin, double ml, double mh)
{
	int i,j;
	int nb,ns;
	double Mc,m1,m2;
	double q,q0;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		Mc = kroupa_IMF(2*ml,mh);
		q0 = ml/(Mc-ml);
		q  = randomz();
		if (q < q0) {
			nb--;
			continue;
		}
		j  = (ns - i)*randomz();
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		m1 = Mc/(1+q);
		m2 = Mc/(1+1./q);
		x[2*i  ].m = m1;
		x[2*i+1].m = m2;
	}
	*nbin = nb;
}

static
void SCP_III(struct star *x, int n, int *nbin, double ml, double mh)
{
	int i,j;
	int nb,ns;
	double Mc,m1,m2;
	double q,q0;
	nb = *nbin;
	ns = n - nb;
	for (i=0;i<nb;++i) {
		Mc = kroupa_IMF(2*ml,mh);
		q0 = ml/(Mc-ml);
		q  = hq(q0,1);
		j  = (ns - i)*randomz();
		SWAP(x[2*i],x[j],struct star);
		x[ns+i] = x[2*i+1];
		m1 = Mc/(1+q);
		m2 = Mc/(1+1./q);
		x[2*i  ].m = m1;
		x[2*i+1].m = m2;
	}
}

void binary_pairing(struct star *x, int n, int *nbin, double ml, double mh, int pairing_type)
{
	if (pairing_type==1) {
		/* PCRP */
		PCRP(x,n,nbin,ml);
	} else if (pairing_type == 2) {
		/* PCP-I */
		PCP_I(x,n,nbin);
	} else if (pairing_type == 3) {
		/* PCP-II */
		PCP_II(x,n,nbin,ml);
	} else if (pairing_type == 4) {
		/* PCP-III */
		PCP_III(x,n,nbin,ml);
	} else if (pairing_type == 5) {
		/* SCP-I */
		SCP_I(x,n,nbin,ml,mh);
	} else if (pairing_type == 6) {
		/* SCP-II */
		SCP_II(x,n,nbin,ml,mh);
	} else if (pairing_type == 7) {
		/* SCP-III */
		SCP_III(x,n,nbin,ml,mh);
	} else {
		/* RP */
		RP(x,n,nbin,ml,mh);
	}
}
