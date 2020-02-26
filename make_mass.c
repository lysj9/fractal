#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "randomz.h"

/* ******** ******** ******** ******** ******** ******** ******** ******** */
/* get star mass from IMF */
/* ******** ******** ******** ******** ******** ******** ******** ******** */

double general_power_law(double ml, double mh, double alpha)
{
	double m;
	double norm;
	if (alpha == 1) {
		norm = log(mh/ml);
		m = ml * exp(norm*randomz());
	} else {
		norm = pow(mh/ml,-alpha+1) - 1;
		m = ml * pow(norm*randomz() + 1, 1/(1 - alpha));
	}
	return m;
}

double general_power_law2(double ml, double alpha, double norm)
{
	double m;
	if (alpha == 1) {
		m = ml * exp(norm*randomz());
	} else {
		m = ml * pow(norm*randomz() + 1, 1/(1 - alpha));
	}
	return m;
}

double get_mass(int s, double *ms, double *as)
{
	/* s sections, ml = ms[0], mh = ms[s], ms[i] < ms[i+1] */
	/* f(m)dm = ki * m^(-ai) , ms[i] <= m < ms[i+1], i=0,1,...,s-1 */
	double m=0;
	double norm;
	double norm0[s];
	double xs[s];
	double ns[s];
	double a1;
	double x;
	double ki;
	int i;
	a1 = 1 - as[0];
	ki = 1;
	norm0[0] = (pow(ms[1],a1) - pow(ms[0],a1)) / a1; /* norm0[0] *= ki */
	for (i=1;i<s;++i) {
		a1 = 1 - as[i];
		ki *= pow(ms[i],as[i]-as[i-1]);
		norm0[i] = ki * (pow(ms[i+1],a1) - pow(ms[i],a1)) / a1;
	}
	norm = 0;
	for (i=0;i<s;++i) norm += norm0[i];
	xs[0] = norm0[0]/norm;
	for (i=1;i<s;++i) xs[i] = xs[i-1] + norm0[i]/norm;
	
	for (i=0;i<s;++i) {
		if (as[i] == 1) {
			ns[i] = log(ms[i+1]/ms[i]);
		} else {
			ns[i] = pow(ms[i+1]/ms[i],-as[i]+1) - 1;
		}
	}

	x = randomz();
	for (i=0;i<s;++i) {
		if (x < xs[i]) {
			m = general_power_law2(ms[i],as[i],ns[i]);
			break;
		}
	}
	return m;
}

#define SECTIONS 16
static int g_s;
static double g_ms[SECTIONS+1];
static double g_as[SECTIONS];
static double g_xs[SECTIONS];
static double g_ns[SECTIONS];

void make_mass_init(int s, double *ms, double *as, double *xs, double *ns)
{
	/* s sections, ml = ms[0], mh = ms[s], ms[i] < ms[i+1] */
	/* f(m)dm = ki * m^(-ai) , ms[i] <= m < ms[i+1], i=0,1,...,s-1 */
	double norm;
	double norm0[s];
	double a1;
	double ki;
	int i;
	a1 = 1 - as[0];
	ki = 1;
	norm0[0] = (pow(ms[1],a1) - pow(ms[0],a1)) / a1; /* norm0[0] *= ki */
	for (i=1;i<s;++i) {
		a1 = 1 - as[i];
		ki *= pow(ms[i],as[i]-as[i-1]);
		norm0[i] = ki * (pow(ms[i+1],a1) - pow(ms[i],a1)) / a1;
	}
	norm = 0;
	for (i=0;i<s;++i) norm += norm0[i];
	xs[0] = norm0[0]/norm;
	for (i=1;i<s;++i) xs[i] = xs[i-1] + norm0[i]/norm;
	
	for (i=0;i<s;++i) {
		if (as[i] == 1) {
			ns[i] = log(ms[i+1]/ms[i]);
		} else {
			ns[i] = pow(ms[i+1]/ms[i],-as[i]+1) - 1;
		}
	}

	g_s = s;
	for (i=0;i<s;++i) {
		g_ms[i] = ms[i];
		g_as[i] = as[i];
		g_xs[i] = xs[i];
		g_ns[i] = ns[i];
	}
	g_ms[s] = ms[s];
}

double get_mass2(int s, double *ms, double *as, double *xs, double *ns)
{
	int i;
	double x;
	double m=0;
	x = randomz();
	for (i=0;i<s;++i) {
		if (x < xs[i]) {
			m = general_power_law2(ms[i],as[i],ns[i]);
			break;
		}
	}
	return m;
}

/* Kroupa IMF */
double kroupa_IMF(double ml, double mh)
{
	int i,j;
	int s=3;
	double ms[4] = {0,0.08,0.5,100};
	double as[3] = {0.3,1.3,2.3};
	for (i=s-1;i>0;--i) {
		if (mh > ms[i]) {
			break;
		}
	}
	s = i + 1;
	ms[s] = mh;

	if (ml < ms[0]) {
		ml = ms[0];
	} else {
		for (i=1;i<s;++i) {
			if (ml < ms[i]) {
				break;
			}
		}
		s = s - i + 1;
		ms[0] = ml;
		if (i>1) {
			for (j=1;j<s+1;++j) {
				ms[j] = ms[j+i-1];
			}
			for (j=0;j<s;++j) {
				as[j] = as[j+i-1];
			}
		}
	}
	return get_mass(s,ms,as);
}

/* Kroupa IMF (Kroupa, 2001, MNRAS 322, 231), a broken PLMF
   of the form

                     _
                    |
                    | m^(-1.3)    if (m_lower < m <= 0.5) 
       MF(m) = A * -|
                    | m^(-2.3)    if (0.5 < m <= m_upper) 
                    |_

   this mass function is a special case of the broken power-law MF
 */
/* star mass: 0.08 <= ml < 0.5 < mh 
 * mass range not change
 */
double kroupa_IMF2(double ml, double mh)
{
	double m;
	double mb=0.5;
	double norm0[2];
	double norm;
	static double alpha[2];
	static double a1[2];
	static double xb;
	static double pow_norm1,pow_norm2;
	static int first=1;
	if (first) {
		if (ml >= 0.5 || mh < 0.5) {
			fprintf(stderr,"kroupa_IMF2 should be used for Kroupa IMF "\
				   	"with mass range 0.08 <= ml < 0.5 < mh\n");
			exit(0);
		}
		alpha[0] = 1.3;
		alpha[1] = 2.3;
		a1[0] = 1 - alpha[0];
		a1[1] = 1 - alpha[1];
		norm0[0] = -(pow(ml/mb, a1[0]) - 1) / a1[0];
		norm0[1] =  (pow(mh/mb, a1[1]) - 1) / a1[1];
		norm = 1./(norm0[0] + norm0[1]);
		xb = norm0[0]*norm;
		pow_norm1 = pow(mb/ml, a1[0]) - 1;
		pow_norm2 = pow(mh/mb, a1[1]) - 1;
		first = 0;
	}
	if (randomz() < xb) {
		m = ml * pow(pow_norm1*randomz()+1, 1/a1[0]);
	} else {
		m = mb * pow(pow_norm2*randomz()+1, 1/a1[1]);
	}
	return m;
}

/* rejection sampling or acceptance-rejection method, inefficient way */
double kroupa_IMF3(double ml, double mh)
{
	double m;
	double m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double upper;
	if (ml<m1) upper = pow(ml,-a1);
	else if (ml>m2) upper = pow(ml,-a3);
	else upper = pow(ml,-a2);
	double dm = mh - ml;
	double temp;
	do {
		m = ml + dm * randomz();
		if (m<m1) {
			temp = pow(m1,a1-a2) * pow(m,-a1);
		} else if (m>m2) {
			temp = pow(m2,a3-a2) * pow(m,-a3);
		} else {
			temp = pow(m,-a2);
		}
	} while (temp > upper*randomz());
	return m;
}

/* ******** ******** ******** ******** ******** ******** ******** ******** */
/* generate mass distribution from IMF */
/* ******** ******** ******** ******** ******** ******** ******** ******** */

void gen_IMF(int n, double *m, int s, double *ms, double *as)
{
	/* s sections, ml = ms[0], mh = ms[s], ms[i] < ms[i+1] */
	/* f(m)dm = ki * m^(-ai) , ms[i] <= m < ms[i+1], i=0,1,...,s-1 */
	double norm;
	double norm0[s];
	double xs[s];
	double ns[s];
	double a1;
	double x;
	double ki;
	int i,j;
	a1 = 1 - as[0];
	ki = 1;
	norm0[0] = (pow(ms[1],a1) - pow(ms[0],a1)) / a1; /* norm0[0] *= ki */
	for (i=1;i<s;++i) {
		a1 = 1 - as[i];
		ki *= pow(ms[i],as[i]-as[i-1]);
		norm0[i] = ki * (pow(ms[i+1],a1) - pow(ms[i],a1)) / a1;
	}
	norm = 0;
	for (i=0;i<s;++i) norm += norm0[i];
	xs[0] = norm0[0]/norm;
	for (i=1;i<s;++i) xs[i] = xs[i-1] + norm0[i]/norm;
	
	for (i=0;i<s;++i) {
		if (as[i] == 1) {
			ns[i] = log(ms[i+1]/ms[i]);
		} else {
			ns[i] = pow(ms[i+1]/ms[i],-as[i]+1) - 1;
		}
	}

	for (j=0;j<n;++j) {
		x = randomz();
		for (i=0;i<s;++i) {
			if (x < xs[i]) {
				m[j] = general_power_law2(ms[i],as[i],ns[i]);
				break;
			}
		}
	}
}

void gen_kroupa_IMF(int n, double *m, double ml, double mh)
{
	double mb=0.5;
	double norm0[2];
	double norm;
	double as[2];
	double a1[2];
	double xb;
	double pow_norm1,pow_norm2;
	int i;

	if (ml >= 0.5 || mh < 0.5) {
		fprintf(stderr,"kroupa_IMF2 should be used for Kroupa IMF "\
			   	"with mass range 0.08 <= ml < 0.5 < mh\n");
		exit(0);
	}
	as[0] = 1.3;
	as[1] = 2.3;
	a1[0] = 1 - as[0];
	a1[1] = 1 - as[1];
	norm0[0] = -(pow(ml/mb, a1[0]) - 1) / a1[0];
	norm0[1] =  (pow(mh/mb, a1[1]) - 1) / a1[1];
	norm = 1./(norm0[0] + norm0[1]);
	xb = norm0[0]*norm;
	pow_norm1 = pow(mb/ml, a1[0]) - 1;
	pow_norm2 = pow(mh/mb, a1[1]) - 1;

	for (i=0;i<n;++i) {
		if (randomz() < xb) {
			m[i] = ml * pow(pow_norm1*randomz()+1, 1/a1[0]);
		} else {
			m[i] = mb * pow(pow_norm2*randomz()+1, 1/a1[1]);
		}
	}
}
