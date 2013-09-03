#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"

#include "randomz.h"

static
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

double makemass(double ml, double mh, double x1, double x2, int s)
{
	double m,m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double x;
	switch (s) {
		case '3':
			x = randomz();
			if (x < x1) {
				m = general_power_law(ml,m1,a1);
			} else if (x < x2) {
				m = general_power_law(m1,m2,a2);
			} else {
				m = general_power_law(m2,mh,a3);
			}
			break;
		case '2':
			if (randomz() < x1) {
				m = general_power_law(ml,m2,a2);
			} else {
				m = general_power_law(m2,mh,a3);
			}
			break;
		case '1':
			general_power_law(ml,mh,a3);
			break;
	}
}

/* Kroupa IMF */
double kroupa_IMF(double ml, double mh)
{
	double m,m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double a11=1-a1,a21=1-a2,a31=1-a3;
	double x,x1,x2;
	double norm1=0,norm2=0,norm3=0,norm;
	if (ml < m1) {
		if (mh < m1) {
			m = general_power_law(ml,mh,a1);
		} else if (mh < m2) {
			norm1 = (pow(m1,a11) - pow(ml,a11))/a11;
			norm2 = (pow(mh,a21) - pow(m1,a21))/a21;
			norm  = 1.0 / (norm1 + norm2);
			x1 = norm * norm1;
			x  = randomz();
			if (x < x1) {
				m = general_power_law(ml,m1,a1);
			} else {
				m = general_power_law(m1,mh,a2);
			}
		} else {
			norm1 = (pow(m1,a11) - pow(ml,a11))/a11;
			norm2 = (pow(m2,a21) - pow(m1,a21))/a21;
			norm3 = (pow(mh,a31) - pow(m2,a31))/a31;
			norm  = 1.0 / (norm1 + norm2 + norm3);
			x1 = norm * norm1;
			x2 = norm * (norm1 + norm2);
			x  = randomz();
			if (x < x1) {
				m = general_power_law(ml,m1,a1);
			} else if (x < x2) {
				m = general_power_law(m1,m2,a2);
			} else {
				m = general_power_law(m2,mh,a3);
			}
		}
	} else if (ml<m2) {
		if (mh < m2) {
			m = general_power_law(ml,mh,a2);
		} else {
			norm2 = (pow(m2,a21) - pow(ml,a21))/a21;
			norm3 = (pow(mh,a31) - pow(m2,a31))/a31;
			norm  = 1.0 / (norm2 + norm3);
			x2 = norm * norm2;
			if (x < x2) {
				m = general_power_law(ml,m2,a2);
			} else {
				m = general_power_law(m2,mh,a3);
			}
		}
	} else {
		m = general_power_law(ml,mh,a3);
	}
	return m;
}

double IMF(int s, double *a, double *m0)
{
	/* s sections, ml = m0[0], mh = m0[s], m0[i] < m0[i+1] */
	/* f(m)dm = ki * m^(-ai) , m0[i] < m < m0[i+1], i=0,1,...,s-1 */
	double m;
	double norm;
	double norm0[s];
	double a1;
	double x0[s];
	double x;
	double ki;
	int i;
	a1 = 1 - a[0];
	ki = 1;
	norm0[0] = (pow(m0[1],a1) - pow(m0[0],a1)) / a1; /* norm0[0] *= ki */
	for (i=0;i<s;++i) {
		a1 = 1 - a[i];
		ki *= pow(m0[i],a[i]-a[i-1]);
		norm0[i] = ki * (pow(m0[i+1],a1) - pow(m0[i],a1)) / a1;
	}
	
	norm = 0;
	for (i=0;i<s;++i) norm += norm0[i];
	x0[i] = norm0[0]/norm;
	for (i=1;i<s;++i) x0[i] = x0[i-1] + norm0[i]/norm;
	x = randomz();
	for (i=0;i<s;++i) {
		if (x < x0[i]) {
			m = general_power_law(m0[i],m0[i+1],a[i]);
			break;
		}
	}
	return m;
}

double make_mass(double mlow, double mhigh)
{
	double m;
	double m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double upper;
	if (mlow<m1) upper = pow(mlow,-a1);
	else if (mlow>m2) upper = pow(mlow,-a3);
	else upper = pow(mlow,-a2);
	double dm = mhigh - mlow;
	double temp;
//	three part pow-law mass function;
	do {
		m = mlow + dm * randomz();
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
