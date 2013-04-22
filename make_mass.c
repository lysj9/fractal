#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "func.h"

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
