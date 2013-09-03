#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"

#include "randomz.h"

void generate_mass(int N, double mlow, double mhigh, double *mass)
{
	int i;
	double m,m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double upper;
	if (mlow<m1) upper = pow(mlow,-a1);
	else if (mlow>m2) upper = pow(mlow,-a3);
	else upper = pow(mlow,-a2);
//	three part pow-law mass function;
	double dm = mhigh - mlow;
	double temp;
	if (mlow>mhigh){
		fprintf(stderr,"lower limit of mass %lf larger than the upper limit %lf!\n",
				mlow,mhigh);
		fprintf(stderr,"swap these two values...\n");
		temp = mlow;
		mlow = mhigh;
		mhigh = temp;
	}
	for (i=0;i<N;++i){
		do {
			m = mlow + dm * randomz();
			if ( m<m1 ) {
				temp = pow( m1,a1-a2 ) * pow( m,-a1 );
			} else if ( m>m2 ) {
				temp = pow( m2,a3-a2 ) * pow( m,-a3 );
			} else {
				temp = pow( m,-a2 );
			}
		} while (upper*randomz() > temp);
		mass[i] = m;
	}
	return;
}
