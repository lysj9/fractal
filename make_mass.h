#ifndef _MAKE_MASS_H_
#define _MAKE_MASS_H_

double makemass(double ml, double mh, double x1, double x2, int s);
double kroupa_IMF(double ml, double mh);
double IMF(int s, double *a, double *m0);
double make_mass(double mlow, double mhigh);

#endif
