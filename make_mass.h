#ifndef _MAKE_MASS_H_
#define _MAKE_MASS_H_

double general_power_law(double ml, double mh, double alpha);
double general_power_law2(double ml, double alpha, double norm);
double get_mass(int s, double *ms, double *as);
void make_mass_init(int s, double *ms, double *as, double *xs, double *ns);
double get_mass2(int s, double *ms, double *as, double *xs, double *ns);
double kroupa_IMF(double ml, double mh);
double kroupa_IMF2(double ml, double mh);
double kroupa_IMF3(double ml, double mh);

void gen_IMF(int n, double *m, int s, double *ms, double *as);
void gen_kroupa_IMF(int n, double *m, double ml, double mh);

#endif
