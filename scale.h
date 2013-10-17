#ifndef _SCALE_H_
#define _SCALE_H_

#include "type.h"

double nbody_scale(int N_node, double q, struct star *s, int *nnbmax_out, double *rs0_out);
void sort_radius(int N_cm, struct star *s, double *r2_sort, int *idx);
double get_radius(int N_cm, struct star *s, double truncate, double *r2_sort, int *idx);

#endif
