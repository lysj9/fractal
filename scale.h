#ifndef _SCALE_H_
#define _SCALE_H_

#include "type.h"

double nbody_scale(int N_node, double q, struct vector_s *star, int *nnbmax_out, double *rs0_out);
void sort_radius(int N_cm, struct vector_s *star, double *r2_sort, int *idx);
double get_radius(int N_cm, struct vector_s *star, double truncate, double *r2_sort, int *idx);

#endif
