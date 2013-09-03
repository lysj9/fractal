#ifndef _QUICK_SORT_H_
#define _QUICK_SORT_H_

void quick_sort_noidx_recursive(double *a, int l, int r);
void quick_sort_widx_recursive(double *a, int *idx, int l, int r);
void quick_sort_noidx_loop(double *a, int n);
void quick_sort_widx_loop(double *a, int *idx, int n);


/* sequential version */
void quick_sort_noidx_seq(double *a, int n);
void quick_sort_widx_seq(double *a, int *idx, int n);

/* openmp version */
void quick_sort_noidx_omp(double *a, int n);
void quick_sort_widx_omp(double *a, int *idx, int n);

/* global interface */
void quick_sort1(double *a, int n);
void quick_sort2(double *a, int *idx, int n);

#endif
