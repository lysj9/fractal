#ifndef _RANDOMZ_H_
#define _RANDOMZ_H_

void randomz_seed(int seed);
float randomz(void);
double randomz_dbl(void);
float normal_rand(void);
float gaussrand(float mean, float sigma);
double normal_rand_dbl(void);
double gaussrand_dbl(double mean, double sigma);

#endif
