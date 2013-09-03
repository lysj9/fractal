#ifndef _RANDQUEUE_H_
#define _RANDQUEUE_H_

#include "type.h"

void randqueue2(int n, int *idx, int randnum, int *randidx, struct vector_s *star, double eps);

void randqueue( int        n,
				int     *idx,
				int  randnum,
				int *randidx );

#endif
