#ifdef __GNUC__
#include <sys/time.h>

#else

#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>

#endif

double get_wtime()
{
#ifdef __GNUC__
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
#else
#ifdef _OPENMP
	return omp_get_wtime();
#else
	clock_t t;
	t = clock();
	return 1.*t/CLOCKS_PER_SEC;
#endif
#endif
}
