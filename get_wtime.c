#ifdef __GNUC__
#include <sys/time.h>
#endif

double get_wtime()
{
#ifdef __GNUC__
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
#else
	clock_t t;
#ifdef _OPENMP
	t = omp_get_wtime();
#else
	t = clock();
#endif
	return 1.*t/CLOCKS_PER_SEC;
#endif
}
