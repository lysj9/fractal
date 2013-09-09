#if defined(__unix__)
#include <stdio.h>
#include <sys/time.h>
#elif defined(_WIN32)
#include <windows.h>
#elif defined(_OPENMP)
#include <omp.h>
#else
#include <time.h>
#endif

double get_wtime(void)
{
#if defined(__unix__)
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
#elif defined(_WIN32)
	LARGE_INTEGER frequency;
	LARGE_INTEGER t;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&t);
	return 1. * t.QuadPart / frequency.QuadPart;
#elif defined(_OPENMP)
	return omp_get_wtime();
#else
	return 1. * clock() / CLOCKS_PER_SEC;
#endif
}
