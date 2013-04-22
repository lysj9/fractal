#include <math.h>
#include <time.h>

#define A 16807
#define M 2147483647
#define Q 127773
#define R 2836

//int z = (int) time(NULL);
static int z=1;
static int ix=245371421,iy=16807;

void randomz_seed( int seed )
{
	double k;
	if ( seed ==0 ){
		z = (int) time(NULL);
	}else{
		z = fabs(seed);
	}

	iy = z;
	k = iy/Q;
	iy = A*(iy-k*Q) - R*k;
	if ( iy<0 ) iy += M;

	ix = iy;
	ix ^= (ix<<13);
	ix ^= (ix>>17);
	ix ^= (ix<<5);
	return;
}

double randomz()
{
	double k,z;
	ix ^= (ix<<13);
	ix ^= (ix>>17);
	ix ^= (ix<<5);

	k = iy/Q;
	iy = A*(iy-k*Q) - R*k;
	if ( iy<0 ) iy += M;

	z = (1.0/M) * ( ( M & ( ix ^ iy ) ) | 1 );
	return z;
}

#ifdef BIG
#undef BIG
#endif
#define BIG 10
double gaussrand1( double mean,double sigma )
{
	double x = -BIG + 2*BIG*randomz();
	while (1){
		if ( randomz() < exp(-0.5*x*x) )
			return mean+sigma*x;
	}
}

#ifdef PI
#undef PI
#endif
#define PI 3.14159265358979323846
double gaussrand2( double mean,double sigma )
{
	double x,y,z;
	x = randomz();
	y = 2*PI*randomz();
//	y = randomz()*acos(0)*4;
	z = sqrt( -2*log(x) ) * cos(y);
	return mean + z*sigma;
}
