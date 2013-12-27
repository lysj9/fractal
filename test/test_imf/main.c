#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "get_wtime.h"
#include "make_mass.h"
#include "randomz.h"

#define MAX_NS_IMF 10

int main()
{
	int nsimf = 3;
	double ms[MAX_NS_IMF+1] = {0.01,0.08,0.5,100.0,100.0};
	double as[MAX_NS_IMF] = {0.3,1.3,2.3};
	int i;
	int n=100000;
	double m;
	double xs[nsimf];
	double ns[nsimf];

	double t_start,t_end,t_cost;

	t_start = get_wtime();
	make_mass_init(nsimf,ms,as,xs,ns);
	for (i=0;i<n;++i) {
		//m = get_mass(nsimf,ms,as);
		//m = get_mass2(nsimf,ms,as,xs,ns);
		m = kroupa_IMF(0.01,100.0);
		printf("%lf\n",log10(m));
	}
	t_end = get_wtime();

	t_cost = t_end - t_start;
	fprintf(stderr,"%lf s\n",t_cost);

	return 0;
}
