#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "constant.h"

#include "get_wtime.h"
#include "make_mass.h"
#include "mempool.h"
#include "randomz.h"
#include "randqueue.h"

#define PART 8

void fractal2(int StarNum, int n_sub, struct star *star_x, double D, int nsimf, double *ms, double *as, double rs_estimated)
{
	unsigned long int LEN = sizeof(struct node);
#ifdef USE_MEMPOOL
	mempool_init(LEN);
#endif
	int N_div = 2;
	// The possibility of a sub-node matures to become a parent node.
	double mp = pow(N_div,D-3);
	fprintf(stderr,"The possibility of a sub-node matures is %lf\n",mp);

	// sub-cluster should contain stars more than sqrt(N)
	int generation_max = log(sqrt(StarNum)) / log(8*mp);

	double m;

	double xs[nsimf];
	double ns[nsimf];
	make_mass_init(nsimf,ms,as,xs,ns);

	int numc=1;
	struct node *child,*parent,*headc,*headp;
	if (NULL == (headc = (struct node*) malloc(LEN))) {
		fprintf(stderr,"fractal: malloc headc failed!\n");
		exit(-1);
	}
#ifdef USE_MEMPOOL
	headc->next = (struct node*) mempool_malloc();
#else
	headc->next = (struct node*) malloc(LEN);
	if (NULL == headc->next){
		fprintf(stderr,"fractal: malloc headc failed!\n");
		exit(-1);
	}
#endif
	headc->next->m  = 1;
	headc->next->x[0] = 0;
	headc->next->x[1] = 0;
	headc->next->x[2] = 0;
	headc->next->x[3] = 0;
	headc->next->x[4] = 0;
	headc->next->x[5] = 0;
	headc->next->next = NULL;

#ifdef USE_MEMPOOL
	headp = (struct node*) mempool_malloc();
#else
	headp = (struct node*) malloc(LEN);
	if (NULL == headc->next){
		fprintf(stderr,"fractal: malloc headp failed!\n");
		exit(-1);
	}
#endif
	headp->next = NULL;

	int i,j;
	int generation=0;
	/* velocity profile */
	int v_profile=1;
	double b = -0.5;		/* pow-index of velocity-radius; */
	/* coherent of velocity/momentum/kinetic */
	int coherent_type=1;
	double mcoefficient=1;	/* coherent coefficient */
	double x,y,z,vx,vy,vz;
	double v,vxy,theta,phi;
	double v0; 				/* velocity dispersiom */
	double r2;
	double delta,rnoise,vnoise;
	delta = 0.5;
	rnoise = 0.1;
	v0 = 1.0;
	vnoise = 1.0;

	int n_cr;
	n_cr = 100*sqrt(StarNum);
//	n_cr = 100*n_sub;

	while (numc < n_cr) {
		parent = headp;
		child  = headc;
		while (NULL != child->next) {
			child = child->next;
			if (NULL == parent->next) {
#ifdef USE_MEMPOOL
				parent->next = (struct node*) mempool_malloc();
#else
				parent->next = (struct node*) malloc(LEN);
				if (NULL == parent->next){
					fprintf(stderr,"fractal: generation=%d\t",generation);
					fprintf(stderr,"malloc parent->next failed!\n");
					exit(-1);
				}
#endif
				parent->next->next = NULL;
			}
			parent = parent->next;

			parent->m  = child->m;
			parent->x[0] = child->x[0];
			parent->x[1] = child->x[1];
			parent->x[2] = child->x[2];
			parent->x[3] = child->x[3];
			parent->x[4] = child->x[4];
			parent->x[5] = child->x[5];
		}

		parent = headp->next;
		child  = headc;
		while (NULL != parent){
			numc--;
			for (j=0;j<PART;++j){
				if (randomz() < mp){

					// generate the mass of sub-nodes:
					m = get_mass2(nsimf,ms,as,xs,ns);

					// generate the position of sub-nodes:
					x = parent->x[0];
					y = parent->x[1];
					z = parent->x[2];
					r2 = x*x + y*y + z*z;
					if ( j%2 < 1 ) x = parent->x[0] + delta;
					else x = parent->x[0] - delta;
					if ( j%4 < 2 ) y = parent->x[1] + delta;
					else y = parent->x[1] - delta;
					if ( j < 4 ) z = parent->x[2] + delta;
					else z = parent->x[2] - delta;

					// generate the velocity of sub-nodes:
					theta = -1.0 + 2.0*randomz();
					phi   = TWO_PI*randomz();

					switch (v_profile) {
						case '1':
							//vnoise = v0;
							//	v_sigma is a constant for every particle everywhere
							break;
						case '2':
							vnoise = v0 * pow(1+r2,b/2);
							//	r=1, v_sigma=v1, so we set v_sigma(r)=<v_sigma>+(1+r^2)^b/2;
							break;
						default:
							//vnoise = v0;
							break;
					}
					v = gaussrand(0.0,vnoise);
					vz  = v*theta;
					vxy = sqrt(v*v - vz*vz);
					switch (coherent_type) {
						case '1':
							mcoefficient = 1;                   /* coherent of velocity */
							break;
						case '2':
							mcoefficient = parent->m / m;       /* coherent of momentum */
							break;
						case '3':
							mcoefficient = sqrt(parent->m / m); /* coherent of kinetic */
							break;
						default:
							mcoefficient = 1;
							break;
					}

					vx  = parent->x[3]*mcoefficient + vxy * cos(phi);
					vy  = parent->x[4]*mcoefficient + vxy * sin(phi);
					vz += parent->x[5]*mcoefficient;

					if (NULL == child->next) {
#ifdef USE_MEMPOOL
						child->next = (struct node*) mempool_malloc();
#else
						child->next = (struct node*) malloc(LEN);
						if (NULL == child->next){
							fprintf(stderr,"fractal: generation=%d\t",generation);
							fprintf(stderr,"malloc child->next failed!\n");
							exit(-1);
						}
#endif
						child->next->next = NULL;
					}
					child = child->next;

					child->m  = m;

					child->x[0]  = gaussrand(x,rnoise);
					child->x[1]  = gaussrand(y,rnoise);
					child->x[2]  = gaussrand(z,rnoise);

					child->x[3] = vx;
					child->x[4] = vy;
					child->x[5] = vz;

					numc++;
//					printf("numc=%d\n",numc);
				}
			}
			parent = parent->next;
		}
		generation++;

		delta  /= 2;
		rnoise /= 2;
		vnoise /= 2;
	}
	fprintf(stderr,"Generating finished!\n");
	fprintf(stderr,"%d generations, %d star candidates\n",generation,numc);
	fprintf(stderr,"delta, rnoise, vnoise: %lf %lf %lf\n",delta,rnoise,vnoise);

	struct star *all_star;
	child = headc->next;
	all_star = (struct star*) malloc (numc*sizeof(struct star));
	if (NULL == all_star){
		fprintf(stderr,"fractal: malloc all_star failed!\n");
		exit(-1);
	}
	for (i=0;i<numc;++i){
		all_star[i].m  = child->m;

		all_star[i].x[0] = child->x[0];
		all_star[i].x[1] = child->x[1];
		all_star[i].x[2] = child->x[2];

		all_star[i].x[3] = child->x[3];
		all_star[i].x[4] = child->x[4];
		all_star[i].x[5] = child->x[5];

		child = child->next;
	}

	int *idx,*sidx;
	if ( NULL == (idx = (int*) malloc (numc*sizeof(int))) ){
		fprintf(stderr,"fractal: malloc idx failed!\n");
		exit(-1);
	}
	if ( NULL == (sidx = (int*) malloc (numc*sizeof(int))) ){
		fprintf(stderr,"fractal: malloc sidx failed!\n");
		exit(-1);
	}
	for (i=0;i<numc;++i){
		idx[i] = i;
		sidx[i] = i;
	}
	fprintf(stderr,"Choose %d particles from %d leaves randomly\n",StarNum,numc);
	double t0,t1;
	double tc;
	t0 = get_wtime();
//	double eps=5e-4;
	double eps=50./206265*rs_estimated; // minimum separation set to 50 AU
//	randqueue(numc,idx,StarNum,sidx);
	randqueue2(numc,idx,StarNum,sidx,all_star,eps);
	for (i=0;i<StarNum;++i) {
		star_x[i] = all_star[ sidx[i] ];
	}

	randqueue(numc,idx,n_sub,sidx);
//	randqueue2(numc,idx,StarNum,sidx,all_star,eps);

	for (i=0;i<n_sub;++i) {
		make_king(nstar,x);
		for (j=0;j<nstar;++j) {
			star_x[idx] = x[j] + xs[i];
		}
	}

	t1 = get_wtime();
	tc = t1 - t0;
	fprintf(stderr,"restrict randomly choosing use %lf ms, %lf s\n",1.e3*tc,tc);
	fprintf(stderr,"...random choose finished!\n");

	free(all_star);
	free(idx);
	free(sidx);

#ifdef USE_MEMPOOL
	mempool_destroy();
#else
	struct node* pf;
	child = headc;
	while (NULL != child->next){
		pf = child->next;
		free(child);
		child = pf;
	}
	free(child);
	parent = headp;
	while (NULL != parent->next){
		pf = parent->next;
		free(parent);
		parent = pf;
	}
	free(parent);
#endif

	return;
}
