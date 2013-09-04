#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "constant.h"

#include "mempool.h"
#include "randomz.h"
#include "randqueue.h"

#define PART 8

void fractal(int StarNum, double D, double mlow, double mhigh, struct vector_s *star)
{
	unsigned long int LEN = sizeof(struct node);
#ifdef USE_MEMPOOL
	mempool_init(LEN);
#endif
	int N_div = 2;
	// The possibility of a sub-node matures to become a parent node.
	double mp = pow(N_div,D-3);
	fprintf(stderr,"The possibility of a sub-node matures is %lf\n",mp);

	double m;
	double temp;

	/* Kroupa's IMF, m_low>=0.08, m_high>0.5, m_break=0.5 */
	int s=2;
	double alpha[s];
	double mb[s+1];
	double xb;
	double norm;
	double norm0[s];
	double a1[s];
	alpha[0] = 1.3;
	alpha[1] = 2.3;
	a1[0] = 1 - alpha[0];
	a1[1] = 1 - alpha[1];
	mb[0] = mlow;
	mb[1] = 0.5;
	mb[2] = mhigh;
	/* check mass range */
	if (mlow > mhigh) {
		fprintf(stderr,"lower limit of mass %lf larger than the upper limit %lf!\n",mlow,mhigh);
		fprintf(stderr,"swap these two values...\n");
		temp = mlow;
		mlow = mhigh;
		mhigh = temp;
	}
	if (mlow < 0.08) {
		fprintf(stderr,"lower limit of mass should >= 0.08 M_sun!\n");
		exit(-1);
	}
	if (mhigh < 0.5) {
		fprintf(stderr,"upper limit of mass should > 0.5 M_sun!\n");
		exit(-1);
	}

	/* normalize and get weight for these two parts - xb */
	norm0[0] = -(pow(mb[0]/mb[1], a1[0]) - 1) / a1[0];
	norm0[1] =  (pow(mb[2]/mb[1], a1[1]) - 1) / a1[1];
	norm = 1./(norm0[0] + norm0[1]);
	xb = norm0[0]*norm;
	/* each part of IMF follows pow law distribution */
	/* normalization factor */
	double pow_norm1,pow_norm2;
	pow_norm1 = pow(mb[1]/mb[0], a1[0]) - 1;
	pow_norm2 = pow(mb[2]/mb[1], a1[1]) - 1;

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
	headc->next->x  = 0;
	headc->next->y  = 0;
	headc->next->z  = 0;
	headc->next->vx = 0;
	headc->next->vy = 0;
	headc->next->vz = 0;
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
	double r2;
	double delta,rnoise,vnoise;
	delta = 0.5;
	rnoise = 0.1;
	vnoise = 1.0;
	while (numc < 100*StarNum) {
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
			parent->x  = child->x;
			parent->y  = child->y;
			parent->z  = child->z;
			parent->vx = child->vx;
			parent->vy = child->vy;
			parent->vz = child->vz;
		}

		parent = headp->next;
		child  = headc;
		while (NULL != parent){
			numc--;
			for (j=0;j<PART;++j){
				if (randomz() < mp){

					// generate the mass of sub-nodes:
					if (randomz() < xb) {
						m = mb[0] * pow(pow_norm1*randomz()+1, 1/a1[0]);
					} else {
						m = mb[1] * pow(pow_norm2*randomz()+1, 1/a1[1]);
					}
					// m = makemass(mlow,mhigh);
						//	may be used some day...

					// generate the position of sub-nodes:
					x = parent->x;
					y = parent->y;
					z = parent->z;
					r2 = x*x + y*y + z*z;
					if ( j%2 < 1 ) x = parent->x + delta;
					else x = parent->x - delta;
					if ( j%4 < 2 ) y = parent->y + delta;
					else y = parent->y - delta;
					if ( j < 4 ) z = parent->z + delta;
					else z = parent->z - delta;

					// generate the velocity of sub-nodes:
					theta = -1.0 + 2.0*randomz();
					phi   = TWO_PI*randomz();

					switch (v_profile) {
						case '1':
							vnoise *= 1;
							//	v_sigma is a constant for every particle everywhere
							break;
						case '2':
							vnoise *= pow(1+r2,b/2);
							//	r=1, v_sigma=v1, so we set v_sigma(r)=<v_sigma>+(1+r^2)^b/2;
							break;
						default:
							vnoise *= 1;
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

					vx  = parent->vx*mcoefficient + vxy * cos(phi);
					vy  = parent->vy*mcoefficient + vxy * sin(phi);
					vz += parent->vz*mcoefficient;

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

					child->x  = gaussrand(x,rnoise);
					child->y  = gaussrand(y,rnoise);
					child->z  = gaussrand(z,rnoise);

					child->vx = vx;
					child->vy = vy;
					child->vz = vz;

					numc++;
//					printf("numc=%d\n",numc);
				}
			}
			parent = parent->next;
		}
		generation++;

		delta  /= 2;
		rnoise /= 2;
//		vnoise /= 2;
	}
	fprintf(stderr,"Generating finished!\n");
	fprintf(stderr,"%d generations, %d star candidates\n",generation,numc);
	fprintf(stderr,"delta, rnoise, vnoise: %lf %lf %lf\n",delta,rnoise,vnoise);

	struct vector_s *all_star;
	child = headc->next;
	all_star = (struct vector_s*) malloc (numc*sizeof(struct vector_s));
	if (NULL == all_star){
		fprintf(stderr,"fractal: malloc all_star failed!\n");
		exit(-1);
	}
	for (i=0;i<numc;++i){
		all_star[i].m  = child->m;

		all_star[i].x  = child->x;
		all_star[i].y  = child->y;
		all_star[i].z  = child->z;

		all_star[i].vx = child->vx;
		all_star[i].vy = child->vy;
		all_star[i].vz = child->vz;

		child = child->next;
	}

	int *idx,*staridx;
	if ( NULL == (idx = (int*) malloc (numc*sizeof(int))) ){
		fprintf(stderr,"fractal: malloc idx failed!\n");
		exit(-1);
	}
	if ( NULL == (staridx = (int*) malloc (numc*sizeof(int))) ){
		fprintf(stderr,"fractal: malloc staridx failed!\n");
		exit(-1);
	}
	for (i=0;i<numc;++i){
		idx[i] = i;
		staridx[i] = i;
	}
	fprintf(stderr,"Choose %d particles from %d leaves randomly\n",StarNum,numc);
	double eps=5e-4;
//	randqueue(numc,idx,StarNum,staridx);
	randqueue2(numc,idx,StarNum,staridx,all_star,eps);
	for (i=0;i<StarNum;++i) {
		star[i] = all_star[ staridx[i] ];
	}
	fprintf(stderr,"...random choose finished!\n");

	free(all_star);
	free(idx);
	free(staridx);

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
