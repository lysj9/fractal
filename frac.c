#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "func.h"
#include "constant.h"

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
	// three part pow-law mass function;
	double m,m1=0.08,m2=0.5;
	double a1=0.3,a2=1.3,a3=2.3;
	double upper;
	double temp;
	if (mlow>mhigh) {
		fprintf(stderr,"lower limit of mass %lf larger than the upper limit %lf!\n",
				mlow,mhigh);
		fprintf(stderr,"swap these two values...\n");
		temp = mlow;
		mlow = mhigh;
		mhigh = temp;
	}
	if (mlow<m1) {
		upper = pow(mlow,-a1);
	} else if (mlow>m2) {
		upper = pow(mlow,-a3);
	} else {
		upper = pow( mlow,-a2 );
	}
	double dm = mhigh - mlow;
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


	int i,j,generation=0;
	double b = -0.5;		//	pow-index of velocity-radius;
	double mcoefficient=1;	//	the coherent of velocity/momentum/kinatic;
	double x,y,z,vx,vy,vz;
	double r2,v,v0,v1,vxy,theta,phi;
//	double radius2;
	double delta,rnoise,vnoise;
	v0 = 1.0;
	v1 = pow(1+1,b/2);
//	radius2 = 1.0*1.0;
	delta = 0.5;
	rnoise = 0.1;
	vnoise = 1.0;
	while (numc<100*StarNum){
		parent = headp;
		child  = headc;
		while (NULL != child->next ){
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
					do {
						m = mlow + dm * randomz();
						if (m<m1){
							temp = pow(m1,a1-a2) * pow(m,-a1);
						} else if (m>m2) {
							temp = pow(m2,a3-a2) * pow(m,-a3);
						} else {
							temp = pow(m,-a2);
						}
					} while (upper*randomz() > temp);
//					m = mass[i];
					// m = makemass( mlow,mhigh );
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
					v0 = 2*v1;
						//	v0 is a constant for every particle everywhere
//					v0 = pow( 1+r2,b/2 )+v1;
						//	r=1, v_sigma=v1, so we set v_sigma(r)=<v_sigma>+(1+r^2)^b/2;
//					vnoise = pow( 1+r2,b/2 );
						//	r=1, v_sigma=v1, so we set v_sigma(r)=<v_sigma>+(1+r^2)^b/2;
					v = gaussrand(0.0,vnoise);
						//	v is the velocity disperion;
					vz  = v*theta;
					vxy = sqrt(v*v - vz*vz);
//					mcoefficient = 1;					//	coherent of velocity
//					mcoefficient = parent->m/m;			//	coherent of momentum
//					mcoefficient = sqrt(parent->m/m);	//	coherent of kinetic

					vx  = parent->vx*mcoefficient + vxy * cos(phi);
					vy  = parent->vy*mcoefficient + vxy * sin(phi);
					vz += parent->vz*mcoefficient;
//					vx  = parent->vx*mcoefficient + randomz(-1,1);
//					vy  = parent->vy*mcoefficient + randomz(-1,1);
//					vz  = parent->vz*mcoefficient + randomz(-1,1);

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

//					child->x  = x+rnoise*randomz(-1,1);
//					child->y  = y+rnoise*randomz(-1,1);
//					child->z  = z+rnoise*randomz(-1,1);
					child->x  = x+gaussrand(0,rnoise);
					child->y  = y+gaussrand(0,rnoise);
					child->z  = z+gaussrand(0,rnoise);

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
//		if ( generation<8 ){
//			delta  /= 2;
//			rnoise /= 2;
//		}

		delta  /= 2;
		rnoise /= 2;
		vnoise /= 2;
	}
	fprintf(stderr,"Generating finished!\n");
	fprintf(stderr,"%d generations\n",generation);
	fprintf(stderr,"delta, rnoise, vnoise: %lf %lf %lf\n",delta,rnoise,vnoise);

	struct vector_s *all_star;
	child = headc->next;
	all_star = (struct vector_s*) malloc ( numc*sizeof(struct vector_s) );
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
	if ( NULL == (idx = (int*) malloc ( numc*sizeof(int) )) ){
		fprintf(stderr,"fractal: malloc idx failed!\n");
		exit(-1);
	}
	if ( NULL == (staridx = (int*) malloc ( numc*sizeof(int) )) ){
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
