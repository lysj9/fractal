#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "type.h"
#include "constant.h"

#include "frac.h"
#include "generate_binary.h"
#include "generate_mass.h"
#include "getopt.h"
#include "make_mass.h"
#include "output.h"
#include "randomz.h"
#include "scale.h"
#include "get_wtime.h"

int main(int argc, char *argv[])
{
	double t_start,t_end; 	// program time used
	double t_cost;

	double t0,t1; 			// function time used
	double tc;

	int i;

	int seed=0; // random seed (=0 choose from system time)
	int N_star=1000; // total star number
	int N_node; // node number (single stars + pairs of binaries)
	int nbin=0; // binary number
	double fbin=0.0; // binary fraction
	double rvir_nbody; // virial radius in NBODY-UNIT
	double r_virial=0; // virial radius (in [pc])
	double v_virial; // v_virial*VSC is the mean velocity (in [m/s])
					 // v_virial = sqrt(Mtot/r_virial);
	double truncate=-1.0;
	double r_tr_nbody=-1.0; // truncate (mass) radius in NBODY-UNIT
	double r_tr=-1.0;       // truncate (mass) mass radius (in [pc])
	truncate=0.5;
	double r_hm_nbody; // half-mass radius in NBODY-UNIT
	double r_hm=-1.0;  // half-mass radius (in [pc])
	// NBODY-UNIT to astrophysics unit (pc, m/s, Myr)
//	double mscale; // mscale = Mtot
	double rscale; // rscale = r_virial/rvir_nbody
	double vscale; // vscale = v_virial*VSC
	double tscale; // tscale = rscale/vscale*(PC/MYR)
				   // tscale = sqrt(r_virial^3/Mtot)*TSC
	double rs_estimated; // estimated rscale;
	FILE *FP=NULL;
	char *outname="a.input";
	char *outlog="a.log";
	double mlow=0.08,mhigh=50.0; // mass range
	double D=2.0; // fractal dimension
	double q=0.5; // virial ratio
	double Mtot;
	double m_mean;
//	double dtadj=1.0;
//	double deltat=1.0;
	int pairing_type=4; // binary pairing type, PCP-III by default
	
	char dic_str[]="s1:n1:r1:R1:t2:f1:l1:u1:D1:b1:B1:q1:p1";
	char arg_name[8];
	char *arg_str[16];
//	char param_str[]="s1n1r1R1t2f1u1l1D1b1B1q1A1O1";
	int c;

/* program starts: */
	t_start = get_wtime();

/*	******** ******** ******** ******** ******** ******** ********	*/
/* input parameters treatment */
/*	******** ******** ******** ******** ******** ******** ********	*/
	while ((c=getopt(argc,argv,arg_name,arg_str,dic_str,"",sizeof(arg_name))) != -1) {
		if (c==0) {
			if (strcmp(arg_name,"s") == 0) {
				seed = atoi(arg_str[0]);
			} else if (strcmp(arg_name,"n") == 0) {
				N_star = atoi(arg_str[0]);
			} else if (strcmp(arg_name,"r") == 0) {
				r_virial = atof(arg_str[0]);
			} else if (strcmp(arg_name,"R") == 0) {
				truncate = 0.5;
				r_tr = atof(arg_str[0]);
	//			r_hm = atof(arg_str[0]);
			} else if (strcmp(arg_name,"t") == 0) {
				truncate = atof(arg_str[0]);
				r_tr = atof(arg_str[1]);
			} else if (strcmp(arg_name,"f") == 0) {
				if ( (FP=fopen(arg_str[0],"w")) == NULL ){
					fprintf(stderr,"can not open %s\n",arg_str[0]);
					exit(-1);
				}
			} else if (strcmp(arg_name,"l") == 0) {
				mlow = atof(arg_str[0]);
			} else if (strcmp(arg_name,"u") == 0) {
				mhigh = atof(arg_str[0]);
			} else if (strcmp(arg_name,"D") == 0) {
				D = atof(arg_str[0]);
			} else if (strcmp(arg_name,"b") == 0) {
				nbin = atoi(arg_str[0]);
			} else if (strcmp(arg_name,"B") == 0) {
				fbin = atof(arg_str[0]);
			} else if (strcmp(arg_name,"q") == 0) {
				q = atof(arg_str[0]);
			} else if (strcmp(arg_name,"p") == 0) {
				pairing_type = atoi(arg_str[0]);
			} else {
				//printf("%d %s\n",c,arg_name);
				printf("show help:\n");
				break;
			}
		} else {
			printf("getopt end_stat = %d, arg_name = \"%s\"\n",c,arg_name);
			printf("show help:\n");
			exit(1);
		}
	}
	if (!FP) {
		FP=fopen("fort.10","w");
		if (NULL==FP){
			fprintf(stderr,"initial file fort.10 open failed...\n");
			exit(-1);
		}
	}
	
	if (N_star<2) {
		fprintf(stderr,"total star number should >= 2\n");
		exit(1);
	}
	if (nbin>N_star/2 || nbin<0){
		fprintf(stderr,"binary number wrong!\n");
		exit(1);
	}
	if (0==nbin) {
		if (fbin<0 || fbin>1) {
			fprintf(stderr,"binary fraction wrong!\n");
			exit(1);
		}
		nbin = 0.5*N_star*fbin;
	}
	rvir_nbody = 2*(1-q);
	if (r_virial<=0) {
		if (truncate<=0 || r_tr<=0) {
			r_virial=1.0;
			r_tr=-1.0;
			truncate=-1.0;
			rs_estimated=r_virial/rvir_nbody;
		}
		rs_estimated=r_tr/rvir_nbody;
	} else {
		rs_estimated=r_virial/rvir_nbody;
	}
	// individual binaries + single stars + c.m. binaries
	struct star *star_x = (struct star*) malloc ( (N_star+nbin)*sizeof(struct star) );
	if (NULL==star_x) {
		fprintf(stderr,"malloc *star_x error\n");
		exit(-1);
	}

	randomz_seed(seed);
/*	******** ******** ******** ******** ******** ******** ********	*/

/*	******** ******** ******** ******** ******** ******** ********	*/
/*	generate fractal distributed star cluster */
/*	******** ******** ******** ******** ******** ******** ********	*/
	t0 = get_wtime();
	N_node = N_star - nbin;
	fractal(N_node,D,mlow,mhigh,star_x,rs_estimated);

	t1 = get_wtime();
	tc = t1 - t0;
	fprintf(stderr,"\ngenerate fractal cluster use %lf seconds...\n",tc);
/*	******** ******** ******** ******** ******** ******** ********	*/

/*	******** ******** ******** ******** ******** ******** ********	*/
/*	generate binaries */
/*	******** ******** ******** ******** ******** ******** ********	*/
	t0 = get_wtime();
	// generate binaries and order data in NBODY data structure
	if (nbin > 0) {
		// position in [pc], velocity in [m/s/VSC]
		// position and velocity are in binary frame
		generate_binaries(star_x,N_star,nbin,mlow,mhigh,pairing_type);
	}

	t1 = get_wtime();
	tc = t1 - t0;
	fprintf(stderr,"\ngenerate binaries use %lf seconds...\n",tc);
/*	******** ******** ******** ******** ******** ******** ********	*/

/*	******** ******** ******** ******** ******** ******** ********	*/
/*	scale to NBODY-UNIT */
/*	******** ******** ******** ******** ******** ******** ********	*/
	t0 = get_wtime();
	// scale (node) to standard NBODY-UNIT,
	// i.e. total energy = -0.25, virial ratio = q
	int NNBMAX;
	double RS0;
	Mtot = nbody_scale(N_node,q,star_x+2*nbin,&NNBMAX,&RS0);
	m_mean = Mtot/N_star;

	// get the scale informations (rscale, vscale & tscale)
	double *r2_sort;
	r2_sort = (double*) malloc(N_node*sizeof(double));
	if (NULL==r2_sort) {
		fprintf(stderr,"malloc *r2_sort error\n");
		exit(-1);
	}
	int *idx;
	idx = (int*) malloc(N_node*sizeof(int));
	if (NULL==idx) {
		fprintf(stderr,"malloc *idx error\n");
		exit(-1);
	}
	sort_radius(N_node,star_x+2*nbin,r2_sort,idx);
	r_tr_nbody = get_radius(N_node,star_x+2*nbin,truncate,r2_sort,idx);
	if (0.5!=truncate) {
		r_hm_nbody = get_radius(N_node,star_x+2*nbin,0.5,r2_sort,idx);
	} else {
		r_hm_nbody = r_tr_nbody;
	}
	if (r_virial>0) {
		rscale = r_virial/rvir_nbody;
		r_hm = r_hm_nbody*rscale;
	} else {
		if (r_tr<=0) {
			r_virial = 1.0;
			rscale = r_virial/rvir_nbody;
			r_hm = r_hm_nbody*rscale;
		} else {
			rscale = r_tr/r_tr_nbody;
			r_virial = rvir_nbody*rscale;
			r_hm = r_hm_nbody*rscale;
		}
	}
	v_virial = sqrt(Mtot/r_virial); // v_virial*VSC in [m/s]
	vscale = v_virial*VSC;
	tscale = rscale/vscale*(PC/MYR);
	/* another way:
	 * tscale = sqrt(r_virial*r_virial*r_virial/Mtot)*TSC;
	 * */

	// scale binaries to cluster frame, in NBODY-UNIT
	for (i=0;i<2*nbin;++i) {
		// mass in [MSUN]
		star_x[i].m /= Mtot;
		// position in [pc]
		star_x[i].x[0] = star_x[N_star+i/2].x[0] + star_x[i].x[0] / rscale;
		star_x[i].x[1] = star_x[N_star+i/2].x[1] + star_x[i].x[1] / rscale;
		star_x[i].x[2] = star_x[N_star+i/2].x[2] + star_x[i].x[2] / rscale;
		// generate_binaries returns velocity in [m/s/VSC]
		star_x[i].x[3] = star_x[N_star+i/2].x[3] + star_x[i].x[3] / v_virial;
		star_x[i].x[4] = star_x[N_star+i/2].x[4] + star_x[i].x[4] / v_virial;
		star_x[i].x[5] = star_x[N_star+i/2].x[5] + star_x[i].x[5] / v_virial;
	}

	t1 = get_wtime();
	tc = t1 - t0;
	fprintf(stderr,"\nscale use %lf seconds...\n",tc);
/*	******** ******** ******** ******** ******** ******** ********	*/
	
/*	******** ******** ******** ******** ******** ******** ********	*/
/*	output */
/*	******** ******** ******** ******** ******** ******** ********	*/
	// output stars
	for (i=0;i<N_star+nbin;++i){
		fprintf(FP,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
				star_x[i].m,
				star_x[i].x[0],star_x[i].x[1],star_x[i].x[2],
				star_x[i].x[3],star_x[i].x[4],star_x[i].x[5]);
	}
	fclose(FP);

	// output NBODY6 input files
//	output_nbody6(outname,N_star,nbin,r2_sort,seed,r_virial,m_mean,mlow,mhigh,q,&NNBMAX,&RS0);
	output_nbody6(outname,N_star,nbin,r2_sort,seed,r_virial,m_mean,mlow,mhigh,q,NNBMAX,RS0);

	// output read-in informations
#ifdef LOG_INFO
	fprintf(stderr,"\n==== input information: ====\n");
	fprintf(stderr,"random seed=%d\n",abs(seed));
	fprintf(stderr,"%d stars in total, include %d single stars and %d pairs of binaries\n",N_star,N_star-2*nbin,nbin);
	fprintf(stderr,"virial ratio = %.4f\n",q);
	fprintf(stderr,"fractal dimension = %.4f\n",D);
	fprintf(stderr,"mass range: [%.4f,%.4f] [M_SUN]\n",mlow,mhigh);
	fprintf(stderr,"==== input information: ====\n");
	// output scaling information (to stderr)
	fprintf(stderr,"\n==== scaling information: ====\n");
	fprintf(stderr,"mass:     NBODY UNIT=1 ==> %e [M_SUN] (total mass)\n",Mtot);
	fprintf(stderr,"position: NBODY UNIT=1 ==> %e [pc] \n",rscale);
	fprintf(stderr,"velocity: NBODY UNIT=1 ==> %e [m/s]\n",vscale);
	fprintf(stderr,"time:     NBODY UNIT=1 ==> %e [Myr]\n",tscale);
	fprintf(stderr,"==== scaling information: ====\n");

	fprintf(stderr,"\n==== cluster information: ====\n");
	fprintf(stderr,"mean mass:        %lf [M_SUN]\n",m_mean);
	fprintf(stderr,"virial radius:    %lf [NBODY UNIT], %lf [pc]\n",rvir_nbody,r_virial);
	fprintf(stderr,"half-mass radius: %lf [NBODY UNIT], %lf [pc]\n",r_hm_nbody,r_hm);
	fprintf(stderr,"==== cluster information: ====\n");

	fprintf(stderr,"\n==== nbody6 information: ====\n");
	fprintf(stderr,"NNBMAX = %d\n",NNBMAX);
	fprintf(stderr,"RS0    = %lf [NBODY UNIT]\n",RS0);
	fprintf(stderr,"==== nbody6 information: ====\n");
#endif

	// output important informations to log file
	FP = fopen(outlog,"w");
	if (NULL == FP){
		fprintf(stderr,"%s open failed...\n",outlog);
		exit(-1);
	}
	fprintf(FP,"\n==== input information: ====\n");
	fprintf(FP,"random seed=%d\n",abs(seed));
	fprintf(FP,"%d stars in total, include %d single stars and %d pairs of binaries\n",N_star,N_star-2*nbin,nbin);
	fprintf(FP,"virial ratio = %.4f\n",q);
	fprintf(FP,"fractal dimension = %.4f\n",D);
	fprintf(FP,"mass range: [%.4f,%.4f] [M_SUN]\n",mlow,mhigh);
	fprintf(FP,"==== input information: ====\n");
	fprintf(FP,"\n==== scaling information: ====\n");
	fprintf(FP,"mass:     NBODY UNIT=1 ==> %e [M_SUN] (total mass)\n",Mtot);
	fprintf(FP,"position: NBODY UNIT=1 ==> %e [pc] \n",rscale);
	fprintf(FP,"velocity: NBODY UNIT=1 ==> %e [m/s]\n",vscale);
	fprintf(FP,"time:     NBODY UNIT=1 ==> %e [Myr]\n",tscale);
	fprintf(FP,"==== scaling information: ====\n");

	fprintf(FP,"\n==== cluster information: ====\n");
	fprintf(FP,"mean mass:        %lf [M_SUN]\n",m_mean);
	fprintf(FP,"virial radius:    %lf [NBODY UNIT], %lf [pc]\n",rvir_nbody,r_virial);
	fprintf(FP,"half-mass radius: %lf [NBODY UNIT], %lf [pc]\n",r_hm_nbody,r_hm);
	fprintf(FP,"==== cluster information: ====\n");

	fprintf(FP,"\n==== nbody6 information: ====\n");
	fprintf(FP,"NNBMAX = %d\n",NNBMAX);
	fprintf(FP,"RS0    = %lf [NBODY UNIT]\n",RS0);
	fprintf(FP,"==== nbody6 information: ====\n");
	fclose(FP);
/*	******** ******** ******** ******** ******** ******** ********	*/

/*	******** ******** ******** ******** ******** ******** ********	*/
/*	memory management */
/*	******** ******** ******** ******** ******** ******** ********	*/
	free(r2_sort);
	free(idx);
	free(star_x);

	t_end = get_wtime();
	t_cost = t_end - t_start;
	fprintf(stderr,"\nin total: use %lf ms, %lf seconds...\n",1.e3*t_cost,t_cost);
/*	******** ******** ******** ******** ******** ******** ********	*/
	return 0;
}
