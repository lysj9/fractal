#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "type.h"
#include "func.h"
#include "constant.h"

int main(int argc, char *argv[])
{
	int seed=0; // random seed (=0 choose from system time)
	int N_star=1000; // total star number
	int N_node; // node number (single stars + pairs of binaries)
	int nbin=0; // binary number
	double fbin=0.0; // binary fraction
	double rvir_nbody; // virial radius in NBODY-UNIT
	double r_virial=-1.0; // virial radius (in [pc])
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
	
	char param_str[]="s1n1r1R1t2f1u1l1D1b1B1q1";
//	char param_str[]="s1n1r1R1t2f1u1l1D1b1B1q1A1O1";
	int c;
	extern char *opt_str;
	extern char *opt_arr[16];
	while ( ( c=getopt(argc,argv,param_str) ) != -1 ){
		switch(c){
			case 's':
				seed = atoi(opt_str);
				break;
			case 'n':
				N_star = atoi(opt_str);
				break;
			case 'r':
				r_virial = atof(opt_str);
				break;
			case 'R':
				truncate = 0.5;
				r_tr = atof(opt_str);
//				r_hm = atof(opt_str);
				break;
			case 't':
				truncate = atof(opt_arr[0]);
				r_tr = atof(opt_arr[1]);
				break;
			case 'f':
				if ( (FP=fopen(opt_str,"w")) == NULL ){
					fprintf(stderr,"can not open %s\n",opt_str);
					exit(-1);
				}
				break;
			case 'l':
				mlow = atof(opt_str);
				break;
			case 'u':
				mhigh = atof(opt_str);
				break;
			case 'D':
				D = atof(opt_str);
				break;
			case 'b':
				nbin = atoi(opt_str);
				break;
			case 'B':
				fbin = atof(opt_str);
				break;
			case 'q':
				q = atof(opt_str);
				break;
			/*
			case 'A':
				dtadj = atof(opt_str);
				break;
			case 'O':
				deltat = atof(opt_str);
				break;
			*/
			default:	// should not happen
				printf("show help:\n");
				break;
		}
	}
	if (!FP){
		FP=fopen("fort.10","w");
		if (NULL==FP){
			fprintf(stderr,"initial file fort.10 open failed...\n");
			exit(-1);
//			FP=stdout;
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
	if (0==nbin){
		if ( fbin<0 || fbin>1 ){
			fprintf(stderr,"binary fraction wrong!\n");
			exit(1);
		}
		nbin = 0.5*N_star*fbin;
	}
	rvir_nbody = 2*(1-q);
	if (r_virial<0) {
		if (truncate<=0 || r_tr<=0) {
			r_virial=1.0;
			r_tr=-1.0;
			truncate=-1.0;
		}
	}
	// individual binaries + single stars + c.m. binaries
	struct vector_s *star = (struct vector_s*) malloc ( (N_star+nbin)*sizeof(struct vector_s) );
	if (NULL==star){
		fprintf(stderr,"malloc *star error\n");
		exit(-1);
	}
	double *mass = (double*) malloc ( (N_star+nbin)*sizeof(double) );
	if (NULL==mass){
		fprintf(stderr,"malloc *mass error\n");
		exit(-1);
	}

	int i,j;
	clock_t t_start,t_end;
	double t_cost;
	t_start=clock();
	randomz_seed(seed);

//	mass_pair(N_star,nbin,mass);
//	sort_mass();
	N_node = N_star - nbin;
	fractal (N_node,D,mlow,mhigh,star);

	t_end=clock();
	t_cost=(double)(t_end-t_start)/CLOCKS_PER_SEC;
	fprintf(stderr,"\ngenerate fractal cluster use %lf seconds...\n",t_cost);
	t_start=clock();

	double Mnode=0;
	double Mbin=0;
	for (i=0;i<N_node;++i) Mnode += star[i].m;
	// order data in NBODY data structure
	if (nbin>0){
		// binaries
		for (i=0;i<nbin;++i){
			// c.m of binaries
			star[i+N_star] = star[i];
		}
		// single stars
		for (i=0;i<N_node-nbin;++i){
			star[N_star-1-i] = star[N_node-1-i];
		}
		for (i=0;i<nbin;++i){
			// individual binaries
			j = nbin-1-i;
			star[2*j+1] = star[j];
			star[2*j] = star[j];
		}
		generate_mass(nbin,mlow,mhigh,mass);
		for (i=0;i<nbin;++i) Mbin += mass[i];
		// position in [pc], velocity in [m/s/VSC]
		// position and velocity are in binary frame
		generate_binaries(N_star,nbin,mass,star);
	}
	fprintf(stderr,"before eigenevolution: Mnode=%lf, Mbin=%lf, Mtot=%lf\n",Mnode,Mbin,Mnode+Mbin);

	// scale (node) to standard NBODY-UNIT,
	// i.e. total energy = -0.25, virial ratio = q
	int NNBMAX;
	double RS0;
	Mtot = nbody_scale(N_node,q,star+2*nbin,&NNBMAX,&RS0);
	m_mean = Mtot/N_star;

	// get the scale informations (rscale, vscale & tscale)
	double *r2_sort;
	r2_sort = (double*) malloc(N_node*sizeof(double));
	int *idx;
	idx = (int*) malloc(N_node*sizeof(int));
	sort_radius(N_node,star+2*nbin,r2_sort,idx);
	/*
	if (r_virial>0) {
		rscale = r_virial/rvir_nbody;
		r_hm_nbody = get_radius(N_star-nbin,star+2*nbin,0.5);
		r_hm = r_hm_nbody*rscale;
	} else {
		if (r_tr<=0) {
			r_virial = 1.0;
			rscale = r_virial/rvir_nbody;
			r_tr = r_tr_nbody*rscale;
		} else {
			rscale = r_tr/r_tr_nbody;
			r_virial = rvir_nbody*rscale;
		}
	}
	*/
	r_tr_nbody = get_radius(N_node,star+2*nbin,truncate,r2_sort,idx);
	if (0.5!=truncate) {
		r_hm_nbody = get_radius(N_node,star+2*nbin,0.5,r2_sort,idx);
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
//	tscale = sqrt(r_virial*r_virial*r_virial/Mtot)*TSC;

	// scale binaries to cluster frame, in NBODY-UNIT
	for (i=0;i<2*nbin;++i) {
		// mass in [MSUN]
		star[i].m /= Mtot;
		// position in [pc]
		star[i].x  = star[N_star+i/2].x + star[i].x/rscale;
		star[i].y  = star[N_star+i/2].y + star[i].y/rscale;
		star[i].z  = star[N_star+i/2].z + star[i].z/rscale;
		// generate_binaries returns velocity in [m/s/VSC]
		star[i].vx = star[N_star+i/2].vx + star[i].vx/v_virial;
		star[i].vy = star[N_star+i/2].vy + star[i].vy/v_virial;
		star[i].vz = star[N_star+i/2].vz + star[i].vz/v_virial;
	}

	t_end=clock();
	t_cost=(double)(t_end-t_start)/CLOCKS_PER_SEC;
	fprintf(stderr,"\nscale use %lf seconds...\n",t_cost);
	// output stars
	for (i=0;i<N_star+nbin;++i){
//		fprintf(FP,"%d %lf %lf %lf %lf %lf %lf %lf\n",
//				i,star[i].m,
		fprintf(FP,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n",
				star[i].m,
				star[i].x,star[i].y,star[i].z,
				star[i].vx,star[i].vy,star[i].vz);
	}
	fclose(FP);

	// output NBODY6 input files
//	output_nbody6(outname,N_star,nbin,r2_sort,seed,r_virial,m_mean,mlow,mhigh,q,&NNBMAX,&RS0);
	output_nbody6(outname,N_star,nbin,r2_sort,seed,r_virial,m_mean,mlow,mhigh,q,NNBMAX,RS0);

	// output read-in informations
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

	// output important informations to log file
	/* */
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
	/* */

	free(star);
	free(mass);

	return 0;
}
