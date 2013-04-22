#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void output_nbody6(char *outname, int N_star, int nbin, double *r2_sort, int seed, double r_virial, double m_mean, double mlow, double mhigh, double q, int *nnbmax_out, double *rs0_out)
{
	FILE *FP;
	if (NULL == (FP=fopen(outname,"w"))){
		fprintf(stderr,"output_nbody6: can not open %s\n",outname);
		exit(-1);
	}
//	int KSTART=1;
//	double TCOMP=5000000.0;
	fprintf(FP,"1 5000000.0\n");

	int N=N_star; // N=Ns+2*Nb;
//	int NFIX=1;
//	int NCRIT=10;
	int NRAND=seed;
	int NNBMAX;
	/**** NNBMAX ****/
	NNBMAX = 1.25*sqrt(N);
	if (NNBMAX<30) NNBMAX=30;
	if (N<=NNBMAX) NNBMAX=0.5*N;
	*nnbmax_out = NNBMAX;
	/**** NNBMAX ****/
//	int NRUN=1;
	fprintf(FP,"%d 1 10 %d %d 1\n",N,NRAND,NNBMAX);
//	fprintf(FP,"%d %d %d %d %d %d\n",N,NFIX,NCRIT,NRAND,NNBMAX,NRUN);

//	double ETAI=0.02;
//	double ETAR=0.02;
	double RS0;
	/**** RS0 ****/
	RS0 = sqrt(r2_sort[NNBMAX]);
	*rs0_out = RS0;
	/**** RS0 ****/
	double DTADJ=1.0;
	double DELTAT=1.0;
	double TCRIT=100.0;
//	double QE=1.0E-03;
	double RBAR=r_virial;
	double ZMBAR=m_mean;
	fprintf(FP,"0.02 0.02 %.4f %.4f %.4f %.4f 1.0E-03 %.8f %.8f\n",RS0,DTADJ,DELTAT,TCRIT,RBAR,ZMBAR);
//	fprintf(FP,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",ETAI,ETAR,RS0,DTADJ,DELTAT,TCRIT,QE,RBAR,ZMBAR);

//	int KZ[50];
//	int kz12=0;
//	int kz14=0;
	fprintf(FP,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(FP,"0 0 0 0 2 1 2 0 3 2\n");
	fprintf(FP,"0 4 2 0 1 2 0 1 0 1\n");
	fprintf(FP,"0 0 0 2 1 0 0 2 0 1\n");
	fprintf(FP,"0 0 0 0 0 0 0 0 0 0\n");

//	double DTMIN;
//	double RMIN;
//	double ETAU;
//	double ECLOSE;
//	double GMIN;
//	double GMAX;
	fprintf(FP,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.01\n");
//	fprintf(FP,"%lf %lf %lf %lf %lf %lf\n",DTMIN,RMIN,ETAU,ECLOSE,GMIN,GMAX);

//	double ALPHAS;
	double BODY1=mhigh;
	double BODYN=mlow;
	int NBIN0=nbin;
//	int NHI0;
	double ZMET=0.02;
//	double EPOCH0;
//	double DTPLOT=100.0;
	fprintf(FP,"2.350000 %.8f %.8f %d 0 %.8f 0.0 100.0\n",BODY1,BODYN,NBIN0,ZMET);
//	fprintf(FP,"%lf %lf %lf %d %d %lf %lf %lf\n",ALPHAS,BODY1,BODYN,NBIN0,NHI0,ZMET,EPOCH0,DTPLOT);

//	double Q=q;
//	double VXROT;
//	double VZROT;
//	double RTIDE;
//	double SMAX;
	fprintf(FP,"%.4f 0.0 0.0 0.00000 0.125\n",q);
//	fprintf(FP,"%lf %lf %lf %lf %lf\n",Q,VXROT,VZROT,RTIDE,SMAX);

	/*
	// output nbody6 information (to stderr)
	fprintf(stderr,"==== nbody6 information: ====\n");
	fprintf(stderr,"NNBMAX = %d\n",NNBMAX);
	fprintf(stderr,"RS0    = %lf [NBODY UNIT]\n",RS0);
	fprintf(stderr,"==== nbody6 information: ====\n");
	*/
	return;
}
