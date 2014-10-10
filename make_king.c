#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358

double densty(double z);

//int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking, double D, int symmetry)
int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking)
{
	
	//ODE variables
	int M = 10001;				//Number of interpolation points
	int KMAX = 10000;			//Maximum number of output steps of the integrator
	
	int i,j,k;
	double h;
	double den;
	double xstart, ystart0, ystart1;
	double x1, x2;
	double xp[KMAX], x[M];	
	int kount = 0;
	double yking[M][2], mass[M];
	double rmin = 0.2;
	
	double **yp;
	yp = (double **)calloc(KMAX,sizeof(double *));
	for (i=0;i<KMAX;i++){
		yp[i] = (double *)calloc(2,sizeof(double));
		if (yp[i] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	double pot = 0.0;
	double totmas;
	double zh;
	
	double ve, cg1;
	double fmass, r2;
	double costh, sinth, phi, r, xstar, ystar, zstar;
	double w1, w, wj1, wj;
	double vmax, speed, fstar, ustar, vstar, wstar, mstar;
	double ri;
	
	double **coord;
	coord = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		coord[j] = (double *)calloc(6,sizeof(double));
		if (coord[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}

	
	
	if (W0>12.0) {
		printf("W0 too large\n");
		return 0;
	} else if (W0 < 0.2) {
		printf("W0 too small\n");
		return 0;
	}
	
	
	
	
	/***************************
	 * INTERPOLATE KING VALUES *
	 ***************************/
	
	h = pow(W0-rmin,2)/(1.0*M-1.0);	//step size for interpolation
	den = densty(W0);			//central density
	
	//x is King's W, y[0] is z**2, y[1] is 2*z*dz/dW, where z is scaled radius
	//so x(y[0]) is energy as function of radius^2 and x(y[1]) is the derivative
	
	xstart = 0.000001;
	ystart0 = 2.0*xstart/3.0;
	ystart1 = -2.0/3.0;
	
	x1 = W0 - xstart;
	x2 = 0.0;
	
	//integrate Poisson's eqn
	printf("Integrating Poisson's equation\n");

	odeint(ystart0, ystart1, x1, x2, den, &kount, xp, yp, M, KMAX);
	
	printf("No of integration steps = %i\n",kount);
	
	
	
	//interpolate yking 
	for (k=0;k<M;k++) {
		x[k] = W0-rmin-sqrt(h*k); 
		
		
		if (x[k] > xp[0]) {
			yking[k][0] = (W0-x[k])*yp[0][0]/(W0 - xp[0]);
			yking[k][1] = yp[0][1];
			printf("+");
		} else {
			
			i = 0;
			
			do {
				if ((x[k]-xp[i])*(x[k]-xp[i+1]) <= 0.0) {
					yking[k][0] = yp[i][0] + (yp[i+1][0]-yp[i][0])*(x[k]-xp[i])/(xp[i+1]-xp[i]);
					yking[k][1] = yp[i][1] + (yp[i+1][1]-yp[i][1])*(x[k]-xp[i])/(xp[i+1]-xp[i]);
					goto jump;
				} else {
					i++;
				}
			} while (i<kount);
		}
		
	jump:
		
		
		if (i >= kount) {
			yking[k][0] = yp[kount][0];
			yking[k][1] = yp[kount][1];
		}
		
		
		//get mass as a function of radius 
		mass[k] = -2.0*pow(yking[k][0],1.5)/yking[k][1];
		
		
		//calculate total potential energy
		if (k == 0) {
			pot = -0.5*pow(mass[k],2)/sqrt(yking[k][0]); 
		} else {
			pot -=  0.5*(pow(mass[k],2) - pow(mass[k-1],2))/(0.5*(sqrt(yking[k][0]) + sqrt(yking[k-1][0])));
		}
	}	
	
	
	
	
	/******************************
	 * DETERMINE BASIC QUANTITIES *
	 ******************************/
	
	//Total mass
	totmas = -2.0*pow(yking[M-2][0],1.5)/yking[M-2][1];
	
	//Half-mass radius
	k=0;
	
	do {
		k++;
	} while (mass[k] < 0.5*totmas);
	
	zh = yking[k][0] - (yking[k][0]-yking[k-1][0])*(mass[k]-0.5*totmas)/(mass[k]-mass[k-1]);
	*rh = sqrt(zh);
	
	//Virial radius and King radius
	*rvir = -pow(totmas,2)/(2.0*pot);
	*rking = sqrt(yking[M-2][0]);
	
	//Central velocity dispersion**2 (3-dimensional) 
	ve = sqrt(2.0*W0);
	cg1 = (-pow(ve,3)*exp(-pow(ve,2)/2.0) - 3.0*ve*exp(-pow(ve,2)/2.0) + 3.0/2.0*sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0) - pow(ve,5)*exp(-pow(ve,2)/2.0)/5.0)/(-ve*exp(-pow(ve,2)/2.0) + sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0)/2.0 - pow(ve,3)*exp(-pow(ve,2)/2.0)/3.0);
	
	printf("\nTheoretical values:\n");
	printf("Total mass (King units) = %g\n", totmas);
	printf("Viriral radius (King units) = %g\n", *rvir);
	printf("Half-mass radius (King units) = %g\t(Nbody units) = %g\n", *rh, *rh/ *rvir);
	printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", *rking, *rking/ *rvir);
	printf("Concentration = %g\n", log10(*rking));
	printf("Core radius (King units) = %g\t(Nbody units) = %g\n", 1.0, 1.0/ *rvir);
	printf("3d velocity dispersion**2: %g (central)\t %g (global)\n", cg1, -pot/totmas);
	
	
	
	/***************************
	 * GENERATE STAR POSITIONS *
	 ***************************/
	
	printf("\nGenerating Stars:\n");	
	
	if (1) {
		for (i=0;i<N;i++) {
		
			if ((i/1000)*1000 == i) printf("Generating orbits #%i\n", i);
		
			fmass = drand48()*mass[M-1];
		
			if (fmass < mass[0]) {
				r2 = pow(fmass/mass[0],2.0/3.0)*yking[0][0];
			} else {
				j = 0;
				do {
					j++;
					if (j>M) printf("WARNING: failing iteration\n");
				} while (mass[j] <= fmass);
			
				r2 = yking[j-1][0] + (fmass-mass[j-1])*(yking[j][0] - yking[j-1][0])/(mass[j]-mass[j-1]);
			}
		
		
			r = sqrt(r2);
			costh = 2.0*drand48()-1.0;
			sinth = sqrt(1.0-pow(costh,2));
			phi = 2.0*PI*drand48();
			xstar = r*sinth*cos(phi);
			ystar = r*sinth*sin(phi);
			zstar = r*costh;

			
			/****************
			 * CHOOSE SPEED *
			 ****************/
						
			if (r < *rking) {
				if (r < sqrt(yking[0][0])) {
					w1 = x[0];
					w = W0 - (r2/yking[0][0])*(W0 - w1);
				} else {
					j = 0;
					do {
						j++;
					} while (r > sqrt(yking[j][0]));
					wj1 = x[j-1];
					wj = x[j];
					w = wj1 + (r2-yking[j-1][0])*(wj-wj1)/(yking[j][0]-yking[j-1][0]);
				}
			} else {
				printf("radius too big\n");
			}
		
			vmax = sqrt(2.0*w);
			do {
				speed = vmax*drand48();
				fstar = pow(speed,2)*(exp(-0.5*pow(speed,2))-exp(-1.0*w));
			} while (fstar < 2.0*drand48()/exp(1.0));
		
			costh = 2.0*drand48()-1.0;
			phi = 2.0*PI*drand48();
			sinth = sqrt(1.0-pow(costh,2));
			ustar = speed*sinth*cos(phi);
			vstar = speed*sinth*sin(phi);
			wstar = speed*costh;
			mstar = star[i][0];
		
			//printf("i: %i\tr=%g\tm=%.5f\tx=%.5f y=%.5f z=%.5f\tvx=%.5f vy=%.5f vz=%.5f\n",i,sqrt(r2),mstar,xstar,ystar,zstar,ustar,vstar,wstar);
		
		
			coord[i][0] = xstar;
			coord[i][1] = ystar;
			coord[i][2] = zstar;
			coord[i][3] = ustar;
			coord[i][4] = vstar;
			coord[i][5] = wstar;
		
		}

		for (i=0;i<N;i++) {
			for (j=0;j<3;j++) {
				star[i][j+1] = 1.0*coord[i][j];
				star[i][j+4] = 1.0*coord[i][j+3];
			}
		}

	}
	
		
	for (j=0;j<KMAX;j++) free (yp[j]);
	free(yp);

	for (j=0;j<N;j++) free (coord[j]);
	free(coord);
	
	return 0;
	

}

double densty(double z){
	double den = -sqrt(z)*(z+1.5)+0.75*sqrt(PI)*exp(z)*erf(sqrt(z));
	return den;
}

int odeint(double ystart0, double ystart1, double x1, double x2, double den, int *kount, double *xp, double **yp, int M, int KMAX) {
	
	double HMIN = 0.0;			//Minimum step size
	double H1 = 0.0001;			//Size of first step
	int MAXSTP = 100000;		//Maximum number of steps for integration
	double TINY = 1.0E-30;		//To avoid certain numbers get zero
	double DXSAV = 0.0001;		//Output step size for integration
	double TOL = 1.0E-12;		//Tolerance of integration
	
	
	double x;
	double h;
	int i,j;
	double y[2];
	double hdid, hnext;
	double xsav;
	double dydx[2], yscal[2];
	
	x = x1;   //King's W parameter
	if (x2-x1 >= 0.0) {
		h = sqrt(pow(H1,2));
	} else {
		h = - sqrt(pow(H1,2));
	}  //step size
	
	y[0] = ystart0;  //z^2
	y[1] = ystart1;  //2*z*dz/dW  where z is scaled radius	
	
	xsav = x-DXSAV*2.0;
	
	for (i=0;i<MAXSTP;i++) {        //limit integration to MAXSTP steps
		
		derivs(x,y,dydx,den); //find derivative
		
		for (j=0;j<2;j++) {
			yscal[j] = sqrt(pow(y[j],2))+sqrt(h*dydx[j])+TINY;  //advance y1 and y2
		}
		
		if (sqrt(pow(x-xsav,2)) > sqrt(pow(DXSAV,2))) {
			if (*kount < KMAX-1) {
				xp[*kount] = x;
				for (j=0;j<2;j++) {
					yp[*kount][j] = y[j];
				}
				*kount = *kount + 1;
				xsav = x;
			}
		}  //store x, y1 and y2 if the difference in x is smaller as the desired output step size DXSAV
		
		if (((x+h-x2)*(x+h-x1)) > 0.0) h = x2-x;
		
		rkqc(y,dydx,&x,&h,den,yscal, &hdid, &hnext, TOL);	//do a Runge-Kutta step
		
		if ((x-x2)*(x2-x1) >= 0.0) {
			ystart0 = y[0];
			ystart1 = y[1];
			
			xp[*kount] = x;
			for (j=0;j<2;j++) {
				yp[*kount][j] = y[j];
			}
			return 0;	
			*kount = *kount +1;
		}
		
		if (sqrt(pow(hnext,2)) < HMIN) {
			printf("Stepsize smaller than minimum.\n");
			return 0;
		}
		
		h = hnext;
	} 
	
	return 0;
}

int derivs(double x, double *y, double *dydx, double den){
	
	double rhox;
	
	if (x >= 0.0) {
		rhox =-sqrt(x)*(x+1.5)+0.75*sqrt(PI)*exp(x)*erf(sqrt(x));
	} else {
		rhox = 0.0;
	}
	
	dydx[0]= y[1];
	dydx[1] = 0.25*pow(y[1],2)*(6.0+9.0*y[1]*rhox/den)/y[0];	
	
	return 0;
}

int rkqc(double *y,double *dydx, double *x, double *h, double den, double *yscal, double *hdid, double *hnext, double TOL){
	
	double safety = 0.9;
	double fcor = 0.0666666667;
	double errcon = 6.0E-4;
    double pgrow = -0.20;
	double pshrnk = -0.25;
	double xsav;
	int i;
	double ysav[2],  dysav[2], ytemp[2];
	double errmax;
	double hh;
	
	xsav = *x;
	
	for (i=0;i<2;i++) {
		ysav[i] = y[i];
		dysav[i] = dydx[i];
	}
	
	do {
		hh = 0.5**h;
		rk4(xsav, ysav, dysav, hh, ytemp, den);
		
		*x = xsav + hh;
		derivs(*x,ytemp,dydx,den); //find derivative
		rk4(*x, ytemp, dydx, hh, y, den);
		
		*x = xsav + *h;
		if (*x  == xsav) {
			printf("ERROR: Stepsize not significant in RKQC.\n");
			return 0;
		}
		rk4(xsav, ysav, dysav, *h, ytemp, den);
		
		errmax = 0.0;
		for (i=0;i<2;i++) {
			ytemp[i] = y[i] - ytemp[i];
			errmax = max(errmax, sqrt(pow(ytemp[i]/yscal[i],2)));
		}
		errmax /= TOL;
		if (errmax > 1.0) *h = safety**h*(pow(errmax,pshrnk)); //if integration error is too large, decrease h
		
	} while (errmax > 1.0);
	
	
	*hdid = *h;
	if (errmax > errcon) {
		*hnext = safety**h*(pow(errmax,pgrow));//increase step size for next step
	} else {
		*hnext = 4.0**h;//if integration error is very small increase step size significantly for next step
	}
	
	for (i=0;i<2;i++) {
		y[i] += ytemp[i]*fcor;
	}
	
	return 0;
	
}

int rk4(double x, double *y, double *dydx, double h, double *yout, double den){
	double hh, h6, xh;
	double yt[2], dyt[2], dym[2];
	int i;
	
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<2;i++) {
		yt[i] = y[i] + hh*dydx[i];
	}
	
	derivs(xh,yt,dyt,den); //find derivative
	for (i=0;i<2;i++) {
		yt[i] = y[i] + hh*dyt[i];
	}
	
	derivs(xh,yt,dym,den); //find derivative
	for (i=0;i<2;i++) {
		yt[i] = y[i] + h*dym[i];
		dym[i] += dyt[i];
	}
	
	derivs(x+h,yt,dyt,den); //find derivative
	for (i=0;i<2;i++) {
		yout[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}
	
	return 0;
}

void quich(int start, int stop) {
}

/*
void quick(int start, int stop) {
	int i, j;
	double temp, median;
	
	median = 0.5 * (star1[start].mass + star1[stop].mass);
	
	i = start;
	j = stop;
	while (i <= j)
	{
		while ((i <= stop) && (star1[i].mass > median)) i++;
		while ((j >= start) && (star1[j].mass < median)) j--;
		if (j > i)
		{
			SWAP(i,j);
			i++;
			j--;
		}
		else if (i == j)
		{
			i++;
			j--;
		};
	};
	
	if (start + 1 < j) quick(start, j);
	else if ((start < j) && (star1[start].mass < star1[j].mass)) SWAP(start, j);
	
	if (i + 1 < stop) quick(i, stop);
	else if ((i < stop) && (star1[i].mass < star1[stop].mass)) SWAP(i, stop)
}
*/
