#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "constant.h"

#include "binary_pairing.h"
#include "make_mass.h"
#include "randomz.h"

static double eigenevolution(double *m1, double *m2, double *P, double *e)
{
	double R_peri;
	double lambda=28,chi=0.75,rou,rou2;
	double e_birth=*e,e_fin,P_birth=*P,P_fin;
	double mtot_birth,mtot_fin;
	double q_birth,q_fin;
	double temp;
	temp = *m1;
	if (*m1 < *m2) {
		*m1 = *m2;
		*m2 = temp;
	}
	mtot_birth = *m1 + *m2;

	R_peri = (1-e_birth) * pow(P_birth*P_birth*((*m1)+(*m2)), (1.0/3.0));
	
	rou = pow((lambda*R_SUN/R_peri), chi);
	e_fin = exp(-rou + log(e_birth));
	*e = e_fin;

	q_birth = *m2 / *m1;
//	if ( q_birth>1.0 ) q_birth = 1.0/q_birth;

	rou2 = (rou<1.0) ? rou : 1.0;

	q_fin = q_birth + (1-q_birth) * rou2;

	*m2 = *m1 * q_fin;

	mtot_fin = *m1 + *m2;

	P_fin = P_birth * sqrt(mtot_birth/mtot_fin) * pow((1-e_birth)/(1-e_fin), 1.5);
	*P = P_fin;

	return mtot_fin;
}

static double orbital_kepler(double M_orbit, double e)
{
	double f,df,dx,xacc;
	double E_orbit=M_orbit;
	int i=0;
	int imax=50;
	xacc = 1E-6;

	for (i=0;i<imax;++i){
		f = E_orbit - e * sin(E_orbit) - M_orbit;
		df = 1 - e * cos(E_orbit);
		dx = f/df;
		E_orbit -= dx;
		if (fabs(dx)<xacc) break;
	}
//	do {
//		f = E_orbit - e * sin(E_orbit) - M_orbit;
//		df = 1 - e * cos(E_orbit);
//		dx = f/df;
//		E_orbit -= dx;
//		++i;
//		fprintf(stderr,"orbital_kepler: i=%d E_orbit=%lf dx=%lf\n",i,E_orbit,dx);
//	} while ( fabs(dx)>xacc && i<imax );
	return E_orbit;
}

void generate_binaries(struct star *star_x, int N_star, int nbin, int nsimf, double *ms, double *as, int pairing_type)
{
	int i,j;
	double a,P,e;
	double log_P,eta,delta,logPmin;
	double x;
	eta=2.5;
	delta=45;
	logPmin=1;
	double r_P,r_Q,v_P,v_Q,r_orbit[3],v_orbit[3],P_matrix[3],Q_matrix[3];
	double mu,sqee;
	double M_orbit,E_orbit,sin_E_orbit,cos_E_orbit;
	double i_orbit,cos_i,sin_i;
	double Omega_longitude_orbit,cos_Omega_longitude,sin_Omega_longitude;
	double omega_periapsis_orbit,cos_omega_periapsis,sin_omega_periapsis;

	double m1,m2,m_cm;
//	double rvir;

	binary_pairing(star_x,N_star,&nbin,nsimf,ms,as,pairing_type);

	for (i=0;i<nbin;++i) {

		m1 = star_x[2*i  ].m;
		m2 = star_x[2*i+1].m;
		do {
			x = randomz();
			log_P = logPmin + sqrt( delta * ( exp(2.0*x/eta) -1.0 ) );

			P = pow(10,log_P);	// days
			P /= 365.25;		// years
			a = pow( (m1+m2)*P*P,(1.0/3.0) );	// AU
		} while (a<10*R_SUN);
		a /= 206264.806247904;				// pc
//		a /= rvir;							// N-body units

		x = randomz();
		e = sqrt(x);		// thermal distribution f(e)=2e

		// Apply Kroupa (1995) eigenevolution
		m_cm = eigenevolution(&m1,&m2,&P,&e);
		// m_cm = m1+m2;
		star_x[2*i  ].m = m1;
		star_x[2*i+1].m = m2;

		// position & velocity in binary frame
		M_orbit = TWO_PI*randomz();
		E_orbit = orbital_kepler(M_orbit,e);
		sin_E_orbit = sin(E_orbit);
		cos_E_orbit = cos(E_orbit);
		sqee = sqrt(1-e*e);

		r_P = a * ( cos_E_orbit - e );
		r_Q = a * sqee * sin_E_orbit;

		mu = sqrt( (m1+m2)/a ) / ( 1-e*cos_E_orbit ); // mu*VSC (in [m/s])
		v_P = -mu * sin_E_orbit;
		v_Q = mu * sqee * cos_E_orbit;

		cos_i = 2.0*randomz()-1.0;
		i_orbit = acos(cos_i);
		sin_i = sin(i_orbit);
		Omega_longitude_orbit = 2.0*PI*randomz();
		sin_Omega_longitude = sin(Omega_longitude_orbit);
		cos_Omega_longitude = cos(Omega_longitude_orbit);
		omega_periapsis_orbit = 2.0*PI*randomz();
		sin_omega_periapsis = sin(omega_periapsis_orbit);
		cos_omega_periapsis = cos(omega_periapsis_orbit);

		P_matrix[0] = cos_Omega_longitude*cos_omega_periapsis - sin_Omega_longitude*sin_omega_periapsis*cos_i;
		P_matrix[1] = sin_Omega_longitude*cos_omega_periapsis + cos_Omega_longitude*sin_omega_periapsis*cos_i;
		P_matrix[2] = sin_omega_periapsis*sin_i;
		Q_matrix[0] = -cos_Omega_longitude*sin_omega_periapsis - sin_Omega_longitude*cos_omega_periapsis*cos_i;
		Q_matrix[1] = -sin_Omega_longitude*sin_omega_periapsis + cos_Omega_longitude*cos_omega_periapsis*cos_i;
		Q_matrix[2] = cos_omega_periapsis*sin_i;

		for (j=0;j<3;++j) {
			r_orbit[j] = r_P*P_matrix[j] + r_Q*Q_matrix[j];
			v_orbit[j] = v_P*P_matrix[j] + v_Q*Q_matrix[j];
		}

		star_x[2*i+1].x[0] =  m1/m_cm * r_orbit[0];
		star_x[2*i+1].x[1] =  m1/m_cm * r_orbit[1];
		star_x[2*i+1].x[2] =  m1/m_cm * r_orbit[2];
		star_x[2*i+1].x[3] =  m1/m_cm * v_orbit[0];
		star_x[2*i+1].x[4] =  m1/m_cm * v_orbit[1];
		star_x[2*i+1].x[5] =  m1/m_cm * v_orbit[2];

		star_x[2*i  ].x[0] = -m2/m_cm * r_orbit[0];
		star_x[2*i  ].x[1] = -m2/m_cm * r_orbit[1];
		star_x[2*i  ].x[2] = -m2/m_cm * r_orbit[2];
		star_x[2*i  ].x[3] = -m2/m_cm * v_orbit[0];
		star_x[2*i  ].x[4] = -m2/m_cm * v_orbit[1];
		star_x[2*i  ].x[5] = -m2/m_cm * v_orbit[2];

//		// position & velocity in cluster frame
//		star_x[2*i+1].x[0] += m1/m_cm * r_orbit[0];
//		star_x[2*i+1].y[1] += m1/m_cm * r_orbit[1];
//		star_x[2*i+1].z[2] += m1/m_cm * r_orbit[2];
//		star_x[2*i+1].v[3] += m1/m_cm * v_orbit[0];
//		star_x[2*i+1].v[4] += m1/m_cm * v_orbit[1];
//		star_x[2*i+1].v[5] += m1/m_cm * v_orbit[2];

//		star_x[2*i  ].x[0] -= m2/m_cm * r_orbit[0];
//		star_x[2*i  ].y[1] -= m2/m_cm * r_orbit[1];
//		star_x[2*i  ].z[2] -= m2/m_cm * r_orbit[2];
//		star_x[2*i  ].v[3] -= m2/m_cm * v_orbit[0];
//		star_x[2*i  ].v[4] -= m2/m_cm * v_orbit[1];
//		star_x[2*i  ].v[5] -= m2/m_cm * v_orbit[2];

		// update centre-mass for binaries
		star_x[N_star+i].m = m_cm;
	}

	return;
}
