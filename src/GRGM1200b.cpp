/*! \mainpage EGM2008.c Calculates gravity according to the EGM 2008 Spherical Harmonic Model
*  AUTHORS:          Austin Probe (abprobe88@gmail.com) and Brent Macomber (brentmacomber@gmail.com)
*  DATE WRITTEN:     October 2014
*  LAST MODIFIED:    May 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Computes Spherical Harmonic gravity for Orbit propagation using MCPI
*
*
*/


/*!
@file EGM2008.c
*/

/*! \file EGM2008.c EGM2008 Function File
* This file contains the fuctions for the EGM2008 Gravity Spherical
* harmonic libriary
*/
#include <math.h>


#include "GRGM1200b.h"
#include "GRGM1200b.cc"
#include "const.h"
#include "eci2ecef.h"
// #include "EphemerisRotation.h"

// Declare Needed Variables
/*!
* This is the allocation for the maximum size associated
* Legendre polynomical matrix
*/
double P_GRGM1200b[(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};
/*!
* This is the allocation for the maximum size associated
* Legendre polynomical matrix scale factor
*/
double scaleFactor_GRGM1200b[(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};
/*!
* \brief Gravity Evaluation
* This is the function that evaluates the spherical harmonic
* series to provide acceleration
*
* \param[in] p 3 element position in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] Gxyz Gravitational Acceloration Output
*/
void GRGM1200b( double* p, double* Gxyz, int DEG)
{

	double r             = {0.0};
	double phic          = {0.0};
	double lambda        = {0.0};
	double slambda       = {0.0};
	double clambda       = {0.0};
	double x = {0.0};
	double y = {0.0};
	double z = {0.0};
	double smlambda[Max_Degree_GRGM1200b+1] = {0.0};
	double cmlambda[Max_Degree_GRGM1200b+1] = {0.0};

	// determine radius of this thread's node
	x = p[0];
	y = p[1];
	z = p[2];

	// Compute geocentric radius
	r = pow( x*x + y*y + z*z , 0.5 );
	// Compute geocentric latitude
	phic  = asin( z / r );
	// Compute lambda
	lambda  = atan2( y, x );
	while (lambda<0)
	lambda = lambda+2*C_PI;
	while (lambda>=2*C_PI)
	lambda = lambda-2*C_PI;


	slambda = sin(lambda);
	clambda = cos(lambda);
	smlambda[0] = 0.0;
	cmlambda[0] = 1.0;
	smlambda[1] = slambda;
	cmlambda[1] = clambda;

	for(int m=2;m<DEG+1;m++){
		smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
		cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
	}

	double P_priv[(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};
	double scaleFactor_priv[(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};
	loc_gravLegendre_GRGM1200b( phic, scaleFactor_priv, P_priv , DEG);

	loc_gravityPCPF_GRGM1200b( p, P_priv, DEG, smlambda, cmlambda, r, scaleFactor_priv, Gxyz );

return;
}

/*!
* \brief Gravity Potential Evaluationc++ co
* This is the function that evaluates the spherical harmonic
* serries to provide gravitational potential
*
* \param[in] p 3 element position in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] Pot Gravitational Potential Output
*/
void GRGM1200bPot( double* p, double* Pot, int DEG)
{
	// determine radius of this thread's node
	double r             = {0.0};
	double phic          = {0.0};
	double lambda        = {0.0};
	double slambda       = {0.0};
	double clambda       = {0.0};
	double smlambda[Max_Degree_GRGM1200b+1] = {0.0};
	double cmlambda[Max_Degree_GRGM1200b+1] = {0.0};

	double x = p[0];
	double y = p[1];
	double z = p[2];
	int m;
	double P [(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};
	double scaleFactor_GRGM1200b [(Max_Degree_GRGM1200b+3)*(Max_Degree_GRGM1200b+3)] = {0.0};

	// Compute geocentric radius
	r = pow( x*x + y*y + z*z , 0.5 );
	// Compute geocentric latitude
	phic  = asin( z / r );
	// Compute lambda
	lambda  = atan2( y, x );
	while (lambda<0){
		lambda = lambda+2*C_PI;
	}
	while (lambda>=2*C_PI){
		lambda = lambda-2*C_PI;
	}

	slambda = sin(lambda);
	clambda = cos(lambda);
	smlambda[0] = 0.0;
	cmlambda[0] = 1.0;
	smlambda[1] = slambda;
	cmlambda[1] = clambda;


	for(m=2;m<DEG+1;m++){
		smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
		cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
	}
	// Compute normalized associated legendre polynomials

	loc_gravLegendre_GRGM1200b( phic, scaleFactor_GRGM1200b, P , DEG);

	loc_gravityPot_GRGM1200b( p, P, DEG, smlambda, cmlambda, r, scaleFactor_GRGM1200b, Pot );


}

/*!
* \brief Legendre Polyniomial Evaluation
* This is the function that computes the normalized associated
* legendre polynomials based on geocentric latitude
*
* \param[in] phi Geocentric latitude
* \param[in] Deg Degree and order of the serries to be used
* \param[out] P associated Legendre polynomial matrix
* \param[out] scaleFactor_GRGM1200b Legendre scale factor
*/
void loc_gravLegendre_GRGM1200b( double phi, double* scaleFactor_GRGM1200b, double* P, int DEG )
{

	int k, p;
	double cphi = {0.0};
	double sphi = {0.0};


	cphi = cos(0.5*C_PI - phi);
	sphi = sin(0.5*C_PI - phi);
	// Seeds for recursion formula
	P[IDX2F(1,1,Max_Degree_GRGM1200b+3)] = 1.0;            // n = 0, m = 0;
	scaleFactor_GRGM1200b[IDX2F(1,1, Max_Degree_GRGM1200b+3)] = 0.0;
	P[IDX2F(2,1, Max_Degree_GRGM1200b+3)] = sqrt(3.0)*cphi ; // n = 1, m = 0;
	scaleFactor_GRGM1200b[IDX2F(2,1,Max_Degree_GRGM1200b+3)]  = 1.0;
	P[IDX2F(2,2,Max_Degree_GRGM1200b+3)] = sqrt(3.0)*sphi; // n = 1, m = 1;
	scaleFactor_GRGM1200b[IDX2F(2,2,Max_Degree_GRGM1200b+3)] = 0.0;

	// Old Method
	for (int nn = 2; nn <= DEG+2;nn++){
		double n = (double)nn;
		k = nn + 1;
		for(int mm=0; mm<=n;mm++) {
			double m = (double)mm;
			p = mm + 1;
			// Compute normalized associated legendre polynomials, P, via recursion relations
			// Scale Factor needed for normalization of dUdphi partial derivative
			if (n == m){
				P[IDX2F(k,k,Max_Degree_GRGM1200b+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree_GRGM1200b+3)];
				scaleFactor_GRGM1200b[IDX2F(k,k,Max_Degree_GRGM1200b+3)] = 0.0;
			}
			else if (m == 0){
				P[IDX2F(k,p,Max_Degree_GRGM1200b+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree_GRGM1200b+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree_GRGM1200b+3)] );
				scaleFactor_GRGM1200b[IDX2F(k,p,Max_Degree_GRGM1200b+3)] = sqrt( (n+1)*(n)/2);
			}
			else {
				P[IDX2F(k,p,Max_Degree_GRGM1200b+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree_GRGM1200b+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree_GRGM1200b+3)] );
				scaleFactor_GRGM1200b[IDX2F(k,p,Max_Degree_GRGM1200b+3)] = sqrt( (n+m+1)*(n-m));
			}
		}
	}
}

void loc_gravityPCPF_GRGM1200b( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor_GRGM1200b, double* Gxyz )
{

	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};

	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double dUdr    = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};



	double x = p[0];
	double y = p[1];
	double z = p[2];
	radu = r;
	//FIXME: 
	rRatio = C_Rmoon/radu;
	rRatio_n = rRatio;
	// summation of gravity in radial coordinates
	for (int n = 2; n <= DEG; n++) {
		k = n+1;
		rRatio_n = rRatio_n*rRatio;
		dUdrSumM      = 0.0;
		dUdphiSumM    = 0.0;
		dUdlambdaSumM = 0.0;
		for (int m = 0; m <= n; m++){
			j = m+1;
			dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree_GRGM1200b+3)] *(C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] + S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
			dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree_GRGM1200b+3)]*scaleFactor_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree_GRGM1200b+3)])*(C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] + S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
			dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree_GRGM1200b+3)]*(S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] - C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
		}
		dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
		dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
		dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
	}

	// gravity in spherical coordinates
	dUdr      = -GM_GRGM1200b/(radu*radu)*dUdrSumN ;
	dUdphi    =  GM_GRGM1200b/radu*dUdphiSumN ;
	dUdlambda =  GM_GRGM1200b/radu*dUdlambdaSumN ;

	//gravity in ECEF coordinates
	Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
	Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
	Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;

	// special case for poles
	/*
	atPole = abs(atan2(p[IDX2F(i+1,3,M)],sqrt(p[IDX2F(i+1,1,M)]^2 + p[IDX2F(i+1,2,M)]^2)))==C_PI/2;
	if any(atPole){
	gx(atPole) = 0;
	gy(atPole) = 0;
	gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
}
*/

}

void loc_gravityPot_GRGM1200b( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor_GRGM1200b, double* Pot )
{

	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};

	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double dUdr  = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};
	double USum = {0.0};


	double x = p[0];
	double y = p[1];
	double z = p[2];
	int n;
	int m;

	radu = r;
	rRatio = C_Rmoon/radu;
	rRatio_n = rRatio;
	// summation of gravity in radial coordinates

	for (n = 2; n <= DEG; n++) {
		k = n+1;
		rRatio_n = rRatio_n*rRatio;
		dUdrSumM      = 0.0;
		dUdphiSumM    = 0.0;
		dUdlambdaSumM = 0.0;
		for (m = 0; m <= n; m++){
			j = m+1;
			dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree_GRGM1200b+3)] *(C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] + S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
			dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree_GRGM1200b+3)]*scaleFactor_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree_GRGM1200b+3)])*(C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] + S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
			dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree_GRGM1200b+3)]*(S_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*cmlambda[j-1] - C_GRGM1200b[IDX2F(k,j,Max_Degree_GRGM1200b)]*smlambda[j-1]);
			// printf("UsumM: %e\n ",dUdrSumM);
		}
		dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n;
		dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
		dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
	}
	// gravity in spherical coordinates
	dUdr      =  GM_GRGM1200b/(radu)*dUdrSumN ;
	dUdphi    =  GM_GRGM1200b/radu*dUdphiSumN ;
	dUdlambda =  GM_GRGM1200b/radu*dUdlambdaSumN ;

	// gravity in potential
	// *Pot = MU * sqrt(dUdr*dUdr + dUdphi*dUdphi + dUdlambda*dUdlambda)/radu;
	*Pot =  sqrt(dUdr*dUdr);


}

void jacobiIntegral_GRGM1200b(double t, double* solN, double* H, int Deg, Orbit &orbit){
	double xI[3]    = {0.0};
	double vI[3]    = {0.0};
	double xF[3] = {0.0};
	double vF[3] = {0.0};

	xI[0] = solN[0];
	xI[1] = solN[1];
	xI[2] = solN[2];
	vI[0] = solN[3];
	vI[1] = solN[4];
	vI[2] = solN[5];

	// Convert from inertial frame to body fixed frame
	//InertialToBodyFixed(xI,vI,xF,vF,t,orbit);
	eci2ecef(t,xI,vI,xF,vF);
	double KE,PE,RotTerm;

	KE = 0.5*(vF[0]*vF[0] + vF[1]*vF[1] + vF[2]*vF[2]);		//Kinetic energy in inertial frame
	GRGM1200bPot(xF, &PE, Deg);								//Potential energy from body fixed position
	PE = -PE;
	RotTerm = 0.5*C_omega*C_omega*(xF[0]*xF[0] + xF[1]*xF[1]);

	*H  = PE + KE - RotTerm; // Hamiltonian

	// printf("KE: %e\tPE: %e\tRT: %e\tSum: %e\n ",KE,PE,RotTerm,*H);
	// getchar();

}
