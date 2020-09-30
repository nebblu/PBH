#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <math.h>       /* pow */
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

// Common functions to all  .cpps


/* LCDM Hubble in 1/s */
double hubble(double om, double orad, double a){
	double a3 = pow(a,3);
	double a4 = a*a3;
	double ol = 1.-om-orad;
	return h0 * sqrt(om/a3+ orad/a4 + ol);
}


// Reheating time as a function of temperature --- Trh should be given in Kelvin
double timeofrh(double Trh){
	return sqrt(5./pow(M_PI,3)/gstar) * pow(10,lgmp)/pow(Trh*keltgev,2) *kgtgev / gevtds ;
}

/* Critical density in kg/m^3 */
double rhoc(double a){
	return 3.*pow(hubble(om, orad, a),2.)/(8.*M_PI*gnewton);
}

/* Omega_x for LCDM */
double omegalcdm(double ox, double a){
	return  ox/pow(a,3.) * pow(h0/hubble(om,orad,a),2.);
}


/* Log10 of Mass spectrum in units 1/kg */
// hi mass - Eq.2.5 of draft
//  Psi = 10^psihilg   as a function of Log10[M] = lgm and omega_pbh
double psihilg(double lgm, double opbh){
	double fpbh = opbh/oc;
	return -12.2764 + log10(fpbh) - 3./2. * lgm;
}

// low mass - Eq.2.2 of draft
//  Psi = 10^psilowlg  as a function of Log10[M] = lgm and peak mass
double psilowlg(double lgm, double peakm){
	double aofmf = (peakm+2.723)/(-0.6576);
	double expt = 2.85*(lgm-peakm);
	return 2.85*(aofmf + lgm) - pow(10,expt)/2.3026;
}

// full spectrum - polychromatic
// Psi = 10^psibroadlg as a function of Log10[M] = lgm and peak mass
double psibroadlg(double lgm, double peakm){
	//  if mass is less than planck mass, set to 0
	if(lgmp>lgm){
		 return -100.;
		}
	// if above peak mass, set to hi spectrum
	else if(lgm>peakm){
		return psihilg(lgm,oc);
		}
	//  if below peak mass, set to low spectrum
	else{
		return psilowlg(lgm,peakm);
	}
}



// calculate volume fration given number density
double epsilon(double lgni, double lgmass){
	double radius = 2. * gnewton * pow(10,lgmass)/pow(sofl,2);
	printf("%e \n", radius);
	return pow(10.,lgni) * 4. * M_PI * pow(radius,3) / 3.;
}

/* Entropy density */
// Temperature given in Kelvin
double entropy(double Trh,  double  a, double arh){
	return 2.*pow(M_PI,2)/45 * gstar * pow(Trh*keltgev*arh/a,3.);
}

// absolute maximum theoretical number density
// Some notes: A bh of mass 10^-5.96 kg has a volume of 1.76 x 10^-98 m^3  and so we can have number densities mathematically up to 10^98 x 0.74 where 0.74 is the close packing ratio (max ratio of cube to volume of packed spheres)
double maxn0(double  lgmass){
	double rsch = 2.*gnewton*pow(10.,lgmass)/pow(sofl,2); // schwarzschild radius in m
	double volume =4./3. * M_PI * pow(rsch,3);
	return 1./volume * 0.74;
}
