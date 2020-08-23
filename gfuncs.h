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

/* Time as a function of scale factor in seconds */
// integrand
double timeofa_int(double ap, void * params) {
  double integrand = 1./(ap*hubble(om,orad,ap));
  return integrand;
}

double timeofa(double a){
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &timeofa_int;
	gsl_integration_qags (&F, ai, a, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
	return result;
}


/* scale factor as a function of time */
// fill array to then spline
void aoftime(arrays_T xxyy){
		for(int i=0; i<1001; i++){
			double a = ai*exp(i*log(0.9999/ai)/(1000.));
			(*xxyy).xx[i] = timeofa(a);
			(*xxyy).yy[i] = a;
		}
		// set final time
		(*xxyy).xx[1001] = 1./hubble(om,orad,1.);
		(*xxyy).yy[1001] = 1.;
}


/* Critical density in kg/m^3 */
double rhoc(double a){
	return 3.*pow(hubble(om, orad, a),2.)/(8.*M_PI*gnewton);
}

/* Omega_x for LCDM */
double omegalcdm(double ox, double a){
	return  ox/pow(a,3.) * pow(h0/hubble(om,orad,a),2.);
}


/* Log[Mass in kg] as a function of scale factor  -- step function for now*/
// double bhmasslg(double lgm, double a){
// 	double logdect = 1./3. * (17.803 + log10(timeofa(a)));
// 	if (lgm>logdect) {
// 		return lgm;
// 	}
// 	else{
// 		return -100;
// 	}
// }



/* Log10[Mass] as a function of scale factor  --  not a step function*/
// M = 10^lgm ,   a is the scale factor , rem= false for no planck mass remnants
double bhmasslg(double lgm, double a, bool rem){
	// If M_i^3 <= 3 K M_p^4 t  then the black hole should have  0 mass ---> Log10[M_i] <= 1/3 (Log10[ 3 K M_p^4] + Log10[t]) = 1/3 ( 17.803 + Log10[t])
	double logdect = 1./3. * (17.803 + log10(timeofa(a)));
	if (lgm>logdect) {
		return log10(pow(pow(10, lgm*3) - decfac*timeofa(a), 1.0/3));
	}
	else{
		if (!rem) {
			return -100;
		}
		else{
			return lgmp;
		}
	}
}

// calculate volume fration given number density
double epsilon(double lgni, double lgmass){
	double radius = 2.*pow(10.,lgmass)/pow(10.,lgmp*2.) /kgtgev;
	return pow(10.,lgni) * 4. * M_PI * pow(radius,3) / 3.  * pow(dmtgev,3);
}
