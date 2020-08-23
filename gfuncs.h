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
	if(lgmp - lgm>0){
		 return -100.;
		}
	// if above peak mass, set to hi spectrum
	else if(lgm-peakm>0){
		return psihilg(lgm,oc);
		}
	//  if below peak mass, set to low spectrum
	else{
		return psilowlg(lgm,peakm);
	}
}


/* Calculate normalisation for psibroad */
// normalisation integrand --- see equation 3.4 for example
// Params are peak mass, time of reheating (unused) and scale factor
// Normalisation = 1/ Integrate[ 10^psilg dM ] =   1/ Integrate[ 10^(psilg + lgm) x Ln[10] dlgm ]

double normlg_int(double lgm, void * params) {
	myparam_type pars = *(myparam_type *)(params);
	double peakm = pars.peakm;
	double trh = pars.trh;
	double a = pars.aval;
	bool rem = pars.rem;
	double lgbhm = bhmasslg(lgm,a, rem);
	double expt = lgm + psibroadlg(lgbhm,peakm);
	double integrand =pow(10.,expt);
  return integrand;
}

// Normalisation = 10^normlg
double normlg(double peakm, double a, bool rem){
	struct myparam_type pars = {peakm, 1., a, rem};

	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &normlg_int;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);

  // Log10[1/Ln[10]] = -0.362216
	return -0.362216 - log10(result);
}



// calculate volume fration given number density
double epsilon(double lgni, double lgmass){
	double radius = 2.*pow(10.,lgmass)/pow(10.,lgmp*2.) /kgtgev;
	return pow(10.,lgni) * 4. * M_PI * pow(radius,3) / 3.  * pow(dmtgev,3);
}

/* Entropy density */
// Temperature given in Kelvin
double entropy(double Trh,  double  a){
	return 2.*pow(M_PI,2)/45 * gstar * pow(Trh*keltgev*ai/a,3.);
}


// absolute maximum theoretical number density
// Some notes: A bh of mass 10^-5.96 kg has a volume of 1.76 x 10^-98 m^3  and so we can have number densities mathematically up to 10^98 x 0.74 where 0.74 is the close packing ratio (max ratio of cube to volume of packed spheres)
double maxn0(double  lgmass){
	double rsch = 2.*gnewton*pow(10.,lgmass)/pow(sofl,2); // schwarzschild radius in m
	double volume =4./3. * M_PI * pow(rsch,3);
	return 1./volume * 0.74;
}
