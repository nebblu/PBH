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


// functions particular to bh_clean.cpp

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


/* Log10[Mass] as a function of scale factor  --  not a step function*/
// M = 10^lgm ,   a is the scale factor , rem= false for no planck mass remnants
double bhmasslg(double lgm, double a, bool rem){
	// If M_i^3 <= 3 K M_p^4 t  then the black hole should have  0 mass ---> Log10[M_i] <= 1/3 (Log10[ 3 K M_p^4] + Log10[t]) = 1/3 ( 17.803 + Log10[t])
	double logdect = 1./3. * (17.803 + log10(timeofa(a)));
	if (lgm>logdect) {
		return log10(pow(pow(10., lgm*3.) - decfac*timeofa(a), 1.0/3.));
	}
	else{
		if (!rem) {
			return -100.;
		}
		else{
			return lgmp;
		}
	}
}


// Initial mass from final mass. Any bh that has decayed completely earlier than a will be mapped to decfac*timeofa(a)^(1/3)
// (There's no way of knowing with certainty what its initial mass is if its final mass is 0)
double bhmasslgi(double lgm, double a, bool rem){
  double logdect = 1./3. * (17.803 + log10(timeofa(a)));
  double mymasslg = log10(pow(10., lgm*3.) + decfac*timeofa(a))/3.;
  if(lgm<logdect ){
    return -100;
  }
  else if(mymasslg>lgmax){
		return -100;
    }
  else{
    return mymasslg;
  }
}




/* Calculate normalisation for psibroad */
// normalisation integrand --- see equation 3.4 for example
// Params are peak mass (as proxy for M_s)
// Normalisation = 1/ Integrate[ 10^psilg dM ] =   1/ Integrate[ 10^(psilg + lgm) x Ln[10] dlgm ]

double normlg_int(double lgm, void * params) {
	myparam_type pars = *(myparam_type *)(params);
	double peakm = pars.peakm;
	double expt = lgm + psibroadlg(lgm,peakm);
	double integrand =pow(10.,expt);
  return integrand;
}

// Normalisation = 10^normlg
double normlg(double peakm){
	struct myparam_type pars = {peakm, 1., ai, false};

	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &normlg_int;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 1e-10, 1e-10, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
	return  log10(1./(log(10.)*result));
}



/* Calculate normalisation for psi_mono */
// normalisation integrand --- see equation 3.4 for example
// Params are peak mass
// Normalisation = 1/ Integrate[ 10^psilg dM ] =   1/ Integrate[ 10^(psilg + lgm) x Ln[10] dlgm ]

double normlg_int_mono(double lgm, void * params) {
	myparam_type pars = *(myparam_type *)(params);
	double peakm = pars.peakm;
	double expt = lgm + psilowlg(lgm,peakm);
	double integrand =pow(10.,expt);
  return integrand;
}

// Normalisation = 10^normlg
double normlg_mono(double peakm){
	struct myparam_type pars = {peakm, 1., ai, false};

	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &normlg_int_mono;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 1e-10, 1e-10, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
	return  log10(1./(log(10.)*result));
}
