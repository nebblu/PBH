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


// has an initial matter domination phase before reheating with matter density given by peak mass x number density
double hubblemd(double a, myparam_type2 pars){
	double peakm = pars.peakm;
	double trh = pars.trh;
	double lgn0 = pars.lgn0;

	double rhoconst =  8. * M_PI * gnewton * pow(10.,lgn0+peakm)/ 3.;
	double arh = pow(3. * sqrt(rhoconst) / 2. * trh + pow(ai,3./2.), 2./3. );


	if (a<arh) {
	double  hubble2 = rhoconst/pow(a,3.);
	return sqrt(hubble2);
	}
	else{
	double a3 = pow(a,3);
	double a4 = a*a3;
	double ol = 1.-om-orad;
	return h0 * sqrt(om/a3+ orad/a4 + ol);
	}
}


/* Time as a function of scale factor in seconds */
// integrand
double timeofa_int(double ap, void * params) {
	myparam_type2 pars = *(myparam_type2 *)(params);
  double integrand = 1./(ap*hubblemd(ap, pars));
  return integrand;
}



double timeofa(double a, myparam_type2 pars){
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.params = &pars;
	F.function = &timeofa_int;
	gsl_integration_qags (&F, ai, a, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
	return result;
}


/* scale factor as a function of time */
// fill array to then spline
void aoftime(arrays_T xxyy, myparam_type2 pars){
		for(int i=0; i<1001; i++){
			double a = ai*exp(i*log(0.9999/ai)/(1000.));
			(*xxyy).xx[i] = timeofa(a,  pars);
			(*xxyy).yy[i] = a;

		}
		// set final time
		(*xxyy).xx[1001] = 1./hubblemd(1., pars);
		(*xxyy).yy[1001] = 1.;
}


/* Log10[Mass] as a function of scale factor  --  not a step function*/
// M = 10^lgm ,   a is the scale factor , rem= false for no planck mass remnants
double bhmasslg(double lgm, void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double a = pars.aval;
	bool rem = pars.rem;
	// If M_i^3 <= 3 K M_p^4 t  then the black hole should have  0 mass ---> Log10[M_i] <= 1/3 (Log10[ 3 K M_p^4] + Log10[t]) = 1/3 ( 17.803 + Log10[t])
	double logdect = 1./3. * (17.803 + log10(timeofa(a, pars)));
	if (lgm>logdect) {
		return log10(pow(pow(10, lgm*3) - decfac*timeofa(a, pars), 1.0/3));
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



/* Calculate normalisation for psibroad */
// normalisation integrand --- see equation 3.4 for example
// Params are peak mass, time of reheating (unused) and scale factor
// Normalisation = 1/ Integrate[ 10^psilg dM ] =   1/ Integrate[ 10^(psilg + lgm) x Ln[10] dlgm ]

double normlg_int(double lgm, void * params) {
	myparam_type2 pars = *(myparam_type2 *)(params);
	double peakm = pars.peakm;
	double lgbhm = bhmasslg(lgm, params);
	double expt = lgm + psibroadlg(lgbhm,peakm);
	double integrand =pow(10.,expt);
  return integrand;
}

// Normalisation = 10^normlg
double normlg(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
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
