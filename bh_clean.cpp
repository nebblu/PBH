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

/* Conversion constants */
// megaparsec to km
const double mptkm = 3.086e19;
// kg to Gev
const double kgtgev = 5.62e26;
// Gev to kg
const double gevtkg = 1.78e-27;
//// 1/seconds to GeV
const double dstgev = 6.58e-25;
// Gev to 1/seconds
const double gevtds = 1.52e24;
// Kelvin to 1/s
const double keltds = 8.62e-14;
// Solar mass in kg
const double msol = 1.989e30;
// 1/m to GeV
const double dmtgev = 1.98e-16;
// Kelvin to Gev
const double keltgev = 8.62e-14;


/* LCDM constants and yield from Planck 2018 1807.06209 */
// Hubble constant in inverse seconds
const double h0 = 67.32/mptkm;
// Cold dark matter fraction from Planck 2018
const double oc = 0.265;
// Baryonic matter fraction from Planck 2018
const double ob = 0.0494;
// Total matter fraction from Planck 2018
const double om = oc+ob;
// Radiation fraction from Planck 2018
const double orad = 0.00006;
// Yield of B-L = n_B-L/ entropy density
const double yield = 1e-9;
// Hooks value for yield 1404.0113
const double yieldhook = 9e-11;

/* Scale factors */
const double ai = 1e-30; // end  of inflation  ---- need to adjust
const double arad = 0.000264; // dm-rad equality
const double acmb = 9e-4; // CMB
const double alam = 0.5; // DM - Lambda equality

/* Mass constants */
// Planck mass in kg
const double mplanck = 2.176e-8;
// log10 of planck mass (in kg)
const double lgmp = -7.66234;
// log10 of maximum mass used in integrations in kg
const double lgmax = 36.;
// horizon mass in kg at time of radiation-matter equality --- 2001.04371
const double meq = 2.8e17;
// Horizon mass when shortest wavelength reenters horizon - choice
const double ms = 1e-7;
// log of 5 times ms
const double lgmf = -6.301;
// proton mass in kg
const double protonm = 1.672621e-27;


/* Decay constants */
// Sum of degrees of freedom of standard model -- see 1404.0113 or our draft
const double gstar = 106.75;

// sum of qi^2 gi
const double sumq2g = 13.;

// 3 x gstar / (30720 pi) * Mplanck^4 in kg^4
const double cc = 7.4397e-34;

// 3 x gstar / (30720 pi) * Mplanck^4 in kg^3/s
const double decfac = 6.35e17; // = cc*gevtds/gevtkg

/* Newton's G in m^3 /kg / s^2 */
const double gnewton = 6.67e-11;

/* Speeed of light in m/s */
const double sofl = 2.99792e8;

// conversion factor from dM to dlog10[M]
const double intf = 2.30259;

/* All parameters should be given in SI base units - conversions are handled in the functions*/

// structure to hold relevant parameters
struct myparam_type2{
	double lambda; // coupling constant
 	double lgn0;  // log10 of final total number density of bh in 1/m^3
	double peakm; // log10  of mass spectrum peak  in kg
	double Trh; // reheating Temperature in Kelvin
	double trh; // time of reheatting in s
	double aval; // scale factor
	bool rem; // remnants or not
	gsl_spline *spline;
	gsl_interp_accel *acc;
};

// Reduced structure for epoch/approximate calculations
struct myparam_type{
	double peakm; // mass spectrum peak
	double trh; // time of reheating
	double aval; // scale factor
	bool rem; // remnants or not
};

// structure for creating a(t) object
typedef struct arrays{
  int count;
  double xx[1002];
  double yy[1002]; } *arrays_T;


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


//Final, normalised log10 of broad spectrum

double psifinallg(double lgm, double peakm, double a,  bool rem){
		//  Log10[M(t)]
		double lgmoft = bhmasslg(lgm,a,rem);
		double expt = psibroadlg(lgmoft, peakm) + normlg(peakm, a, rem) ;
		return expt;
}


// Avg mass x 10^lgn0 - multiplied by  10^lgn0 to reduce magnitude for numerical stability
double avgm_int(double lgm, void * params){
	myparam_type pars = *(myparam_type *)(params);
	double peakm = pars.peakm;
	double a = pars.aval;
	bool rem = pars.rem;
	double expt = 2.*lgm +  psifinallg(lgm, peakm, a, rem);
	return pow(10.,expt);
}



// average mass (see avgm_int function) in kg
double avgm(double peakm, double a, bool rem){
	struct myparam_type pars = {peakm, 1., a, rem};
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &avgm_int;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return intf * result ;
}


/* Critical density in kg/m^3 */
double rhoc(double a){
	return 3.*pow(hubble(om, orad, a),2.)/(8.*M_PI*gnewton);
}

/* Omega_pbh */
// lgn0 = Log10[n0] where n0 is the number density of bh today
double omegapbh(double lgn0, double peakm, double a, bool rem){
	double avgmass = avgm(peakm, a, rem);
	return  pow(10.,lgn0) * avgmass / (rhoc(a) * pow(a,3.));
}

/* Omega_x for LCDM */
double omegalcdm(double ox, double a){
	return  ox/pow(a,3.) * pow(h0/hubble(om,orad,a),2.);
}


/* Epoch treatment - probably won't present these results in final draft */
/* We consider 3 epochs : ai to arad , arad to alam,  alam to 1 */
/*  The B-L number density integrand */
double nbl_int_ep(double lgm, void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double lambda = pars.lambda;
	double peakm = pars.peakm;
	double trh = pars.aval;
	double a = pars.aval;
	bool rem = pars.rem;
	gsl_spline *myaoft = pars.spline;
	gsl_interp_accel *acc = pars.acc;

	// determine the decay time
	double dect = pow(10., (3.*lgm-17.803));
	// so spline doesn't run into issues, set deca = 1 if the decay time exceeds lifetime of the universe....
	double deca;
	if (dect>=1./hubble(om,orad,1)) {
		deca = 1.;
	}
	else if(dect<=trh){
			return 0.; // no chemical potential before reheatingg time so no baryons from BH that have already decayed
		}
	else{
		deca = gsl_spline_eval (myaoft, dect, acc);
	}

	// reheating scale factor
	double arh = gsl_spline_eval (myaoft, trh, acc);

	// chem  pot and q prefactor
	 double prefac = -9.*sumq2g/(96.*M_PI *pow(10.,2.*lgmp));

	 /*Add various charge contributions */
 	double lgqr, lgqm;
 	double exptq1,exptq2,exptq3,exptq4;
 	double w0r, w0m, hubr, hubm, qfac;

	if (a<=arad) {
		double af = gsl_min(deca,a);
		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
		w0r = (1.-1e-3)/3.;
		double eost = (1.+w0r)*(1.-3.*w0r);
		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
		double a2 = pow(af,2);
		double a3 = af*a2;
		double ai2 = pow(arh,2);
		double ai3 = arh*a2;
		double ol = 1.-om-orad;
		hubr = 1./6. * pow(h0,2) * ( 3.*om * ((a2-ai2)/(a2*ai2))  + 2.*orad * (a3-ai3/(a3*ai3)) + 6.*ol*(af-arh));
		// //printf("%e \n", lghubr);
		// double lgmur = 15.0783  +  lghubr + lglam; // 15.0783 = log10[9/Mp^2 * 19 / 96/pi]
		qfac = lambda*prefac*eost*hubr; //pow(10.,lgmur);
	}
	// same thing but if we run into matter domination
	else{
		/* Radiation terms */
		double af = gsl_min(arad,deca);

		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
		w0r = (1.-1e-3)/3.;
		double eost = (1.+w0r)*(1.-3.*w0r);
		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
		double a2 = pow(af,2);
		double a3 = af*a2;
		double ai2 = pow(arh,2);
		double ai3 = arh*a2;
		double ol = 1.-om-orad;
		hubr =eost*(1./6. * pow(h0,2) * ( 3.*om * ((a2-ai2)/(a2*ai2))  + 2.*orad * (a3-ai3/(a3*ai3)) + 6.*ol*(af-arh)));

		/*Matter terms */
		double afm;
		// check if bh has decayed before mat-de equality
		afm = gsl_min(a,alam);
		// check if bh has decayed before final time
		afm = gsl_min(afm,deca);
		// check if black hole has decayed before arad
		afm = gsl_max(afm,arad);
		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
		w0m = 0.;
		double eostm = (1.+w0m)*(1.-3.*w0m);
		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
		double a2m = pow(af,2);
		double a3m = af*a2m;
		double ai2m = pow(arad,2);
		double ai3m = arad*ai2m;
		hubm = eostm*(1./6. * pow(h0,2) * ( 3.*om * ((a2m-ai2m)/(a2m*ai2m))  + 2.*orad * (a3m-ai3m/(a3m*ai3m)) + 6.*ol*(afm-arad)));

		qfac = lambda*prefac*(hubr + hubm);
	}

	// psi term
	double exptpsi = psifinallg(lgm, peakm, ai, rem);

	return  intf * pow(10.,exptpsi) * qfac;
}


double nbl_ep(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double peakm = pars.peakm;
	double lgn0 = pars.lgn0;
	double trh = pars.trh;
	double a = pars.aval;
	bool rem = pars.rem;

	// can't create asymmetry before reheating
	double time = timeofa(a);
	if (time<=trh) {
		return  0.;
	}

	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &nbl_int_ep;
	F.params = &pars;

	double mavgi = avgm(peakm, ai, rem);
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return  pow(10.,lgn0) * mavgi * result / pow(a,3);
}


/* Omega_baryon for epoch treatment */
double omegab_ep(void * params){
		myparam_type2 pars = *(myparam_type2 *)(params);
		double a = pars.aval;
		double rhob = protonm * nbl_ep(params)/2.;
		return rhob /rhoc(a);
}



/* LCDM background treatment */
// we assume Hubble  = LCDM with ob,ocdm,orad
/*  The B-L number density integrand */
double nbl_int_lcdm(double lgm, void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double lambda = pars.lambda;
	double peakm = pars.peakm;
	double trh = pars.trh;
	double a = pars.aval;
	bool rem = pars.rem;
	gsl_spline *myaoft = pars.spline;
	gsl_interp_accel *acc = pars.acc;

	// determine the decay time
	double dect = pow(10., (3.*lgm-17.803));

	// time at which to stop integrating the charge, RHS of eq.14 of 1404.0113
	double deca;
	// so spline doesn't run into issues, set deca = 1 if the decay time exceeds lifetime of the (LCDM) universe....
	if (dect>=1./hubble(om,orad,1)) {
		deca = 1.;
	}
	// if the decay time is less than the time of reheating we should return 0.
	// no chemical potential before reheatingg time so no baryons from BH that have already decayed
	else if(dect<=trh){
		return 0.;
	}
	// otherwise we find the scale factor at which to stop integrating
	else{
		deca = gsl_spline_eval (myaoft, dect, acc);
	}

	// We start integrating after reheating
	// reheating scale factor
	double arh = gsl_spline_eval (myaoft, trh, acc);

	// We stop integrating when BH has decayed  or we reach the target time
	// final time
	double af = gsl_min(deca,a);

	// the result of integrating H^3((1+w)(1-3w)+w') actually has very short analytic result: -h0^2 Omega_m/(3a^3)
	//as well as conversion of Hubble^3 to kg and dt to 1/kg
	double hubint =  pow(h0,2) * om * (pow(af,3)-pow(arh,3))/(3.*pow(af*arh,3)); // integral in units of h0 = 1/s^2

	double planck2ds = pow(10.,2.*lgmp)*pow(kgtgev*gevtds,2); // planck mass ^2 in 1/s^2

	double chempot = -9*lambda*hubint/planck2ds; // unitless

	double qbl = (sumq2g/96./M_PI) * chempot; // unitless

	// psi term
	double exptpsi = psifinallg(lgm, peakm, ai, rem);
	double psi =  pow(10.,exptpsi); // units 1/kg

	return  intf * psi * qbl; // units 1/kg
}


double nbl_lcdm(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double peakm = pars.peakm;
	double lgn0 = pars.lgn0;
	double trh = pars.trh;
	double a = pars.aval;
	bool rem = pars.rem;

	// What's the time?
	double time = timeofa(a);
	// is time before reheating time? If yes, number density should be 0!
	// can't create asymmetry before reheating
	if (time<=trh) {
		return  0.;
	}

	// integral of Qbl x psi / M  ---> units 1/kg
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double qpsiint;
	double error;
	gsl_function F;
	F.function = &nbl_int_lcdm;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &qpsiint, &error);
	gsl_integration_workspace_free (w);

	// average initial mass times number density today in kg/m^3
	double mavgi = avgm(peakm, ai, rem);
	double rhot =   pow(10.,lgn0) * mavgi / pow(a,3) ;

	return rhot * qpsiint  ; // units of 1/m^3
}



/* Entropy density */
// Temperature given in Kelvin
double entropy(double Trh){
	return 2.*pow(M_PI,2)/45 * gstar * pow(Trh*keltgev,3.);
}


/* Yield */
// Params:
	// double lambda; // coupling constant - unitless
 	// double lgn0;  // Log10 of final total number density of bh in 1/m^3
	// double peakm; // log10 of mass spectrum peak in kg
	// double Trh; // reheating Temperature in Kelvin
	// double trh; // time of reheating in s
	// double aval; // scale factor
	// bool rem; // remnants or not
	// gsl_spline *spline;
	// gsl_interp_accel *acc;

/* Omega_baryon - lcdm background*/
double omegab(void * params){
		myparam_type2 pars = *(myparam_type2 *)(params);
		double a = pars.aval;
		double rhob = protonm * nbl_lcdm(params)/2.; // units:  kg/m^3
		return rhob /rhoc(a);
}

double yieldbl(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	// Temperature in Kelvin
	double Trh = pars.Trh;
	// entropy density in GeV^3
	double s = entropy(Trh);
	// nbl density converted to GeV^3
	double nblv = pow(dmtgev,3)*nbl_lcdm(params);
	// yield - unitless
	return nblv/s;
}



/* Hook's sec 3 B treatment - consistency test! */
/* We consider 1 epoch in which all black holes decay - matter domination! */
double nbl_hook(double  lgni, double lambda, double lgmass){
	// determine the decay time
	double dect = pow(10., (3.*lgmass-17.803)) ; // s

	// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
	double w0m = 0.;
	double eostm = (1.+w0m)*(1.-3.*w0m); // =1
	double rhobh = pow(10.,lgni+lgmass);
	double hubble2 = 8.*M_PI*gnewton/3. * rhobh;

	double hub3 = eostm * pow(hubble2, 3./2.); // 1/s^3

	double chempot = 9.*lambda*eostm* (hub3/pow(10.,2.*lgmp));

	double qfac = 13 * chempot * dect / 96 / M_PI * pow(dstgev/kgtgev,2);
	double nbh = pow(10.,lgni) * pow(dmtgev,3.); // BH number density in gev

	return   nbh *  qfac; // GeV^3
}
// calculate volume fration given number density
double epsilon(double lgni, double lgmass){
	double radius = 2.*pow(10.,lgmass)/pow(10.,lgmp*2.) /kgtgev;
	return pow(10.,lgni) * 4. * M_PI * pow(radius,3) / 3.  * pow(dmtgev,3);
}
// max allowed volume fraction (eq. 28 of 1404.0113)
double epsilonmax(double lgmass){
	return pow(0.08*pow(10.,lgmp)/pow(10.,lgmass),4);
}
/* Entropy density */
// Temperature given in GeV
// entropy in gev^3
double entropy_hook(double Trh){
	return 2.*pow(M_PI,2)/45 * gstar * pow(Trh,3.);
}
// Maximum reheating temperatture  in GeV
// Temperature of reheating (eq. 29 of 1404.0113)
double temp_hook(double epsilon, double lgmass){
	return pow(45.*epsilon*pow(10.,lgmp*6.)*pow(kgtgev,4.)/16./pow(M_PI,3.)/gstar/pow(10.,lgmass*2),1./4.);
}
//nbl: explicit (RHS) nbl from eq. 31 of 1404.0113
// in gev^3
double nblh(double lambda, double epsilon, double lgmass){
	double pmass = pow(10.,lgmp);
	double mass = pow(10.,lgmass);
	return 45./4./M_PI/gstar * 13 * lambda * pow(epsilon,5./2) * pow(pmass,6.)*pow(kgtgev/mass,3.);
}
// // explicit yield from eq. 32 of 1404.0113
double yblh(double epsilon, double lambda, double lgmass){
 	return  pow(3.,5./2.)*pow(5.,5./4.)/pow(M_PI,3./4.)*13./pow(gstar,5./4.) * lambda * pow(epsilon,7./4.) * pow(pow(10.,lgmp-lgmass),3./2.);
 }



// absolute maximum theoretical number density
double maxn0(double  lgmass){
	double rsch = 2.*gnewton*pow(10.,lgmass)/pow(sofl,2); // schwarzschild radius in m
	double volume =4./3. * M_PI * pow(rsch,3);
	return 1./volume * 0.74;
}


int main(int argc, char* argv[]) {

/* Hook [1404.0113] calculations - consistency tests */
// Some notes: A bh of mass 10^-5.96 kg has a volume of 1.76 x 10^-98 m^3  and so we can have number densities mathematically up to 10^98 x 0.74 where 0.74 is the close packing ratio (max ratio of cube to volume of packed spheres)
double lgmass = -5.966; // = 50 x planck mass
double lgni = 87.;
double lambda = pow(10.,12);

// // max allowed volume fraction (eq. 28 of 1404.0113)
// double epimax = epsilonmax(lgmass);
// std::cout<< "Max epsilon explicit: " << epimax << std::endl;
// std::cout<< "Max epsilon x lambda explicit: " << epimax*lambda  << std::endl;
// std::cout<< "" << std::endl;
//
// // calculate volume fration given number density
// double myep = epsilon(lgni, lgmass);
// std::cout<< "Epsilon: " << myep << std::endl;
// std::cout<< "Epsilon x lambda: " << myep*lambda << std::endl;
// std::cout<< "" << std::endl;
//
// // Temperature of reheating (eq. 29 of 1404.0113)
// double temph =  temp_hook(epimax, lgmass);
// std::cout<< "Max Temp of reheating: " << temph << std::endl;
// std::cout<< "" << std::endl;
//
// // entropy density (eq. 30 of 1404.0113)
// double ent = entropy_hook(temph);
//
// // nbl x nbh  from eq. 31 of 1404.0113
// double nbl = nbl_hook(lgni, lambda, lgmass);
// // yield
// double yb = nbl/ent;
//
// // explicit (RHS) nbl from eq. 31 of 1404.0113
// double nbl2 = nblh(lambda, epimax ,lgmass);
// double yb2 = nbl2/ent;
//
// // // explicit yield from eq. 32 of 1404.0113
//  //double ybe = yblh(epimax, lambda, lgmass);
//
// std::cout<< "Yield from our own calculation of Hook: " << yb << std::endl;
// std::cout<< "Yield from Hook  eq.32 : " << yb2 << std::endl;
// std::cout<< "" << std::endl;
//

// //double timeofa(3e-36);
// double hubcheck = pow(hubble(om,orad,3e-28),2);
// double hubhook  =  8.*M_PI*gnewton/3. *  pow(10.,lgni+lgmass);
//
// std::cout<< "LCDM Hubble:  " << hubcheck << std::endl;
// std::cout<< "Hook Hubble:  " << hubhook << std::endl;
// Hub^3 x t_dec ~ 9.351185e+41  s^-2 for reference


/* Our LCDM background calculations  */
// create spline of a(t)
arrays_T myxxyy = (arrays_T)malloc( sizeof(struct arrays) );
// populate array
aoftime(myxxyy);
// ceate spline
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *myspline    = gsl_spline_alloc (gsl_interp_cspline, 1002);
gsl_spline_init (myspline, (*myxxyy).xx, (*myxxyy).yy, 1002);

// Reset params but check consistency
// Reheating temp in GeV - check units in timeofrh function
// Can vary between ~ 10^14 GeV and 10 GeV -https://medium.com/predict/linking-cosmic-inflation-and-the-big-bang-with-reheating-period-1eb3f81526a1
double Trh = pow(10.,4.);  // 10 < Trh < 10^14

lgmass = 13.; // -7.67 < lgmass < 36
std::cout<< "Maximum number density: " << maxn0(lgmass) << std::endl;
std::cout<< "" << std::endl;

lgni = 39. ; // -100(none) < lgni < 98 (max from sphere-in-cube fitting problem)
lambda =  pow(10.,114); // any value

double myep = epsilon(lgni, lgmass);
std::cout<< "Epsilon x lambda: " << myep*lambda << std::endl;
std::cout<< "" << std::endl;
// reheating time can't be before inflation - bound on temperature
// set reheating time
double trh = timeofrh(Trh/keltgev); // reheating time
double arh = gsl_spline_eval (myspline, trh, acc); // reheating scale factor
std::cout<< "Reheating Temperature:  " << Trh << std::endl;
std::cout<< "Reheating time:  " << trh << std::endl;
std::cout<< "Reheating scale factor:  " << arh << std::endl;
std::cout<< "" << std::endl;

std::cout << "Results at CMB  " << std::endl;
double afin = acmb;
bool rem = false; // remnants or not
double lgn0 = 3.*log10(ai)  + lgni; // number density today
double omegabhh = omegapbh(lgn0, lgmass, afin, rem); // bh density fraction at CMB
struct myparam_type2 pars = {-lambda,lgn0,lgmass,Trh,trh, afin,rem,myspline,acc};
double omegabar =  omegab(&pars); // baryon density fraction at CMB

std::cout << "BH number density today:  " << pow(10.,lgn0) << std::endl;
std::cout << "BH density at CMB:  " << omegabhh << std::endl;
std::cout << "CDM density at CMB: " << omegalcdm(oc, afin) << std::endl;
std::cout <<  "Baryon density at CMB:  " << omegabar << std::endl;
std::cout << "Baryon density at CMB: " << omegalcdm(ob, afin) << std::endl;

double myy = yieldbl(&pars);
std::cout <<  "Yield at CMB :  " << myy << std::endl;
std::cout << "" << std::endl;


std::cout << "Results today  " << std::endl;
afin = 0.999;
omegabhh = omegapbh(lgn0, lgmass, afin, rem); // bh density fraction at CMB

struct myparam_type2 pars2 = {-lambda,lgn0,lgmass,Trh,trh, afin,rem,myspline,acc};
omegabar =  omegab(&pars2); // baryon density fraction at CMB

std::cout << "BH number density today:  " << pow(10.,lgn0) << std::endl;
std::cout << "BH density today:  " << omegabhh << std::endl;
std::cout << "CDM density today: " << omegalcdm(oc, afin) << std::endl;
std::cout <<  "Baryon density today:  " << omegabar << std::endl;
std::cout << "Baryon density today: " << omegalcdm(ob, afin) << std::endl;


myy = yieldbl(&pars2);
std::cout <<  "Yield today :  " << myy << std::endl;

return 0;
}



// double p1,p2,p3,p4,p5,p6;
// const char* output = "myoutput.dat";
// FILE* fp = fopen(output, "w");
//
// double lgn0,peakm,lam;
// // number of time points to output
// double amax =1.;
// double amin = acmb;
// double amax0 =acmb;
// double amin0 = ai;
// int na = 100;
// bool rem = false;
//
//
//   std::cout<<normlg(15, 0.1) << std::endl;
// 	std::cout<< pow(10,psifinallg(20, 20, 0.1, false)) << std::endl;
// 	std::cout<< avgm(0., 0.001, false) << std::endl;



//
// lgn0 = -62.01; // log10 of number density
// peakm = lgmf; // log10 of peak mass
// std::cout << omegapbh(lgn0, peakm, acmb)/0.695 << std::endl;
// lam = -1e60;
// //std::cout << Trh << std::endl;
// struct myparam_type2 pars = {lam,lgn0,peakm,trh,acmb,myspline,acc};
// std::cout << omegab(&pars)/0.129 << std::endl;
//
// struct myparam_type2 pars1 = {lam,lgn0,peakm,trh,1,myspline,acc};
// double ppar = yieldbl(&pars);
// std::cout << ppar << std::endl;


// for(int i=0; i<na; i++){
// 	double a = amin * exp(i*log(amax/amin)/(na-1));
//
// 			lgn0 = -49.794; // log10 of number density
// 			peakm = lgmf; // log10 of peak mass
// 			lam = -pow(10,54)/1.0555;
// 			struct myparam_type2 pars = {lam,lgn0,peakm,trh,a,myspline,acc};
// 			p1 = omegab(&pars);
//
// 			lgn0 = -50.366; // log10 of number density
// 			peakm = lgmf+18; // log10 of peak mass
// 			lam = -pow(10,42)/3.35;
// 			struct myparam_type2 pars2 = {lam,lgn0,peakm,trh,a,myspline,acc};
// 			p2 = omegab(&pars2);
//
// 			lgn0 = -51.357; // log10 of number density
// 			peakm = lgmf+20; // log10 of peak mass
// 			lam = -pow(10,44)/3.44;
// 			struct myparam_type2 pars3 = {lam,lgn0,peakm,trh,a,myspline,acc};
// 			p3 = omegab(&pars3);
//
//
// 		printf("%e %e %e %e  \n", a, p1,p2,p3);
//  		fprintf(fp,"%e %e %e %e \n",a, p1, p2, p3);
//  }
//

// free spline objects
// gsl_spline_free (myspline);
// free(myxxyy);
// gsl_interp_accel_free (acc);

/*close output file*/
 //fclose(fp);
