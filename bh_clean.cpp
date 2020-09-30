
#include <constants.h>
#include <gfuncs.h>
#include <clean.h>


//Final, normalised log10 of broad spectrum
double psifinallg(double lgm, double peakm){
		double expt = psibroadlg(lgm, peakm) + normlg(peakm) ;
		return expt;
}

// Average mass
double avgm_int(double lgm, void * params){
	myparam_type pars = *(myparam_type *)(params);

	double peakm = pars.peakm;
	double a = pars.aval;
	bool rem = pars.rem;

	double lgmf =  bhmasslg(lgm,a,rem);

	double expt = lgmf + lgm +  psifinallg(lgm, peakm);

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

	double lgmaximum = bhmasslg(lgmax,a,rem);
	double logdect = 1./3. * (17.803 + log10(timeofa(a)));
	double lgminimum = gsl_max(lgmp, logdect);

	gsl_integration_qags (&F, lgminimum, lgmaximum, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);
	return intf * result ;
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
	double exptpsi = psifinallg(lgm, peakm);
	double psi =  pow(10.,exptpsi); // units 1/kg

	return  intf * psi * qbl; // units 1/kg
}
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
	// matd : True - includes initial matter domination, false - we exclude initial phase and begin integrating at end of reheating
double nbl_lcdm(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double peakm = pars.peakm;
	double lgn0 = pars.lgn0;
	double trh = pars.trh;
	double a = pars.aval;
	bool rem = pars.rem;

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


/* Omega_pbh */
// lgn0 = Log10[n0] where n0 is the number density of bh today
double omegapbh(double lgn0, double peakm, double a, bool rem){
	double avgmass = avgm(peakm, a, rem);
	return  pow(10.,lgn0) * avgmass / (rhoc(a) * pow(a,3.));
}

/* Omega_baryon - lcdm background*/
//  nbl  in m^-3
double omegab(double nbl, double a){
		double rhob = protonm * nbl/2.; // units:  kg/m^3
		return rhob /rhoc(a);
}


/* Yield */
// Trh in  Kelvin and  nbl  in m^-3
double yieldbl(double Trh, double nbl, double a, double arh){
	// entropy density in GeV^3
	double s = entropy(Trh, a, arh);
	// nbl density converted to GeV^3
	double nblv = pow(dmtgev,3)*nbl;
	// yield - unitless
	return nblv/s;
}



int main(int argc, char* argv[]) {
/* Our LCDM background calculations  */

/* Set key parameters */

double lgpeak = log10(1e5); // log10 of mass function peak in kg
bool rem = false; // remnants or not
double Trhgev = pow(10.,8);  // reheating temperature 1 < Trh < 10^24

/* create spline of a(t) */
arrays_T myxxyy = (arrays_T)malloc( sizeof(struct arrays) );
// populate array
aoftime(myxxyy);
// ceate spline
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *myspline    = gsl_spline_alloc (gsl_interp_cspline, 1002);
gsl_spline_init (myspline, (*myxxyy).xx, (*myxxyy).yy, 1002);

double  myep, trh, arh, Trh, lgn0, lgncmb, lgni, omegabh, nbl, omegabar, yield, avgmass, cdmfrac, barfrac, lambda;

Trh = Trhgev/keltgev;
trh = timeofrh(Trh); // reheating time
arh = gsl_spline_eval (myspline, trh, acc); // reheating scale factor

/* Set number density and lambda to match CMB observations */

avgmass = log10(avgm(lgpeak, acmb , rem));
/*  What number density gives thhe correct PBH fraction at CMB (= CDM fraction) */
cdmfrac = omegalcdm(oc, acmb);
lgni = log10(cdmfrac*rhoc(acmb)) + 3.*(log10(acmb)-log10(ai)) - avgmass; // -100(none) < lgni < 98 (max from sphere-in-cube fitting problem if planck  mass avg)
lgncmb = lgni + 3.*(log10(ai)-log10(acmb)); // number density at acmb
lgn0 = lgni + 3.*log10(ai); // number density today


// Take lambda to be the value needed by CMB Omega_b constraint and chosen T_RH
lambda =  1.;
struct myparam_type2 pars = {-lambda,lgn0,lgpeak,Trh,trh,acmb,rem,myspline,acc};

nbl = nbl_lcdm(&pars);
omegabar = omegab(nbl,  acmb); // baryon density fraction at CMB
barfrac = omegalcdm(ob, acmb);

/* reset lambda to needed value */
lambda = barfrac/omegabar;

/* Volume fraction of bH at CMB */
myep = epsilon(lgncmb, avgmass);

/* Parameter output */
std::cout<< "Log10 peak mass :  " << lgpeak << std::endl;
std::cout<< "Remnants?  " << rem << std::endl;
std::cout<< "Reheating Temperature [GeV]:  " << Trhgev << std::endl;
std::cout<< "Reheating time:  " << trh << std::endl;
std::cout<< "Reheating scale factor:  " << arh << std::endl;
std::cout<< "" << std::endl;
std::cout<< "Number density at initial time : " << pow(10,lgni) << std::endl;
std::cout<< "Number density at CMB : " << pow(10,lgncmb) << std::endl;
std::cout<< "Maximum number density at CMB: " << maxn0(avgmass) << std::endl;
std::cout << "BH number density today:  " << pow(10.,lgn0) << std::endl;
std::cout<< "" << std::endl;
std::cout<< "lambda: " << lambda << std::endl;
std::cout<< "volume fraction epsilon at CMB: " << myep << std::endl;
std::cout<< "Epsilon x lambda at CMB: " << myep*lambda << std::endl;
std::cout<< "" << std::endl;
std::cout<< "" << std::endl;
std::cout<< "" << std::endl;
std::cout<< "Temp x scale factor @ reheating:" << arh *  Trh << std::endl;
std::cout<< "" << std::endl;
std::cout<< "" << std::endl;

double avgmassi = avgm(lgpeak, ai , rem);
double avgmassf = avgm(lgpeak, acmb , rem);
double bhloss = pow(10.,lgn0)*(avgmassi - avgmassf);
double bargain = nbl * protonm /2. * pow(acmb,3) ;
std::cout << " energy lost by BH" << bhloss << std::endl;
std::cout << " energy gained by baryons" << bargain << std::endl;
std::cout << "Difference " << bhloss-bargain << std::endl;
std::cout << "Minimum radiation energy density @ z=0 is " << (bhloss-bargain)/rhoc(1.) << std::endl;

std::cout<< "" << std::endl;
std::cout<< "" << std::endl;
std::cout<< "Time of matter/radiation equality:  " << gsl_spline_eval (myspline, 10, acc)  << std::endl;


/* Output to file */
const char* output = "data/massspec_peak_1e5.dat";
FILE* fp = fopen(output, "w");
int nout = 10000; // number of linearly spaced output points

double ainit = arh;//1e-6; // initial scale factor to start output
double afin = 1.; // final scale factor

for(int i=0; i< nout; i++){
	/* MASS FUNCTION */
	double lgmass = lgmp + i*(lgmax-lgmp)/(nout-1.); //lgmp * exp(i*log(lgmax/lgmp)/(nout-1.));

	double p1 = psifinallg(lgmass, lgpeak);

	double p2 = psilowlg(lgmass,lgpeak);
	double p3 = psihilg(lgmass, oc);

	printf("%e %e %e %e  \n", lgmass, p1,p2,p3);
	fprintf(fp,"%e %e %e %e  \n", lgmass, p1,p2,p3);

	/* Density fractions and Yield  */
	// double scalef = ainit* exp(i*log(afin/ainit)/(nout-1.));
	//
	// omegabh = omegapbh(lgn0, lgpeak, scalef, rem); // bh density fraction
	//
	// struct myparam_type2 pars = {-lambda,lgn0,lgpeak,Trh,trh,scalef,rem,myspline,acc};
	//
	// nbl = nbl_lcdm(&pars);
	// omegabar =  omegab(nbl,  scalef); // baryon density fraction at CMB
	//
	// yield = yieldbl(Trh, nbl, scalef, arh);
	//
	// printf("%e %e %e %e %e %e   \n", scalef, omegabh, omegabar, yield, lambda, pow(10,lgn0));
	// fprintf(fp,"%e %e %e %e %e %e \n", scalef, omegabh, omegabar, yield, lambda, pow(10,lgn0));

	/* Average Mass */
	// double avgmi = avgm(lgmass, ai , rem);
	// double avgmf = avgm(lgmass, acmb , rem);

	 // printf("%e %e %e %e \n", lgmass, avgmi,avgmf, avgmi-avgmf);
	 // fprintf(fp,"%e %e %e \n", lgmass, avgmi,avgmf);

}


return 0;
}



 // Unused epoch treatment
 /* Epoch treatment - probably won't present these results in final draft */
 /* We consider 3 epochs : ai to arad , arad to alam,  alam to 1 */
 /*  The B-L number density integrand */
 // double nbl_int_ep(double lgm, void * params){
 // 	myparam_type2 pars = *(myparam_type2 *)(params);
 // 	double lambda = pars.lambda;
 // 	double peakm = pars.peakm;
 // 	double trh = pars.aval;
 // 	double a = pars.aval;
 // 	bool rem = pars.rem;
 // 	gsl_spline *myaoft = pars.spline;
 // 	gsl_interp_accel *acc = pars.acc;
 //
 // 	// determine the decay time
 // 	double dect = pow(10., (3.*lgm-17.803));
 // 	// so spline doesn't run into issues, set deca = 1 if the decay time exceeds lifetime of the universe....
 // 	double deca;
 // 	if (dect>=1./hubble(om,orad,1)) {
 // 		deca = 1.;
 // 	}
 // 	else if(dect<=trh){
 // 			return 0.; // no chemical potential before reheatingg time so no baryons from BH that have already decayed
 // 		}
 // 	else{
 // 		deca = gsl_spline_eval (myaoft, dect, acc);
 // 	}
 //
 // 	// reheating scale factor
 // 	double arh = gsl_spline_eval (myaoft, trh, acc);
 //
 // 	// chem  pot and q prefactor
 // 	 double prefac = -9.*sumq2g/(96.*M_PI *pow(10.,2.*lgmp));
 //
 // 	 /*Add various charge contributions */
 //  	double lgqr, lgqm;
 //  	double exptq1,exptq2,exptq3,exptq4;
 //  	double w0r, w0m, hubr, hubm, qfac;
 //
 // 	if (a<=arad) {
 // 		double af = gsl_min(deca,a);
 // 		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
 // 		w0r = (1.-1e-3)/3.;
 // 		double eost = (1.+w0r)*(1.-3.*w0r);
 // 		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
 // 		double a2 = pow(af,2);
 // 		double a3 = af*a2;
 // 		double ai2 = pow(arh,2);
 // 		double ai3 = arh*a2;
 // 		double ol = 1.-om-orad;
 // 		hubr = 1./6. * pow(h0,2) * ( 3.*om * ((a2-ai2)/(a2*ai2))  + 2.*orad * (a3-ai3/(a3*ai3)) + 6.*ol*(af-arh));
 // 		// //printf("%e \n", lghubr);
 // 		// double lgmur = 15.0783  +  lghubr + lglam; // 15.0783 = log10[9/Mp^2 * 19 / 96/pi]
 // 		qfac = lambda*prefac*eost*hubr; //pow(10.,lgmur);
 // 	}
 // 	// same thing but if we run into matter domination
 // 	else{
 // 		/* Radiation terms */
 // 		double af = gsl_min(arad,deca);
 //
 // 		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
 // 		w0r = (1.-1e-3)/3.;
 // 		double eost = (1.+w0r)*(1.-3.*w0r);
 // 		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
 // 		double a2 = pow(af,2);
 // 		double a3 = af*a2;
 // 		double ai2 = pow(arh,2);
 // 		double ai3 = arh*a2;
 // 		double ol = 1.-om-orad;
 // 		hubr =eost*(1./6. * pow(h0,2) * ( 3.*om * ((a2-ai2)/(a2*ai2))  + 2.*orad * (a3-ai3/(a3*ai3)) + 6.*ol*(af-arh)));
 //
 // 		/*Matter terms */
 // 		double afm;
 // 		// check if bh has decayed before mat-de equality
 // 		afm = gsl_min(a,alam);
 // 		// check if bh has decayed before final time
 // 		afm = gsl_min(afm,deca);
 // 		// check if black hole has decayed before arad
 // 		afm = gsl_max(afm,arad);
 // 		// equation of state during radiation domination --- equation 19 of 1404.0113. w'=0
 // 		w0m = 0.;
 // 		double eostm = (1.+w0m)*(1.-3.*w0m);
 // 		// Integrating eq.14 of 1404.0113 with constant eos over scale factor (just an integral of hubble^2)
 // 		double a2m = pow(af,2);
 // 		double a3m = af*a2m;
 // 		double ai2m = pow(arad,2);
 // 		double ai3m = arad*ai2m;
 // 		hubm = eostm*(1./6. * pow(h0,2) * ( 3.*om * ((a2m-ai2m)/(a2m*ai2m))  + 2.*orad * (a3m-ai3m/(a3m*ai3m)) + 6.*ol*(afm-arad)));
 //
 // 		qfac = lambda*prefac*(hubr + hubm);
 // 	}
 //
 // 	// psi term
 // 	double exptpsi = psifinallg(lgm, peakm, ai, rem);
 //
 // 	return  intf * pow(10.,exptpsi) * qfac;
 // }
 //
 //
 // double nbl_ep(void * params){
 // 	myparam_type2 pars = *(myparam_type2 *)(params);
 // 	double peakm = pars.peakm;
 // 	double lgn0 = pars.lgn0;
 // 	double trh = pars.trh;
 // 	double a = pars.aval;
 // 	bool rem = pars.rem;
 //
 // 	// can't create asymmetry before reheating
 // 	double time = timeofa(a);
 // 	if (time<=trh) {
 // 		return  0.;
 // 	}
 //
 // 	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
 // 	double result, error;
 // 	gsl_function F;
 // 	F.function = &nbl_int_ep;
 // 	F.params = &pars;
 //
 // 	double mavgi = avgm(peakm, ai, rem);
 // 	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &result, &error);
 // 	gsl_integration_workspace_free (w);
 //
 // 	return  pow(10.,lgn0) * mavgi * result / pow(a,3);
 // }
 //
 //
 // /* Omega_baryon for epoch treatment */
 // double omegab_ep(void * params){
 // 		myparam_type2 pars = *(myparam_type2 *)(params);
 // 		double a = pars.aval;
 // 		double rhob = protonm * nbl_ep(params)/2.;
 // 		return rhob /rhoc(a);
 // }
