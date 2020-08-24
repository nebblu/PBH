
#include <constants.h>
#include <gfuncs.h>
#include <md.h>


//Final, normalised log10 of broad spectrum
double psifinallg(double lgm, void * params){
		myparam_type2 pars = *(myparam_type2 *)(params);
		double peakm = pars.peakm;
		//  Log10[M(t)]
		double lgmoft = bhmasslg(lgm, params);
		double expt = psibroadlg(lgmoft, peakm) + normlg(params) ;
		return expt;
}

// Average mass
double avgm_int(double lgm, void * params){
	double expt = 2.*lgm +  psifinallg(lgm, params);
	return pow(10.,expt);
}

// average mass (see avgm_int function) in kg
double avgm(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &avgm_int;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);
	return intf * result ;
}


/* LCDM background treatment */
// we assume Hubble  = LCDM with ob,ocdm,orad after reheating, and matter dom before reheating
/*  The B-L number density integrand */
double nbl_int_lcdm_rh(double lgm, void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double lambda = pars.lambda;
	double peakm = pars.peakm;
	double trh = pars.trh;
	double lgn0 = pars.lgn0;
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

	double rhobh, hubble2, aftrh, befrh, mavg;
	// if decay time is after reheating we add contributions from  initial matter phase and normal LCDM expansion
	if (dect>trh) {
		//integrate before reheating:
		// average initial mass at reheating time times number density today in kg/m^3
		//mavg = pow(10.,peakm);//avgm(params);
		// matter density in kg/m^3
		rhobh = pow(10.,lgn0+peakm)/pow(arh,3) ;
		// hubble squared
		hubble2 = 8.*M_PI*gnewton/3. * rhobh;
		befrh = pow(hubble2, 3./2.) * trh;

		// integrate after reheating :
		// the result of integrating H^3((1+w)(1-3w)+w') actually has very short analytic result: -h0^2 Omega_m/(3a^3)
		//as well as conversion of Hubble^3 to kg and dt to 1/kg
		aftrh =  pow(h0,2) * om * (pow(af,3)-pow(arh,3))/(3.*pow(af*arh,3)); // integral in units of h0 = 1/s^2
	}
	// if decay time is less than reheating, we only have matter domination contribution
	else{
		// average initial mass times number density today in kg/m^3
		//mavg = avgm(params);
		// matter density in kg/m^3
		rhobh = pow(10.,lgn0+peakm)/pow(deca,3) * mavg;
		// hubble squared
		hubble2 = 8.*M_PI*gnewton/3. * rhobh;
		// integration assuming constant hubble
		befrh = pow(hubble2, 3./2.) * dect;
		aftrh = 0.;
	}


	double planck2ds = pow(10.,2.*lgmp)*pow(kgtgev*gevtds,2); // planck mass ^2 in 1/s^2

	double chempot = -9*lambda*(befrh + aftrh)/planck2ds; // unitless

	double qbl = (sumq2g/96./M_PI) * chempot; // unitless

	// psi term
	double exptpsi = psifinallg(lgm, params);
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
	// initial matter domination phase included
// LCDM expansion with initial matter domination era
double nbl_lcdm_rh(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double lgn0 = pars.lgn0;
	double a = pars.aval;

	// integral of Qbl x psi / M  ---> units 1/kg
	gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
	double qpsiint;
	double error;
	gsl_function F;
	F.function = &nbl_int_lcdm_rh;
	F.params = &pars;
	gsl_integration_qags (&F, lgmp, lgmax, 0, 1e-7, 1000, w, &qpsiint, &error);
	gsl_integration_workspace_free (w);

	// average initial mass times number density today in kg/m^3
	double mavgi = avgm(params);
	double rhot =   pow(10.,lgn0) * mavgi / pow(a,3) ;

	return rhot * qpsiint  ; // units of 1/m^3
}



/* Omega_pbh */
// lgn0 = Log10[n0] where n0 is the number density of bh today
double omegapbh(void * params){
	myparam_type2 pars = *(myparam_type2 *)(params);
	double lgn0 = pars.lgn0;
	double a = pars.aval;
	double avgmass = avgm(params);
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
double yieldbl(double Trh, double nbl, double a){
	// entropy density in GeV^3
	double s = entropy(Trh, a);
	// nbl density converted to GeV^3
	double nblv = pow(dmtgev,3)*nbl;
	// yield - unitless
	return nblv/s;
}


int main(int argc, char* argv[]) {
/* Our LCDM background calculations  with matter domination pre reheating */
//  Input parameters
double lgmass = 12.; // -7.67 < lgmass < 36
double lgni = 39.48 ; // -100(none) < lgni < 98 (max from sphere-in-cube fitting problem if planck  mass avg)
double lambda =  pow(10.,73)/4./3; // any value
double Trhgev = pow(10.,5);  // 10 < Trh < 10^24  GeV https://medium.com/predict/linking-cosmic-inflation-and-the-big-bang-with-reheating-period-1eb3f81526a1
double afin = acmb;
bool rem = true; // remnants or not

double  trh, Trh, myep, arh, lgn0, omegabhh, nbl, omegabar, myy;

std::cout<< "Maximum number density: " << maxn0(lgmass) << std::endl;
std::cout<< "" << std::endl;

myep = epsilon(lgni, lgmass);
std::cout<< "Epsilon x lambda: " << myep*lambda << std::endl;
std::cout<< "" << std::endl;


Trh = Trhgev/keltgev;
trh = timeofrh(Trh); // reheating time
lgn0 = 3.*log10(ai)  + lgni; // number density today

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *myspline    = gsl_spline_alloc (gsl_interp_cspline, 1002);

// Load  parameters into structure
struct myparam_type2 pars = {-lambda,lgn0,lgmass,Trh,trh, afin,rem,myspline,acc};

//create spline of a(t) which depends on peak mass, number density  and reheating temperature/time
arrays_T myxxyy = (arrays_T)malloc( sizeof(struct arrays) );
// populate array
aoftime(myxxyy, pars);

// ceate spline
gsl_spline_init (myspline, (*myxxyy).xx, (*myxxyy).yy, 1002);

arh = gsl_spline_eval (myspline, trh, acc); // reheating scale factor
std::cout<< "Reheating Temperature [K]:  " << Trh << std::endl;
std::cout<< "Reheating Temperature [GeV]:  " << Trhgev << std::endl;
std::cout<< "Reheating time:  " << trh << std::endl;
std::cout<< "Reheating scale factor:  " << arh << std::endl;
std::cout<< "" << std::endl;


// for(int i=0; i<1000; i++){
// 	double a = ai*exp(i*log(0.9999/ai)/(1000.));
// 	std::cout << a << " " << hubblemd(a,&pars)/hubble(om, orad, a) << std::endl;
//
// }

std::cout << "Results at CMB  " << std::endl;

std::cout << "BH number density today:  " << pow(10.,lgn0) << std::endl;

omegabhh = omegapbh(&pars); // bh density fraction at CMB
std::cout << "BH density at CMB:  " << omegabhh << std::endl;
std::cout << "CDM density at CMB: " << omegalcdm(oc, afin) << std::endl;
std::cout<< "" << std::endl;

nbl = nbl_lcdm_rh(&pars);
omegabar =  omegab(nbl, afin); // baryon density fraction at CMB
std::cout <<  "Baryon density at CMB:  " << omegabar << std::endl;
std::cout << "Baryon density at CMB: " << omegalcdm(ob, afin) << std::endl;

// myy = yieldbl(Trh, nbl, afin);
// std::cout <<  "Yield at CMB :  " << myy << std::endl;
// std::cout << "" << std::endl;

return 0;
}
