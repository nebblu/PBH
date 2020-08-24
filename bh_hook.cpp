
#include <constants.h>
#include <gfuncs.h>

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

		printf("%e \n", sqrt(hubble2));

	double hub3 = eostm * pow(hubble2, 3./2.); // 1/s^3

	double chempot = 9.*lambda * (hub3/pow(10.,2.*lgmp));

	double qfac = 13 * chempot * dect / 96 / M_PI * pow(dstgev/kgtgev,2);
	double nbh = pow(10.,lgni) * pow(dmtgev,3.); // BH number density in gev

	return   nbh *  qfac; // GeV^3
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



int main(int argc, char* argv[]) {

/* Hook [1404.0113] calculations - consistency tests */
// Some notes: A bh of mass 10^-5.96 kg has a volume of 1.76 x 10^-98 m^3  and so we can have number densities mathematically up to 10^98 x 0.74 where 0.74 is the close packing ratio (max ratio of cube to volume of packed spheres)
double lgmass = -5.966; // = 50 x planck mass
double lgni = 86.59;
double lambda = pow(10.,12);

// max allowed volume fraction (eq. 28 of 1404.0113)
double epimax = epsilonmax(lgmass);
std::cout<< "Max epsilon explicit: " << epimax << std::endl;
std::cout<< "Max epsilon x lambda explicit: " << epimax*lambda  << std::endl;
std::cout<< "" << std::endl;

// calculate volume fration given number density
double myep = epsilon(lgni, lgmass);
std::cout<< "Epsilon: " << myep << std::endl;
std::cout<< "Epsilon x lambda: " << myep*lambda << std::endl;
std::cout<< "" << std::endl;

// Temperature of reheating (eq. 29 of 1404.0113)
double temph =  temp_hook(epimax, lgmass);
std::cout<< "Max Temp of reheating: " << temph << std::endl;
std::cout<< "" << std::endl;

// entropy density (eq. 30 of 1404.0113)
double ent = entropy_hook(temph);

// nbl x nbh  from eq. 31 of 1404.0113
double nbl = nbl_hook(lgni, lambda, lgmass);
// yield
double yb = nbl/ent;

// explicit (RHS) nbl from eq. 31 of 1404.0113
double nbl2 = nblh(lambda, epimax ,lgmass);
double yb2 = nbl2/ent;

// // explicit yield from eq. 32 of 1404.0113
 double ybe = yblh(epimax, lambda, lgmass);

std::cout<< "Yield from our own calculation of Hook: " << yb << std::endl;
std::cout<< "Yield from Hook  eq.31 : " << yb2 << std::endl;
std::cout<< "Yield from Hook  eq.32 : " << ybe << std::endl;
std::cout<< "" << std::endl;


return 0;
}
