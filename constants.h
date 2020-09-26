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
// Radiation fraction from Planck 2018 + massless neutrinos
const double orad = 9.2367e-5;//0.00006;
// Yield of B-L = n_B-L/ entropy density
const double yield = 1e-10;
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
const double lgmax =36.;
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
} ;

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
