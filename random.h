/****************************  random.h   *********************************
* Author:        Agner Fog
* Date created:  2023-07-05
* Last modified: 2023-12-31
* Version:       3.00.00
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This header file defines random number generators with uniform and 
* non-uniform distributions
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#pragma once
#include <inttypes.h>


const int MAXCOLORS = 10;    // Maximum number of colors in multivariate distributions

/***********************************************************************
*  Define random number generator classes
***********************************************************************/

/***********************************************************************
*  Mersenne twister. Has very long cycle length
***********************************************************************/
class RandomMersenne {
// Choose which version of Mersenne Twister you want:
protected:
    enum par {
    #if false 
        // Define constants for type MT11213A:
        MERS_N = 351,
        MERS_M = 175,
        MERS_R = 19,
        MERS_U = 11,
        MERS_S = 7,
        MERS_T = 15,
        MERS_L = 17,
        MERS_A = 0xE4BD75F5,
        MERS_B = 0x655E5280,
        MERS_C = 0xFFD58000
    #else    
        // or constants for type MT19937:
        MERS_N = 624,
        MERS_M = 397,
        MERS_R = 31,
        MERS_U = 11,
        MERS_S = 7,
        MERS_T = 15,
        MERS_L = 18,
        MERS_A = 0x9908B0DF,
        MERS_B = 0x9D2C5680,
        MERS_C = 0xEFC60000
    #endif
    };
public:
   RandomMersenne(int32_t seed) {      // Constructor
      randomInit(seed); lastInterval = 0;}
   void randomInit(int32_t seed);      // Re-seed
   void randomInitByArray(int32_t const seeds[], int numSeeds); // Seed by more than 32 bits
   int32_t iRandom (int32_t min, int32_t max); // Output random integer
   int32_t iRandomX(int32_t min, int32_t max); // Output random integer, exact
   double random();                    // Output random float
   uint32_t bRandom();                 // Output random bits
private:
   void init0(int32_t seed);           // Basic initialization procedure
   uint32_t mt[par::MERS_N];           // State vector
   int mti;                            // Index into mt
   uint32_t lastInterval;              // Last interval length for iRandomX
   uint32_t rLimit;                    // Rejection limit used by iRandomX
};


/***********************************************************************
* Mother-of-all multiply-with-carry random number generator
* Has high bifurcation
***********************************************************************/

class RandomMother {                   // Encapsulate random number generator
public:
   RandomMother(int seed) {            // Constructor
      randomInit(seed);}
   void randomInit(int seed);          // Initialization
   int iRandom(int min, int max);      // Get integer random number in desired interval
   double random();                    // Get floating point random number
   uint32_t bRandom();                 // Output random bits
protected:
   uint32_t x[5];                      // History buffer
};


/***********************************************************************
*  Combination of two generators. Combines advantages of both
***********************************************************************/

class RandomCombined : protected RandomMersenne, protected RandomMother {
public:
   RandomCombined(int seed) : RandomMersenne(seed), RandomMother(seed+1) {
       randomixInterval = randomixLimit = 0;
   };
   void randomInit(int seed) {                   // re-seed
      RandomMersenne::randomInit(seed);
      RandomMother::randomInit(seed+1);
   }
   uint32_t bRandom() {                          // Output random bits
       return RandomMersenne::bRandom() + RandomMother::bRandom();
   }
   double random();                              // Output random float, double precision
   float  randomf();                             // Output random float, single precision
   int32_t iRandom (int32_t min, int32_t max);   // Output random integer
   int32_t iRandomX(int32_t min, int32_t max);   // Output random integer, exact
protected:
    uint32_t randomixInterval;                   // used by function IRandomX
    uint32_t randomixLimit;                      // used by function IRandomX
};


/**************************************************************************
  Non-uniform random variate generator
This can generate variates with several different probability distributions
**************************************************************************/

//  choose which uniform random number generator to use as base for non-uniform generator:
#define RANDOM_BASE RandomCombined

class RandomVariates : public RANDOM_BASE {
public:
    // constructor
    RandomVariates(int32_t seed = 0);            // constructor
    void randomInit(int32_t seed) {              // re-seed. Must override RANDOM_BASE::randomInit
        RANDOM_BASE::randomInit(seed);           // initialize underlying generator
        // discard any pending normal_x2 to make sure sequence is reproducible
        normal_x2 = 0.; normal_x2_valid = false;
    }
    void reportError(const char * e);            // Print error message if parameters out of range
    //const char * test(int testcase);           // Test random variate generator and print results
    bool    bernoulli(double p);                 // Bernoulli distribution
    double  normal(double m, double s);          // Normal distribution
    double  normalTrunc(double m, double s, double limit1, double limit2); // Truncated normal distribution
    int32_t poisson (double L);                  // Poisson distribution
    int32_t binomial (int32_t n, double p);      // Binomial distribution
    int32_t hypergeometric (int32_t n, int32_t m, int32_t N); // Hypergeometric distribution
    int32_t walleniusNCHyp (int32_t n, int32_t m, int32_t N, double odds); // Wallenius noncentral hypergeometric distribution
    int32_t complWalleniusNCHyp (int32_t n, int32_t m, int32_t N, double odds); // Complementary Wallenius noncentral hypergeometric distribution
    int32_t fishersNCHyp (int32_t n, int32_t m, int32_t N, double odds); // Fisher's noncentral hypergeometric distribution
    void    multinomial (int32_t * destination, double * p, int32_t n, int colors); // Multinomial distribution
    void    multiHypergeometric(int32_t * destination, int32_t * source, int32_t n, int colors); // Multivariate hypergeometric distribution
    void    multiWalleniusNCHyp(int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Multivariate Wallenius noncentral hypergeometric distribution
    void    complMultiWalleniusNCHyp(int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Complementary Multivariate Wallenius noncentral hypergeometric distribution
    void    multiFishersNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors); // Multivariate Fisher's noncentral hypergeometric distribution
    void    shuffle(int32_t * list, int32_t n1, int32_t n2);// make a list of numbers from n1 to n2 in random order
private:
    // Variables used by function normal
    bool   normal_x2_valid;
    double normal_x2;
    // Variables used by function poissonInver
    double pois_L_last;                          // previous value of L
    double pois_f0;                              // value at x=0 or at mode
    // Variables used by function poissonRatioUniforms
    double pois_a;                               // hat center
    double pois_h;                               // hat width
    double pois_g;                               // ln(L)
    int32_t pois_bound;                          // upper bound
    // Variables used by binomialRatioOfUniforms
    int32_t bino_n_last;                         // last n
    double bino_p_last;                          // last p
    int32_t bino_mode;                           // mode
    int32_t bino_bound;                          // upper bound
    double bino_a;                               // hat center
    double bino_h;                               // hat width
    double bino_g;                               // value at mode
    double bino_r1;                              // p/(1-p) or ln(p/(1-p))
    // Variables used by hypRatioOfUnifoms
    int32_t  hyp_n_last, hyp_m_last, hyp_N_last; // Last values of parameters
    int32_t  hyp_mode, hyp_mp;                   // Mode, mode+1
    int32_t  hyp_bound;                          // Safety upper bound
    double hyp_a;                                // hat center
    double hyp_h;                                // hat width
    double hyp_fm;                               // Value at mode
    // Variables used by walleniusNCHypRatioOfUnifoms
    int32_t wnc_n_last, wnc_m_last, wnc_N_last;  // previous parameters
    double wnc_o_last;
    int32_t wnc_bound1, wnc_bound2;              // lower and upper bound
    int32_t wnc_mode;                            // mode
    double wnc_a;                                // hat center
    double wnc_h;                                // hat width
    double wnc_k;                                // probability value at mode
    bool useChopDown;                            // use chop down inversion instead
    // Variables used by fishersNCHyp and its subfunctions:
    int32_t fnc_n_last, fnc_m_last, fnc_N_last;  // last values of parameters
    int32_t fnc_bound;                           // upper bound
    double fnc_o_last;
    double fnc_f0, fnc_scale;
    double fnc_a;                                // hat center
    double fnc_h;                                // hat width
    double fnc_lfm;                              // ln(f(mode))
    double fnc_logb;                             // ln(odds)

protected:
    double  accuracy;                                           // desired accuracy of calculations for walleniusNCHyp
    int32_t poissonLow(double L);                               // Subfunction used by poisson function
    int32_t poissonInver(double L);                             // Subfunction used by poisson function
    int32_t poissonRatioUniforms(double L);                     // Subfunction used by poisson function
    int32_t binomialInver (int32_t n, double p);                // Subfunction used by binomial function
    int32_t binomialRatioOfUniforms (int32_t n, double p);      // Subfunction used by binomial function
    int32_t hypInversionMod (int32_t n, int32_t M, int32_t N);  // Subfunction used by hypergeometric function
    int32_t hypRatioOfUnifoms (int32_t n, int32_t M, int32_t N);// Subfunction used by hypergeometric function
    int32_t walleniusNCHypUrn(int32_t n, int32_t m, int32_t N, double odds); // Subfunction used by walleniusNCHyp
    int32_t walleniusNCHypInversion(int32_t n, int32_t m, int32_t N, double odds); // Subfunction used by walleniusNCHyp
    int32_t walleniusNCHypRatioOfUnifoms(int32_t n, int32_t m, int32_t N, double odds); // Subfunction used by walleniusNCHyp
    int32_t fishersNCHypInversion (int32_t n, int32_t m, int32_t N, double odds); // Subfunction used by fishersNCHyp
    int32_t fishersNCHypRatioOfUnifoms (int32_t n, int32_t m, int32_t N, double odds); // Subfunction used by fishersNCHyp
};


/***********************************************************************
    Class WalleniusNCHypergeometric
This class contains methods for calculating the univariate
Wallenius' noncentral hypergeometric probability function.
***********************************************************************/

class WalleniusNCHypergeometric {
public:
   WalleniusNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy=1.E-6); // constructor
   void setParameters(int32_t n, int32_t m, int32_t N, double odds); // change parameters
   double probability(int32_t x);                          // calculate probability function
   //int32_t MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                                      // approximate mean
   double variance(void);                                  // approximate variance (poor approximation)
   int32_t mode(void);                                     // calculate mode
   //double moments(double * mean, double * var);          // calculate exact mean and variance
   int bernouilliH(int32_t x, double h, double rh, RandomVariates *sto); // used by rejection method

   // implementations of different calculation methods
protected:
   double recursive(void);                                 // recursive calculation
   double binoexpand(void);                                // binomial expansion of integrand
   double laplace(void);                                   // Laplace's method with narrow integration interval
   double integrate(void);                                 // numerical integration

   // other subfunctions
   double lnbico(void);                                    // natural log of binomial coefficients
   void findpars(void);                                    // calculate r, w, E
   double integrate_step(double a, double b);              // used by integrate()
   double search_inflect(double t_from, double t_to);      // used by integrate()
   void reportError(const char* text);                     // report an error message

   // parameters
   double omega;                                           // Odds
   int32_t n, m, N, x;                                     // Parameters
   int32_t xmin, xmax;                                     // Minimum and maximum x
   double accuracy;                                        // Desired precision
   // parameters used by lnbico
   int32_t xLastBico;
   double bico, mFac, xFac;
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   int32_t xLastFindpars;
};


/***********************************************************************
Class MultiWalleniusNCHypergeometric
This class contains methods for calculating the multivariate
Wallenius' noncentral hypergeometric probability function.
This is similar to the univariate distribution, but with more than two
different kinds of individuals with different fitnesses.
***********************************************************************/

class MultiWalleniusNCHypergeometric {
    // This class encapsulates the different methods for calculating the
    // multivariate Wallenius noncentral hypergeometric probability function
public:
    MultiWalleniusNCHypergeometric(int32_t n, int32_t * m, double * odds, int colors, double accuracy = 1.E-8); // constructor
    void setParameters(int32_t n, int32_t * m, double * odds, int colors); // change parameters
    double probability(int32_t * x);                        // calculate probability function
    void mean(double * mu);                                 // calculate approximate mean

    // implementations of different calculation methods
protected:
    double binoexpand(void);                                // binomial expansion of integrand
    double laplace(void);                                   // Laplace's method with narrow integration interval
    double integrate(void);                                 // numerical integration

    // other subfunctions
    double lnbico(void);                                    // natural log of binomial coefficients
    void findpars(void);                                    // calculate r, w, E
    double integrate_step(double a, double b);              // used by integrate()
    double search_inflect(double t_from, double t_to);      // used by integrate()
    void reportError(const char* text);                     // report an error message

    // parameters
    double * omega;                                         // odds
    double accuracy;                                        // desired accuracy
    int32_t n;                                              // sample size
    int32_t N;                                              // total items in urn
    int32_t * m;                                            // items of each color in urn
    int32_t * x;                                            // items of each color sampled
    int32_t colors;                                         // number of different colors or types
    // parameters generated by findpars and used by probability, laplace, integrate:
    double r, rd, w, wr, E, phi2d;
    // generated by lnbico
    double bico;
};


/***********************************************************************
Class MultiWalleniusNCHypergeometricMoments
***********************************************************************/

class MultiWalleniusNCHypergeometricMoments : public MultiWalleniusNCHypergeometric {
    // This class calculates the exact mean and variance of the multivariate
    // Wallenius noncentral hypergeometric distribution by calculating all the 
    // possible x-combinations with probability < accuracy
public:
    MultiWalleniusNCHypergeometricMoments(int32_t n, int32_t * m, double * odds, int colors, double accuracy = 1.E-8)
        : MultiWalleniusNCHypergeometric(n, m, odds, colors, accuracy) {};
    double moments(double * mean, double * variance, int32_t * combinations = 0);

protected:
    // functions used internally
    double loop(int32_t n, int c);                          // recursive loops
    // data
    int32_t xi[MAXCOLORS];                                  // x vector to calculate probability of
    int32_t xm[MAXCOLORS];                                  // rounded approximate mean of x[i]
    int32_t remaining[MAXCOLORS];                           // number of balls of color > c in urn
    double sx[MAXCOLORS];                                   // sum of x*f(x)
    double sxx[MAXCOLORS];                                  // sum of x^2*f(x)
    int32_t sn;                                             // number of combinations
};


/***********************************************************************
Class FishersNCHypergeometric
***********************************************************************/

class FishersNCHypergeometric {
    // This class contains methods for calculating the univariate Fisher's
    // noncentral hypergeometric probability function
public:
    FishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy = 1E-8); // constructor
    double probability(int32_t x);                          // calculate probability function
    double mean(void);                                      // calculate approximate mean
    double variance(void);                                  // approximate variance
    int32_t mode(void);                                     // calculate mode (exact)
    double moments(double * mean, double * var);            // calculate exact mean and variance
protected:
    double lng(int32_t x);                                  // natural log of proportional function
    void reportError(const char* text);                     // report an error message

    // parameters
    double odds;                                            // odds ratio
    double logodds;                                         // ln odds ratio
    double accuracy;                                        // accuracy
    int32_t n, m, N;                                        // Parameters
    int32_t xmin, xmax;                                     // minimum and maximum of x
    // parameters used by subfunctions
    int32_t xLast;
    double mFac, xFac;                                      // log factorials
    double scale;                                           // scale to apply to lng function
    double rsum;                                            // reciprocal sum of proportional function
    bool parametersChanged;
};


/***********************************************************************
Class MultiFishersNCHypergeometric
This class contains functions for calculating the multivariate
Fisher's noncentral hypergeometric probability function and its mean and 
variance. Warning: the time consumption for first call to 
probability or moments is proportional to the total number of
possible x combinations, which may be extreme!
***********************************************************************/

class MultiFishersNCHypergeometric {
public:
    MultiFishersNCHypergeometric(int32_t n, int32_t * m, double * odds, int colors, double accuracy = 1E-9); // constructor
    double probability(int32_t * x);                       // calculate probability function
    void mean(double * mu);                                // calculate approximate mean
    void variance(double * var, double * mean = 0);        // calculate approximate variance and mean
    double moments(double * mean, double * variance, int32_t * combinations = 0); // calculate exact mean and variance

protected:
    void   mean1(double * mu);                             // calculate approximate mean, except for unused colors
    double lng(int32_t * x);                               // natural log of proportional function
    void   sumOfAll(void);                                 // calculates sum of proportional function for all x combinations
    double loop(int32_t n, int c);                         // recursive loops used by SumOfAll
    void reportError(const char* text);                    // report an error message

    int32_t n, N;                                          // copy of parameters
    int32_t Nu;                                            // number of balls in urn with nonzero weight
    int32_t m[MAXCOLORS];                                  // copy of all nonzero m
    double odds[MAXCOLORS];                                // copy of all nonzero odds
    double logodds[MAXCOLORS];                             // log odds
    int    colors;
    double mFac;                                           // sum of log m[i]!
    double scale;                                          // scale to apply to lng function
    double rsum;                                           // reciprocal sum of proportional function
    double accuracy;                                       // accuracy of calculation
    int    reduced;                                        // bit 0: some colors have m=0 or odds=0.
    // bit 1: all nonzero odds are equal
    int    usedcolors;                                     // number of colors with m > 0 and odds > 0
    int    nonzero[MAXCOLORS];                             // colors for which m and odds are not zero

    // data used by used by sumOfAll:
    int32_t xi[MAXCOLORS];                                 // x vector to calculate probability of
    int32_t xm[MAXCOLORS];                                 // rounded approximate mean of x[i]
    int32_t remaining[MAXCOLORS];                          // number of balls of color > c in urn
    double sx[MAXCOLORS];                                  // sum of x*f(x) or mean
    double sxx[MAXCOLORS];                                 // sum of x^2*f(x) or variance
    int32_t sn;                                            // number of possible combinations of x
};


/***********************************************************************
 Various common functions and constants
***********************************************************************/

double lnFac(int32_t n);                                   // Log factorial function
const int FAK_LEN = 1024;                                  // length of factorial table used by lnFac function
double fc_lnpk(int32_t k, int32_t L, int32_t m, int32_t n);// subfunction used by hypergeometric distribution
double pow2_1(double q, double * y0 = 0);                  // calculate 2^q and (1-2^q) without loss of precision.
double log1mx(double x, double x1);                        // Calculate log(1-x) without loss of precision when x is small.
double log1pow(double q, double x);                        // calculate log((1-e^q)^x) without loss of precision.
double fallingFactorial(double a, double b);               // calculates ln(a*(a-1)*(a-2)* ... * (a-b+1))
double lnFacr(double x);                                   // log factorial of non-integer x
int32_t floorLog2(float x);                                // This function calculates floor(log2(x)) for positive x.

const double ln2 = 0.693147180559945309417;                // log(2)

// The following tables are tables of residues of a certain expansion
// of the error function. These tables are used in the Laplace method
// for calculating Wallenius' noncentral hypergeometric distribution.
// There are ERFRES_N tables covering desired precisions from
// 2^(-ERFRES_B) to 2^(-ERFRES_E). Only the table that matches the
// desired precision is used. The tables are defined in erfres.h which
// is included in wnchyppr.cpp.

// constants for ErfRes tables:
const int ERFRES_B = 16;                                   // begin: -log2 of lowest precision
const int ERFRES_E = 40;                                   // end:   -log2 of highest precision
const int ERFRES_S =  2;                                   // step size from begin to end
const int ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1;       // number of tables
const int ERFRES_L = 48;                                   // length of each table

// tables of error function residues:
extern double erfRes [ERFRES_N][ERFRES_L];

// number of std. deviations to include in integral to obtain desired precision:
extern double numSDev[ERFRES_N];
