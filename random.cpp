/****************************  random.cpp   *****************************
* Author:        Agner Fog
* Date created:  2023-07-05
* Last modified: 2024-09-20
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines a pseudo random number generator and generation of
* random variates with the following distributions:
* - uniform, integer and float
* - bernoulli
* - normal, truncated normal
* - poisson
* - binomial
* - hypergeometric
* Wallenius' and Fisher's noncentral hypergeometric distributions are defined
* in a separate file wallenius.cpp
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"

/**************************************************************************
* Mersenne twister. Members of class CRandomMersenne
* Parameters are defined in file random.h
**************************************************************************/

void RandomMersenne::init0(int32_t seed) {
    // Seed generator
    const uint32_t factor = 1812433253U;
    mt[0] = seed;
    for (mti = 1; mti < MERS_N; mti++) {
        mt[mti] = (factor * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
    }
}

void RandomMersenne::randomInit(int32_t seed) {
    // Initialize and seed
    init0(seed);
    // Randomize some more
    for (int i = 0; i < 37; i++) bRandom();
}


void RandomMersenne::randomInitByArray(int32_t const seeds[], int NumSeeds) {
    // Seed by more than 32 bits
    int i, j, k;

    // Initialize
    init0(19650218);
    if (NumSeeds <= 0) return;

    // Randomize mt[] using whole seeds[] array
    i = 1;  j = 0;
    k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
        i++; j++;
        if (i >= MERS_N) { mt[0] = mt[MERS_N - 1]; i = 1; }
        if (j >= NumSeeds) j = 0;
    }
    for (k = MERS_N - 1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i;
        if (++i >= MERS_N) { mt[0] = mt[MERS_N - 1]; i = 1; }
    }
    mt[0] = 0x80000000U;  // MSB is 1; assuring non-zero initial array

    // Randomize some more
    mti = 0;
    for (int i = 0; i <= MERS_N; i++) bRandom();
}


uint32_t RandomMersenne::bRandom() {
    // Generate 32 random bits
    uint32_t y;

    if (mti >= MERS_N) {
        // Generate MERS_N words at one time
        const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;     // Lower MERS_R bits
        const uint32_t UPPER_MASK = 0xFFFFFFFFU << MERS_R;   // Upper (32 - MERS_R) bits (overflow intended)
        static const uint32_t mag01[2] = { 0, uint32_t(MERS_A) };

        int kk;
        for (kk = 0; kk < MERS_N - MERS_M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MERS_M] ^ (y >> 1) ^ mag01[y & 1];
        }

        for (; kk < MERS_N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MERS_M - MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
        }

        y = (mt[MERS_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MERS_N - 1] = mt[MERS_M - 1] ^ (y >> 1) ^ mag01[y & 1];
        mti = 0;
    }
    y = mt[mti++];
    // Tempering (May be omitted):
    y ^= y >> MERS_U;
    y ^= (y << MERS_S) & MERS_B;
    y ^= (y << MERS_T) & MERS_C;
    y ^= y >> MERS_L;
    return y;
}


double RandomMersenne::random() {
    // Output random float number in the interval 0 <= x < 1
    // Multiply by 2^(-32)
    return (double)bRandom() * (1. / (65536. * 65536.));
}


int32_t RandomMersenne::iRandom(int32_t min, int32_t max) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int32_t r = int32_t((double)(uint32_t)(max - min + 1) * random() + min);
    if (r > max) r = max;
    return r;
}


int32_t RandomMersenne::iRandomX(int32_t min, int32_t max) {
    // Output random integer in the interval min <= x <= max
    // Each output value has exactly the same probability.
    // This is obtained by rejecting certain bit values so that the number
    // of possible bit values is divisible by the interval length
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
    // 64 bit integers available. Use multiply and shift method
    uint32_t interval;                 // Length of interval
    uint64_t longran;                  // Random bits * interval
    uint32_t iran;                     // Longran / 2^32
    uint32_t remainder;                // Longran % 2^32

    interval = uint32_t(max - min + 1);
    if (interval != lastInterval) {
        // Interval length has changed. Must calculate rejection limit
        // Reject when remainder >= 2^32 / interval * interval
        // RLimit will be 0 if interval is a power of 2. No rejection then
        rLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
        lastInterval = interval;
    }
    do { // Rejection loop
        longran = (uint64_t)bRandom() * interval;
        iran = (uint32_t)(longran >> 32);
        remainder = (uint32_t)longran;
    } while (remainder > rLimit);
    // Convert back to signed and return result
    return (int32_t)iran + min;
}


/**************************************************************************
* Random Number generator of type 'Mother-Of-All generator'.
*
* This is a multiply-with-carry type of random number generator
* invented by George Marsaglia.  The algorithm is:
* S = 2111111111*X[n-4] + 1492*X[n-3] + 1776*X[n-2] + 5115*X[n-1] + C
* X[n] = S modulo 2^32
* C = floor(S / 2^32)
**************************************************************************/

// Output random bits
uint32_t RandomMother::bRandom() {
    uint64_t sum;
    sum = (uint64_t)2111111111UL * (uint64_t)x[3] +
        (uint64_t)1492 * (uint64_t)(x[2]) +
        (uint64_t)1776 * (uint64_t)(x[1]) +
        (uint64_t)5115 * (uint64_t)(x[0]) +
        (uint64_t)x[4];
    x[3] = x[2];  x[2] = x[1];  x[1] = x[0];
    x[4] = (uint32_t)(sum >> 32);                // Carry
    x[0] = (uint32_t)sum;                        // Low 32 bits of sum
    return x[0];
}


// returns a random number between 0 and 1:
double RandomMother::random() {
    return (double)bRandom() * (1. / (65536. * 65536.));
}


// returns integer random number in desired interval:
int32_t RandomMother::iRandom(int32_t min, int32_t max) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
    // Assume 64 bit integers supported. Use multiply and shift method
    uint32_t interval;                 // Length of interval
    uint64_t longran;                  // Random bits * interval
    uint32_t iran;                     // Longran / 2^32

    interval = (uint32_t)(max - min + 1);
    longran = (uint64_t)bRandom() * interval;
    iran = (uint32_t)(longran >> 32);
    // Convert back to signed and return result
    return int32_t(iran + min);
}


// this function initializes the random number generator:
void RandomMother::randomInit(int32_t seed) {
    int i;
    uint32_t s = seed;
    // make random numbers and put them into the buffer
    for (i = 0; i < 5; i++) {
        s = s * 29943829 - 1;
        x[i] = s;
    }
    // randomize some more
    for (i = 0; i < 19; i++) bRandom();
}


/**************************************************************************
* Combined random number generator
*
* Combines the Mersenne twister which gives a long cycle lengt with the
* multiply-with-carry generator which gives good diffusion or bifurcation.
*
* The combined generator is good for large, demanding applications that
* may simulate events with very low probability
*
**************************************************************************/

double RandomCombined::random() {
    // output random double precision float
    union Uqd {
        uint64_t q;
        double d;
    };
    Uqd u1;
    uint64_t r = bRandom();                      // get 64 random bits
    r |= (uint64_t)bRandom() << 32;
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2,
    // then subtract 1.0 to get the interval 0 <= x < 1.
    // A resolution of 2^-52 is sufficient. No need to put in an extra bit
    u1.q = (r >> 12) | 0x3FF0000000000000;       // bit 12 - 63
    return u1.d - 1.0;
}

float RandomCombined::randomf() {
    // output random single precision float
    union Uif {
        uint32_t i;
        float f;
    };
    Uif u1, u2;
    uint32_t r = bRandom();                      // get 32 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2
    // Subtract 1.0 if next bit is 0, or 1.0 - 2^-24 = 0.99999994f if next bit is 1
    u1.i = 0x3F800000 - ((r >> 8) & 1);          // bit 8
    u2.i = (r >> 9) | 0x3F800000;                // bit 9 - 31
    return u2.f - u1.f;
}


int32_t RandomCombined::iRandom(int32_t min, int32_t max) {
    // Output random integer in specified interval
    if (max <= min) {
        if (max == min) {
            return min;                          // Interval has only one value
        }
        else {
            return 0x80000000;                   // Error: interval length is negative
        }
    }
    uint32_t interval;                           // Length of interval
    uint64_t longran;                            // Random bits * interval
    uint32_t iran;                               // Longran / 2^32
    uint32_t ranbits = bRandom();                // Random bits
    interval = (uint32_t)(max - min) + 1u;       // Length of interval
    if (interval == 0) return ranbits;           // interval overflows
    longran = (uint64_t)ranbits * interval;
    iran = (uint32_t)(longran >> 32);
    return int32_t(iran + min);                  // Convert back to signed and return result
}

int32_t RandomCombined::iRandomX(int32_t min, int32_t max) {
    // Output random integer in the interval min <= x <= max
    // Each output value has exactly the same probability.
    // This is obtained by rejecting certain bit values so that the number
    // of possible bit values is divisible by the interval length
    if (max <= min) {
        if (max == min) {
            return min;                          // Interval has only one value
        }
        else {
            return 0x80000000;                   // Error: interval length is negative
        }
    }
    // Assume that 64 bit integers are available. Use multiply and shift method
    uint32_t interval;                           // Length of interval
    uint64_t longran;                            // Random bits * interval
    uint32_t iran;                               // longran / 2^32
    uint32_t remainder;                          // longran % 2^32

    interval = uint32_t(max - min) + 1u;         // length of interval
    if (interval == 0) return bRandom();         // avoid division by 0
    if (interval != randomixInterval) {
        // Interval length has changed. Must calculate rejection limit.
        // Reject when remainder >= 2^32 / interval * interval
        // randomixLimit will be 0 if interval is a power of 2. No rejection then
        randomixLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
        randomixInterval = interval;
    }
    do { // Rejection loop
        longran = (uint64_t)bRandom() * interval;
        iran = (uint32_t)(longran >> 32);
        remainder = (uint32_t)longran;
    } while (remainder > randomixLimit);
    return int32_t(iran + min);                  // Convert back to signed and return result
}

/*
const char * RandomCombined::test(int testcase) {
    // Test random number generator and print results
    const int N = 10000;                         // number of samples
    static char text[1024];                      // result string
    double r, sum = 0, ssum = 0, min, max;
    for (int i = 0; i < N; i++) {
        switch (testcase) {
        case 0: default:
            r = random();
            break;
        case 1:
            r = randomf();
            break;
        case 2:
            r = (int32_t)bRandom();
            break;
        case 3:
            r = iRandom(0, 1000);
            break;
        case 4:
            r = iRandomX(0, 1000000);
            break;
        }
        if (i == 0) {
            min = max = r;
        }
        else {
            if (r < min) min = r;
            if (r > max) max = r;
        }
        sum += r;
        ssum += r*r;
    }
    // output results
    sprintf_s(text, "Test case %i with %i samples:\n"
    "mean = %f, s.d. = %f\nmin = %f, max = %f",
        testcase, N, sum / N,
        sqrt((ssum - sum*sum/N) / (N-1)), min, max);
    return text;
}*/


/***********************************************************************
*  Non-uniform random variate generators
***********************************************************************/

// Constructor
RandomVariates::RandomVariates(int32_t seed) : RANDOM_BASE(seed) {
    normal_x2 = 0.; normal_x2_valid = false;
    pois_L_last = pois_f0 = 0.;
    pois_a = pois_h = pois_g = 0.; pois_bound = 0;
    bino_n_last = bino_mode = bino_bound = 0;
    bino_p_last = bino_a = bino_h = bino_g = bino_r1 = 0.;
    hyp_n_last = hyp_m_last = hyp_N_last = hyp_mode = hyp_mp = hyp_bound = 0;
    hyp_a = hyp_h = hyp_fm = 0.;
    wnc_n_last = wnc_m_last = wnc_N_last = wnc_bound1 = wnc_bound2 = wnc_mode = 0;
    wnc_o_last = wnc_a = wnc_h = wnc_k = 0;
    useChopDown = false;
    accuracy = 1.E-5;   // arbitrary accuracy
}

void RandomVariates::reportError(const char * text) {
    // report an error message via global object
    errors.reportError(text);
}

/*
const char * RandomVariates::test(int testcase) {
    // Test random number generator and print results
    if (testcase < 5) return RANDOM_BASE::test(testcase);  // test case for uniform generator

    const int N = 100000;                        // number of samples
    static char text[1024];                      // result string
    double r, sum = 0, ssum = 0, min, max;
    for (int i = 0; i < N; i++) {
        switch (testcase) {
        case 5: default:
            r = bernoulli(0.1);
            break;
        case 6:
            r = normal(100, 10);
            break;
        case 7:
            r = normalTrunc(100, 10, 15);
            break;
        case 8:
            r = poisson(100);
            break;
        case 9:
            r = binomial(100, 0.4);
            break;
        case 10:
            r = hypergeometric(10, 20, 100);
            //r = hypergeometric(80, 90, 200);
            break;
        case 11:
            r = walleniusNCHyp(10, 20, 100, 0.25);
            break;
        case 12: {// multi wallenius
            const int ncolors = 4;
            int32_t m[ncolors] = {200,300,400,200};
            double weights[ncolors] = {2,1,4,5};
            int n = 100;
            int32_t x[ncolors];
            multiWalleniusNCHyp(x, m, weights, n, ncolors);
            r = x[3];
            }
            break;
        }

        if (i == 0) {
            min = max = r;
        }
        else {
            if (r < min) min = r;
            if (r > max) max = r;
        }
        sum += r;
        ssum += r*r;
    }
    // output results
    sprintf_s(text, "Test case %i with %i samples:\n"
    "mean = %G, s.d. = %G\nmin = %G, max = %G",
        testcase, N, sum / N,
        sqrt((ssum - sum*sum/N) / (N-1)), min, max);
    return text;
}*/


bool RandomVariates::bernoulli(double p) {
    // Bernoulli distribution with parameter p. This function returns 
    // 0 or 1 with probability (1-p) and p, respectively.
    if (p < 0. || p > 1.) reportError("Parameter out of range in bernoulli function");
    return random() < p;
}

double RandomVariates::normal(double m, double s) {
    // Make normal distribution with mean m and standard deviation s
    // Use Box-Muller transform
    double normal_x1;                            // first random coordinate (normal_x2 is member of class)
    double w;                                    // radius
    if (normal_x2_valid) {                       // we have a valid result from last call
        normal_x2_valid = false;
        return normal_x2 * s + m;
    }
    // make two normally distributed variates by Box-Muller transformation
    do {
        normal_x1 = 2. * random() - 1.;
        normal_x2 = 2. * random() - 1.;
        w = normal_x1 * normal_x1 + normal_x2 * normal_x2;
    } while (w >= 1. || w < 1E-30);
    w = sqrt(log(w) * (-2. / w));
    normal_x1 *= w;  normal_x2 *= w;             // normal_x1 and normal_x2 are independent normally distributed variates
    normal_x2_valid = true;                      // save normal_x2 for next call
    return normal_x1 * s + m;                    // set mean and standard deviation
}

double RandomVariates::normalTrunc(double m, double s, double limit1, double limit2) {
    // Truncated normal distribution with mean m, standard deviation s.
    // The tails are cut off so that the output
    // is in the interval from limit1 to limit2
    if (limit1 > m || limit2 < m || limit1 >= limit2) {
        reportError("limits out of range in normalTrunc function");
    }
    double x;
    do {
        x = normal(m, s);
    } while (x < limit1 || x > limit2);          // reject if outside limits
    return x;
}

void RandomVariates::shuffle(int32_t * list, int32_t n1, int32_t n2) {
    // make a list of numbers from n1 to n2 in random order
    int32_t n = n2 - n1 + 1;  // size of list
    // first make ordered list
    for (int32_t i = 0; i < n; i++) list[i] = i + n1;
    // shuffle list
    for (int32_t j = 0; j < n - 1; j++) {
        int32_t k = iRandom(j, n - 1);
        // swap elements j and k
        int32_t temp = list[j];  list[j] = list[k];  list[k] = temp;
    }
}

int32_t RandomVariates::poisson(double L) {
    // This function generates a random variate with the poisson distribution.
    // Uses inversion by chop-down method for L < 17,
    // ratio-of-uniforms method for L >= 17.
    // For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.

    // Choose calculation method
    if (L < 17) {
        if (L < 1.E-6) {
            if (L == 0) return 0;
            if (L < 0) reportError("Parameter negative in poisson function");
            // For extremely small L we calculate the probabilities of x = 1
            // and x = 2 (ignoring higher x). The reason for using this 
            // method is to prevent numerical inaccuracies in other methods.
            return poissonLow(L);
        }
        else {
            // inversion method
            // The computation time for this method grows with L.
            // Gives overflow for L > 80
            return poissonInver(L);
        }
    }
    else {
        if (L > 2.E9) reportError("Parameter too big in poisson function");
        // ratio-of-uniforms method
        // The computation time for this method does not depend on L.
        // Use this where other methods would be slower.
        return poissonRatioUniforms(L);
    }
}

// Subfunctions used by poisson:
int32_t RandomVariates::poissonLow(double L) {
    // This subfunction generates a random variate with the poisson 
    // distribution for extremely low values of L.
    // The method is a simple calculation of the probabilities of x = 1
    // and x = 2. Higher values are ignored.
    // The reason for using this method is to avoid the numerical inaccuracies 
    // in other methods.
    double r, p0, p1, p2;
    r = random();
    p0 = 1. - L;                       // exp(-L) ~= 1-L
    p1 = p0 * L;
    p2 = p1 * L * 0.5;
    if (r < p2) return 2;
    else if (r < p1) return 1;
    else return 0;
}


int32_t RandomVariates::poissonInver(double L) {
    // This subfunction generates a random variate with the poisson
    // distribution using inversion by the chop down method (PIN).
    // Execution time grows with L. Gives overflow for L > 80.
    // The value of bound must be adjusted to the maximal value of L.
    const int bound = 130;             // safety bound. Must be > L + 8*sqrt(L).
    double r;                          // uniform random number
    double f;                          // function value
    int32_t x;                         // return value

    if (L != pois_L_last) {            // set up
        pois_L_last = L;
        pois_f0 = exp(-L);             // f(0) = probability of x=0
    }
    while (true) {
        r = random();  x = 0;  f = pois_f0;
        do {                           // recursive calculation: f(x) = f(x-1) * L / x
            r -= f;
            if (r <= 0.) return x;
            x++;
            f *= L;
            r *= x;                    // instead of f /= x
        } while (x <= bound);
        return x;
    }
}

// constants used by poissonRatioUniforms
const double SHAT1 = 2.943035529371538573;       // 8/e
const double SHAT2 = 0.8989161620588987408;      // 3-sqrt(12/e)

int32_t RandomVariates::poissonRatioUniforms(double L) {
    // This subfunction generates a random variate with the poisson 
    // distribution using the ratio-of-uniforms rejection method (PRUAt).

    // Execution time does not depend on L, except that it matters whether L
    // is within the range where ln(n!) is tabulated.

    // Reference: E. Stadlober: "The ratio of uniforms approach for generating
    // discrete random variates". Journal of Computational and Applied Mathematics,
    // vol. 31, no. 1, 1990, pp. 181-189.
    double u;                                    // uniform random
    double lf;                                   // ln(f(x))
    double x;                                    // real sample
    int32_t k;                                   // integer sample

    if (pois_L_last != L) {
        pois_L_last = L;                         // Set-up
        pois_a = L + 0.5;                        // hat center
        int32_t mode = (int32_t)L;               // mode
        pois_g = log(L);
        pois_f0 = mode * pois_g - lnFac(mode);   // value at mode
        pois_h = sqrt(SHAT1 * (L + 0.5)) + SHAT2;// hat width
        pois_bound = (int32_t)(pois_a + 6.0 * pois_h); // safety-bound
        if ((uint32_t)pois_bound > (uint32_t)INT32_MAX) pois_bound = INT32_MAX; // prevent overflow
    }
    while (true) {
        u = random();
        if (u == 0.) continue;                   // avoid division by 0
        x = pois_a + pois_h * (random() - 0.5) / u;
        if (x < 0 || x >= pois_bound) continue;  // reject if outside valid range
        k = int32_t(x);
        lf = k * pois_g - lnFac(k) - pois_f0;
        if (lf >= u * (4.0 - u) - 3.0) break;    // quick acceptance
        if (u * (u - lf) > 1.0) continue;        // quick rejection
        if (2.0 * log(u) <= lf) break;           // final acceptance
    }
    return k;
}

int32_t RandomVariates::binomial(int32_t n, double p) {
    // This function generates a random variate with the binomial distribution.
    // Uses inversion by chop-down method for n*p < 30, and ratio-of-uniforms method for n*p >= 30.
    // For n*p < 1.E-6 numerical inaccuracy is avoided by poisson approximation.
    int inv = 0;                                 // invert
    int32_t x;                                   // result
    double np = n * p;

    if (p > 0.5) {                               // faster calculation by inversion
        p = 1. - p;  inv = 1;
    }
    if (n <= 0 || p <= 0) {
        if (n == 0 || p == 0) {
            return inv * n;                      // only one possible result
        }
        // error exit
        reportError("Parameter out of range in binomial function");
    }
    // choose method
    if (np < 30.) {
        if (np < 1.E-6) {
            // Poisson approximation for extremely low np
            x = poissonLow(np);
        }
        else {
            // inversion method, using chop-down search from 0
            x = binomialInver(n, p);
        }
    }
    else {
        // ratio of uniforms method
        x = binomialRatioOfUniforms(n, p);
    }
    if (inv) {
        x = n - x;  // undo inversion
    }
    return x;
}

// Subfunctions used by binomial

int32_t RandomVariates::binomialInver(int32_t n, double p) {
    // Subfunction for Binomial distribution. Assumes p < 0.5.
    // Uses inversion method by search starting at 0.
    // Gives overflow for n*p > 60.
    // This method is fast when n*p is low. 
    double f0, f, q;
    double pn, r, rc;
    int32_t bound;
    int32_t x, n1, i;

    // f(0) = probability of x=0 is (1-p)^n
    // fast calculation of (1-p)^n
    f0 = 1.;  pn = 1. - p;  n1 = n;
    while (n1) {
        if (n1 & 1) f0 *= pn;
        pn *= pn;  n1 >>= 1;
    }
    // calculate safety bound
    rc = (n + 1) * p;
    bound = (int32_t)(rc + 11.0 * (sqrt(rc) + 1.0));
    if (bound > n) bound = n;
    q = p / (1. - p);

    while (true) {
        r = random();
        // recursive calculation: f(x) = f(x-1) * (n-x+1)/x*p/(1-p)
        f = f0;  x = 0;  i = n;
        do {
            r -= f;
            if (r <= 0) return x;
            x++;
            f *= q * i;
            r *= x;       // it is faster to multiply r by x than dividing f by x
            i--;
        } while (x <= bound);
    }
}


int32_t RandomVariates::binomialRatioOfUniforms(int32_t n, double p) {
    // Subfunction for Binomial distribution. Assumes p < 0.5.
    // Uses the Ratio-of-Uniforms rejection method.

    // The computation time hardly depends on the parameters, except that it matters
    // a lot whether parameters are within the range where the LnFac function is tabulated.

    // Reference: E. Stadlober: "The ratio of uniforms approach for generating
    // discrete random variates". Journal of Computational and Applied Mathematics,
    // vol. 31, no. 1, 1990, pp. 181-189.
    double u;                                    // uniform random
    double q1;                                   // 1-p
    double np;                                   // n*p
    double var;                                  // variance
    double lf;                                   // ln(f(x))
    double x;                                    // real sample
    int32_t k;                                   // integer sample

    if (bino_n_last != n || bino_p_last != p) {  // Set_up
        bino_n_last = n;
        bino_p_last = p;
        q1 = 1.0 - p;
        np = n * p;
        bino_mode = (int32_t)(np + p);           // mode
        bino_a = np + 0.5;                       // hat center
        bino_r1 = log(p / q1);
        bino_g = lnFac(bino_mode) + lnFac(n - bino_mode);
        var = np * q1;                           // variance
        bino_h = sqrt(SHAT1 * (var + 0.5)) + SHAT2;   // hat width
        bino_bound = (int32_t)(bino_a + 6.0 * bino_h);// safety-bound
        if (bino_bound > n) bino_bound = n;      // safety-bound
    }

    while (true) {                               // rejection loop
        u = random();
        if (u == 0.) continue;                   // avoid division by 0
        x = bino_a + bino_h * (random() - 0.5) / u;
        if (x < 0. || x > bino_bound) continue;  // reject, avoid overflow
        k = (int32_t)x;                          // truncate
        lf = (k - bino_mode) * bino_r1 + bino_g - lnFac(k) - lnFac(n - k); // ln(f(k))
        if (u * (4.0 - u) - 3.0 <= lf) break;    // lower squeeze accept
        if (u * (u - lf) > 1.0) continue;        // upper squeeze reject
        if (2.0 * log(u) <= lf) break;           // final acceptance
    }
    return k;
}

void RandomVariates::multinomial(int32_t * destination, double * p, int32_t n, int colors) {
    /*
    This function generates a vector of random variates, each with the binomial distribution.

    The multinomial distribution is the distribution you get when drawing
    balls from an urn with more than two colors, with replacement, with or without bias.

    Parameters:
    destination:   An output array to receive the number of balls of each color.
                   Must have space for at least 'colors' elements.
    p:             An input array containing the probability or fraction of each color
                   in the urn. Must have 'colors' elements.
                   All elements must be non-negative. The sum does not have to be 1
    n:             The number of balls drawn from the urn.
    colors:        The number of possible colors.
    */
    double s, sum;
    int32_t x;
    int i;
    if (n < 0 || colors < 0) reportError("Parameter negative in multinomial function");
    if (colors == 0) return;

    // compute sum of probabilities
    for (i = 0, sum = 0.; i < colors; i++) {
        s = p[i];
        if (s < 0) reportError("Parameter negative in multinomial function");
        sum += s;
    }
    if (sum == 0 && n > 0) reportError("Sum of p[] is zero in multinomial function");

    for (i = 0; i < colors - 1; i++) {
        // generate output by calling binomial (colors-1) times
        s = p[i];
        if (sum <= s) {
            // this fixes two problems:
            // 1. prevent division by 0 when sum = 0
            // 2. prevent s/sum getting bigger than 1 in case of rounding errors
            x = n;
        }
        else {
            x = binomial(n, s / sum);
        }
        n -= x; sum -= s;
        destination[i] = x;
    }
    // get the last one
    destination[i] = n;
}


int32_t RandomVariates::hypergeometric(int32_t n, int32_t m, int32_t N) {
    // This function generates a random variate with the hypergeometric
    // distribution. This is the distribution you get when drawing balls without 
    // replacement from an urn with two colors. n is the number of balls you take,
    // m is the number of red balls in the urn, N is the total number of balls in 
    // the urn, and the return value is the number of red balls you get.

    // This function uses inversion by chop-down search from the mode when
    // parameters are small, and the ratio-of-uniforms method when the former
    // method would be too slow or would give overflow.
    int32_t fak, addd;                 // used for undoing transformations
    int32_t x;                         // result

    // check if parameters are valid
    if (n > N || m > N || n < 0 || m < 0) {
        reportError("Parameter out of range in hypergeometric function");
    }
    // symmetry transformations
    fak = 1;  addd = 0;
    if (m > N / 2) {
        // invert m
        m = N - m;
        fak = -1;  addd = n;
    }
    if (n > N / 2) {
        // invert n
        n = N - n;
        addd += fak * m;  fak = -fak;
    }
    if (n > m) {
        // swap n and m
        x = n;  n = m;  m = x;
    }
    // cases with only one possible result end here
    if (n == 0)  return addd;

    // choose method
    if (N > 680 || n > 70) {
        // use ratio-of-uniforms method
        x = hypRatioOfUnifoms(n, m, N);
    }
    else {
        // inversion method, using chop-down search from mode
        x = hypInversionMod(n, m, N);
    }
    // undo symmetry transformations  
    return x * fak + addd;
}

// Subfunctions used by hypergeometric
int32_t RandomVariates::hypInversionMod(int32_t n, int32_t m, int32_t N) {
    // Subfunction for Hypergeometric distribution. Assumes 0 <= n <= m <= N/2.
    // Overflow protection is needed when N > 680 or n > 75.

    // Hypergeometric distribution by inversion method, using down-up 
    // search starting at the mode using the chop-down technique.
    // This method is faster than the rejection method when the variance is low.

    // Sampling 
    int32_t       I;                             // Loop counter
    int32_t       L = N - m - n;                 // Parameter
    double        modef;                         // mode, float
    double        Mp, np;                        // m + 1, n + 1
    double        p;                             // temporary
    double        U;                             // uniform random
    double        c, d;                          // factors in iteration
    double        divisor;                       // divisor, eliminated by scaling
    double        k1, k2;                        // float version of loop counter
    double        L1 = L;                        // float version of L

    Mp = (double)(m + 1);
    np = (double)(n + 1);

    if (N != hyp_N_last || m != hyp_m_last || n != hyp_n_last) {
        // set-up when parameters have changed
        hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;

        p = Mp / (N + 2.);
        modef = np * p;                          // mode, real
        hyp_mode = (int32_t)modef;               // mode, integer
        if (hyp_mode == modef && p == 0.5) {
            hyp_mp = hyp_mode--;
        }
        else {
            hyp_mp = hyp_mode + 1;
        }
        // mode probability, using log factorial function
        // (may read directly from fac_table if N < FAK_LEN)
        hyp_fm = exp(lnFac(N - m) - lnFac(L + hyp_mode) - lnFac(n - hyp_mode)
            + lnFac(m) - lnFac(m - hyp_mode) - lnFac(hyp_mode)
            - lnFac(N) + lnFac(N - n) + lnFac(n));

        // safety bound - guarantees at least 17 significant decimal digits
        // bound = min(n, (int32_t)(modef + k*c'))
        hyp_bound = (int32_t)(modef + 11. * sqrt(modef * (1. - p) * (1. - n / (double)N) + 1.));
        if (hyp_bound > n) hyp_bound = n;
    }

    // loop until accepted
    while (true) {
        U = random();                            // uniform random number to be converted

        // start chop-down search at mode
        if ((U -= hyp_fm) <= 0.) return(hyp_mode);
        c = d = hyp_fm;

        // alternating down- and upward search from the mode
        k1 = hyp_mp - 1;  k2 = hyp_mode + 1;
        for (I = 1; I <= hyp_mode; I++, k1--, k2++) {
            // Downward search from k1 = hyp_mp - 1
            divisor = (np - k1) * (Mp - k1);
            // Instead of dividing c with divisor, we multiply U and d because 
            // multiplication is faster. This will give overflow if N > 800
            U *= divisor;  d *= divisor;
            c *= k1 * (L1 + k1);
            if ((U -= c) <= 0.)  return(hyp_mp - I - 1); // = k1 - 1

            // Upward search from k2 = hyp_mode + 1
            divisor = k2 * (L1 + k2);
            // re-scale parameters to avoid time-consuming division
            U *= divisor;  c *= divisor;
            d *= (np - k2) * (Mp - k2);
            if ((U -= d) <= 0.)  return(hyp_mode + I);  // = k2
            // Values of n > 75 or N > 680 may give overflow if you leave out this..
            // overflow protection
            // if (U > 1.E100) {U *= 1.E-100; c *= 1.E-100; d *= 1.E-100;}
        }

        // Upward search from k2 = 2*mode + 1 to bound
        for (k2 = I = hyp_mp + hyp_mode; I <= hyp_bound; I++, k2++) {
            divisor = k2 * (L1 + k2);
            U *= divisor;
            d *= (np - k2) * (Mp - k2);
            if ((U -= d) <= 0.)  return(I);
            // more overflow protection
            // if (U > 1.E100) {U *= 1.E-100; d *= 1.E-100;}
        }
    }
}

int32_t RandomVariates::hypRatioOfUnifoms(int32_t n, int32_t m, int32_t N) {
    // Subfunction for Hypergeometric distribution using the ratio-of-uniforms rejection method.

    // This code is valid for 0 < n <= m <= N/2.

    // The computation time hardly depends on the parameters, except that it matters
    // a lot whether parameters are within the range where the LnFac function is tabulated.

    // Reference: E. Stadlober: "The ratio of uniforms approach for generating
    // discrete random variates". Journal of Computational and Applied Mathematics,
    // vol. 31, no. 1, 1990, pp. 181-189.
    int32_t L;                                             // N-m-n
    int32_t mode;                                          // mode
    int32_t k;                                             // integer sample
    double x;                                              // real sample
    double rNN;                                            // 1/(N*(N+2))
    double my;                                             // mean
    double var;                                            // variance
    double u;                                              // uniform random
    double lf;                                             // ln(f(x))

    L = N - m - n;
    if (hyp_N_last != N || hyp_m_last != m || hyp_n_last != n) {
        hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;  // Set-up
        rNN = 1. / ((double)N * (N + 2));                  // make two divisions in one
        my = (double)n * m * rNN * (N + 2);                // mean = n*m/N
        mode = (int32_t)(double(n + 1) * double(m + 1) * rNN * N); // mode = floor((n+1)*(m+1)/(N+2))
        var = (double)n * m * (N - m) * (N - n) / ((double)N * N * (N - 1)); // variance
        hyp_h = sqrt(SHAT1 * (var + 0.5)) + SHAT2;         // hat width
        hyp_a = my + 0.5;                                  // hat center
        hyp_fm = fc_lnpk(mode, L, m, n);                   // maximum
        hyp_bound = (int32_t)(hyp_a + 4.0 * hyp_h);        // safety-bound
        if (hyp_bound > n) hyp_bound = n;
    }
    while (true) {
        u = random();                                      // uniform random number
        if (u == 0.) continue;                             // avoid division by 0
        x = hyp_a + hyp_h * (random() - 0.5) / u;          // generate hat distribution
        if (x < 0. || x > 2.E9) continue;                  // reject, avoid overflow
        k = (int32_t)x;
        if (k > hyp_bound) continue;                       // reject if outside range
        lf = hyp_fm - fc_lnpk(k, L, m, n);                 // ln(f(k))
        if (u * (4.0 - u) - 3.0 <= lf) break;              // lower squeeze accept
        if (u * (u - lf) > 1.0) continue;                  // upper squeeze reject
        if (2.0 * log(u) <= lf) break;                     // final acceptance
    }
    return k;
}


void RandomVariates::multiHypergeometric(int32_t * destination, int32_t * source, int32_t n, int colors) {
    // Multivariate hypergeometric distribution
    int i;                                                 // loop counter
    int32_t m;                                             // number of items of one color
    int32_t N = 0;                                         // total number of items
    // check validity of parameters
    if (n < 0 || colors < 0) reportError("Parameter out of range in function multiHypergeometric");
    if (colors == 0) return;
    if (n == 0) {
        for (i = 0; i < colors; i++) destination[i] = 0; return;
    }

    // check validity of array parameters
    for (i = 0; i < colors; i++) {
        m = source[i];
        if (m < 0) reportError("Parameter negative in function multiHypergeometric");
        N += m;
    }
    if (n > N) reportError("Taking more items than there are in function multiHypergeometric");

    // sample by multiple calls to hypergeometric
    for (i = 0; i < colors - 1; i++) {
        m = source[i];
        int32_t x = hypergeometric(n, m, N);
        destination[i] = x;
        n -= x;
        N -= m;
    }
    destination[i] = n;
}


double lnFac(int32_t n) {
    // log factorial function. gives natural logarithm of n!

    // define constants
    static const double                // coefficients in Stirling approximation     
        C0 = 0.918938533204672722,     // ln(sqrt(2*pi))
        C1 = 1. / 12.,
        C3 = -1. / 360.;
    //  C5 =  1./1260.,                // use r^5 term if FAK_LEN < 50
    //  C7 = -1./1680.;                // use r^7 term if FAK_LEN < 20
    // static variables
    static double fac_table[FAK_LEN];  // table of ln(n!):
    static bool initialized = false;   // remember if fac_table has been initialized

    if (n < FAK_LEN) {
        if (n <= 1) {
            if (n < 0) {
                // report an error message via global object
                errors.reportError("parameter negative in function lnFac");
            }
            return 0;
        }
        if (!initialized) {            // first time. Must initialize table
            // make table of ln(n!)
            double sum = fac_table[0] = 0.;
            for (int i = 1; i < FAK_LEN; i++) {
                sum += log(double(i));
                fac_table[i] = sum;
            }
            initialized = 1;
        }
        return fac_table[n];
    }
    // not found in table. use Stirling approximation
    double  n1, r;
    n1 = n;  r = 1. / n1;
    return (n1 + 0.5) * log(n1) - n1 + C0 + r * (C1 + r * r * C3);
}

double fc_lnpk(int32_t k, int32_t L, int32_t m, int32_t n) {
    // subfunction used by hypergeometric and Fisher's noncentral hypergeometric distribution
    return(lnFac(k) + lnFac(m - k) + lnFac(n - k) + lnFac(L + k));
}
