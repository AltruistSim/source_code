/****************************  wallenius.cpp   *****************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2024-05-10
* Version:       3.002
* Project:       Altruist: Simulation of evolution in structured populations
* Description:
* This C++ file defines Wallenius' and Fisher's noncentral hypergeometric distributions
* including functions for calculating probabilities and generating
* random variates with these distributions
*
* Theory:
* To simulate positive selection, use Wallenius distribution.
* Example: The population of an animal species is limited by food resources.
* A mutant variant is better at finding food. The total number of survivors is
* fixed by the available food. The distribution of mutants after selection follows
* the Wallenius distribution. If there are more than two variants, use the
* multivariate distribution.
*
* To simulate negative selection, use the complementary Wallenius distribution.
* Example: The population is limited by predation. A mutant variant is better at
* avoiding predators.The total number of survivors is fixed by a predator-prey equilibrium.
* The distribution of mutants after selection follows the complementary Wallenius
* distribution. If there are more than two variants, use the multivariate distribution.
*
* The Wallenius and complementary Wallenius distributions both assume that the fates
* of individual animals are interdependent so that the total number of survivors is fixed.
* If the fates of each animal is independent of what happens to other conspecifics, then
* the distribution of survivors of each variant follows a binomial distribution.
*
* We may model the independent situation and yet fix the total number of survivors if
* the population is limited by some factor unrelated to the selection process in question,
* or simply for comparison with the positive and negative selection models. The
* Fisher's noncentral hypergeometric distribution can be used in this case. It models
* the distribution of two or more independent binomial variates, given their sum.
*
*
* The methods for calculaing probabilities in Wallenius' noncentral
* hypergeometric distribution is described in
* Fog, Agner:
* Calculation methods for Wallenius' noncentral hypergeometric distribution,
* Communications in Statistics - Simulation and Computation, 37(2), 258-273, 2008
* https://doi.org/10.1080/03610910701790269
* The methods for generating random variates from these distributions are described in
* Fog, Agner:
* Sampling methods for Wallenius' and Fisher's noncentral hypergeometric distributions,
* Communications in Statistics - Simulation and Computation, 37(2), 241-257, 2008
* https://doi.org/10.1080/03610910701790236
*
* (c) Copyright 2023-2024 Agner Fog.
* GNU General Public License, version 3.0 or later
******************************************************************************/

#include "stdafx.h"


/******************************************************************************
*  Non-uniform random variate generator for Wallenius distribution,
*  in class RandomVariates, univariate
******************************************************************************/

int32_t RandomVariates::walleniusNCHyp(int32_t n, int32_t m, int32_t N, double odds) {
    /*
    This function generates a random variate with Wallenius noncentral hypergeometric distribution.

    Wallenius noncentral hypergeometric distribution is the distribution you get when drawing balls
    without replacement from an urn containing red and white balls, with bias.

    We define the weight of the balls so that the probability of taking a
    particular ball is proportional to its weight. The value of odds is the
    normalized odds ratio: odds = weight(red) / weight(white).
    If all balls have the same weight, i.e. odds = 1, then we get the hypergeometric distribution.

    n is the number of balls you take,
    m is the number of red balls in the urn,
    N is the total number of balls in the urn,
    odds is the odds ratio,
    and the return value is the number of red balls you get.

    This function models positive selection with two phenotypes or genotypes

    Four different calculation methods are implemented. This function decides
    which method to use, based on the parameters.
    */

    // check parameters
    if (n >= N || m >= N || n <= 0 || m <= 0 || odds <= 0.) {
        // trivial cases
        if (n == 0 || m == 0) return 0;
        if (m == N) return n;
        if (n == N) return m;
        if (odds == 0.) {
            if (n > N - m) reportError("Not enough items with nonzero weight in function WalleniusNCHyp");
            return 0;
        }
        // illegal parameter    
        reportError("Parameter out of range in function WalleniusNCHyp");
    }

    if (odds == 1.) {
        // use hypergeometric function if odds == 1
        return hypergeometric(n, m, N);
    }

    if (n < 30) {
        return walleniusNCHypUrn(n, m, N, odds);
    }
    /* Table method useful only if called many times with the same parameters
    if (double(n) * N < 10000) {
        return walleniusNCHypTable(n, m, N, odds);
    } */

    return walleniusNCHypRatioOfUnifoms(n, m, N, odds);
    // the decision to use NoncentralHypergeometricInversion is
    // taken inside WalleniusNCHypRatioOfUnifoms based
    // on the calculated variance.
}

int32_t RandomVariates::complWalleniusNCHyp(int32_t n, int32_t m, int32_t N, double odds) {
    // Complementary Wallenius noncentral hypergeometric distribution.
    // This function models negative selection with two phenotypes or genotypes
    if (odds == 0.) return 0;  // avoid diding by 0
    return m - walleniusNCHypRatioOfUnifoms(N - n, m, N, 1. / odds);
}


// Subfunctions for WalleniusNCHyp:

int32_t RandomVariates::walleniusNCHypUrn(int32_t n, int32_t m, int32_t N, double odds) {
    // sampling from Wallenius noncentral hypergeometric distribution  by simulating urn model
    int32_t x;                         // sample
    int32_t m2;                        // items of color 2 in urn
    double mw1, mw2;                   // total weight of balls of color 1 or 2
    x = 0;  m2 = N - m;
    mw1 = m * odds;  mw2 = m2;
    do {
        if (random() * (mw1 + mw2) < mw1) {
            x++;  m--;
            if (m == 0) break;
            mw1 = m * odds;
        }
        else {
            m2--;
            if (m2 == 0) {
                x += n - 1; break;
            }
            mw2 = m2;
        }
    } while (--n);
    return x;
}

int32_t RandomVariates::walleniusNCHypRatioOfUnifoms(int32_t n, int32_t m, int32_t N, double odds) {
    // sampling from Wallenius noncentral hypergeometric distribution 
    // using ratio-of-uniforms rejection method.
    int32_t xmin, xmax;                // x limits
    double mean;                       // mean
    double variance;                   // variance
    double x;                          // real sample
    int32_t xi;                        // integer sample
    int32_t x2;                        // limit when searching for mode
    double u;                          // uniform random
    double f, f2;                      // probability function value
    double s123;                       // components 1,2,3 of hat width
    double s4;                         // component 4 of hat width
    double r1, r2;                     // temporaries
    const double rsqrt2pi = 0.3989422804014326857; // 1/sqrt(2*pi)

    // Make object for calculating mean and probability.
    WalleniusNCHypergeometric wnch(n, m, N, odds, accuracy);

    xmin = m + n - N; if (xmin < 0) xmin = 0;    // calculate limits
    xmax = n;         if (xmax > m) xmax = m;

    if (n != wnc_n_last || m != wnc_m_last || N != wnc_N_last || odds != wnc_o_last) {
        // set-up: This is done only when parameters have changed
        wnc_n_last = n;  wnc_m_last = m;  wnc_N_last = N;  wnc_o_last = odds;

        // find approximate mean
        mean = wnch.mean();

        // find approximate variance from Fisher's noncentral hypergeometric approximation
        r1 = mean * (m - mean); r2 = (n - mean) * (mean + N - n - m);
        variance = N * r1 * r2 / ((N - 1) * (m * r2 + (N - m) * r1));
        useChopDown = variance < 4.;             // use chop-down method if variance is low

        if (!useChopDown) {
            // find mode (same code in WalleniusNCHypergeometric::mode)
            wnc_mode = (int32_t)(mean);  f2 = 0.;
            if (odds < 1.) {
                if (wnc_mode < xmax) wnc_mode++;
                x2 = xmin;
                if (odds > 0.294 && N <= 10000000) {
                    x2 = wnc_mode - 1;
                } // search for mode can be limited
                for (xi = wnc_mode; xi >= x2; xi--) {
                    f = wnch.probability(xi);
                    if (f <= f2) break;
                    wnc_mode = xi; f2 = f;
                }
            }
            else {
                if (wnc_mode < xmin) wnc_mode++;
                x2 = xmax;
                if (odds < 3.4 && N <= 10000000) {
                    x2 = wnc_mode + 1;
                }  // search for mode can be limited
                for (xi = wnc_mode; xi <= x2; xi++) {
                    f = wnch.probability(xi);
                    if (f <= f2) break;
                    wnc_mode = xi; f2 = f;
                }
            }
            wnc_k = f2;                          // value at mode

            // find approximate variance from normal distribution approximation
            variance = rsqrt2pi / wnc_k;  variance *= variance;

            // find center and width of hat function
            wnc_a = mean + 0.5;
            s123 = 0.40 + 0.8579 * sqrt(variance + 0.5) + 0.4 * fabs(mean - wnc_mode);
            s4 = 0.;
            r1 = xmax - mean - s123;  r2 = mean - s123 - xmin;
            if (r1 > r2) r1 = r2;
            if ((odds > 5. || odds < 0.2) && r1 >= -0.5 && r1 <= 8.) {
                // s4 correction needed
                if (r1 < 1.) r1 = 1.;
                s4 = 0.029 * pow(double(N), 0.23) / (r1 * r1);
            }
            wnc_h = 2. * (s123 + s4);

            // find safety bounds
            wnc_bound1 = (int32_t)(mean - 4. * wnc_h);
            if (wnc_bound1 < xmin) wnc_bound1 = xmin;
            wnc_bound2 = (int32_t)(mean + 4. * wnc_h);
            if (wnc_bound2 > xmax) wnc_bound2 = xmax;
        }
    }

    if (useChopDown) { // for small variance, use chop down inversion
        return walleniusNCHypInversion(n, m, N, odds);
    }

    // use ratio-of-uniforms rejection method
    while (true) {                               // rejection loop
        u = random();
        if (u == 0.) continue;                   // avoid division by 0
        x = wnc_a + wnc_h * (random() - 0.5) / u;
        if (x < 0. || x > 2E9) continue;         // reject, avoid overflow
        xi = int32_t(x);                       // truncate
        if (xi < wnc_bound1 || xi > wnc_bound2) {
            continue;                            // reject if outside safety bounds
        }
#if false // use rejection in x-domain
        if (xi == wnc_mode) break;               // accept      
        f = wnch.probability(xi);                // function value
        if (f > wnc_k * u * u) {
            break;                               // acceptance
        }
#else // use rejection in t-domain (this is faster)
        double hx, s2, xma2;                     // compute h(x)
        s2 = wnc_h * 0.5;  s2 *= s2;
        xma2 = xi - (wnc_a - 0.5);
        xma2 *= xma2;
        hx = (s2 >= xma2) ? 1. : s2 / xma2;
        // rejection in t-domain implemented in WalleniusNCHypergeometric::bernouilliH
        if (wnch.bernouilliH(xi, hx * wnc_k * 1.01, u * u * wnc_k * 1.01, this)) {
            break;                               // acceptance
        }
#endif      
    }                                            // rejection
    return xi;
}


int32_t RandomVariates::walleniusNCHypInversion(int32_t n, int32_t m, int32_t N, double odds) {
    // sampling from Wallenius noncentral hypergeometric distribution 
    // using down-up search starting at the mean using the chop-down technique.
    // This method is faster than the rejection method when the variance is low.
    int32_t x1s, x1, x2;                     // search values
    int32_t xmin, xmax;                          // x limits
    double   u;                                  // uniform random number to be converted
    double   f;                                  // probability function value
    double   accura;                             // absolute accuracy
    int      updown;                             // 1 = search down, 2 = search up, 3 = both

    // Make objects for calculating mean and probability.
    // It is more efficient to have two identical objects, one for down search
    // and one for up search, because they are obtimized for consecutive x values.
    WalleniusNCHypergeometric wnch1(n, m, N, odds, accuracy);
    WalleniusNCHypergeometric wnch2(n, m, N, odds, accuracy);

    accura = accuracy * 0.01;
    if (accura > 1E-7) accura = 1E-7;            // absolute accuracy

    x1s = (int32_t)(wnch1.mean());               // start at floor x1 and ceiling x2 of mean
    xmin = m + n - N; if (xmin < 0) xmin = 0;    // calculate limits
    xmax = n;     if (xmax > m) xmax = m;

    while (true) {                               // loop until accepted (normally executes only once)
        x1 = x1s;
        x2 = x1s + 1;
        updown = 3;                              // start searching both up and down
        u = random();                            // uniform random number to be converted
        while (updown) {                         // search loop
            if (updown & 1) {                    // search down
                if (x1 < xmin) {
                    updown &= ~1;                // stop searching down
                }
                else {
                    f = wnch1.probability(x1);
                    u -= f;                      // subtract probability until 0
                    if (u <= 0.) return x1;
                    x1--;
                    if (f < accura) updown &= ~1;// stop searching down
                }
            }
            if (updown & 2) {                    // search up
                if (x2 > xmax) {
                    updown &= ~2;                // stop searching up
                }
                else {
                    f = wnch2.probability(x2);
                    u -= f;                      // subtract probability until 0
                    if (u <= 0.) return x2;
                    x2++;
                    if (f < accura) updown &= ~2;// stop searching down
                }
            }
        }
    }
}


/***********************************************************************
* Members of class WalleniusNCHypergeometric for calculating
* probabilities in Wallenius' noncentral hypergeometric distribution
***********************************************************************/

// constructor
WalleniusNCHypergeometric::WalleniusNCHypergeometric(int32_t n_, int32_t m_, int32_t N_, double odds_, double accuracy_) {
    accuracy = accuracy_;
    setParameters(n_, m_, N_, odds_);
}


void WalleniusNCHypergeometric::setParameters(int32_t n_, int32_t m_, int32_t N_, double odds) {
    // change parameters
    if (n_ < 0 || n_ > N_ || m_ < 0 || m_ > N_ || odds < 0) {
        reportError("Parameter out of range in WalleniusNCHypergeometric");
    }
    n = n_; m = m_; N = N_; omega = odds;        // set parameters
    xmin = m + n - N;  if (xmin < 0) xmin = 0;   // calculate xmin
    xmax = n;  if (xmax > m) xmax = m;           // calculate xmax
    xLastBico = xLastFindpars = -99;             // indicate last x is invalid
    r = 1.;                                      // initialize
}


double WalleniusNCHypergeometric::mean(void) {
    // find approximate mean for Wallenius' distribution
    int iter;                                    // number of iterations
    double a, b;                                 // temporaries in calculation of first guess
    double mean, mean1;                          // iteration value of mean
    double m1r, m2r;                             // 1/m, 1/m2
    double e1, e2;                               // temporaries
    double g;                                    // function to find root of
    double gd;                                   // derivative of g
    double omegar;                               // 1/omega

    if (omega == 1.) { // simple hypergeometric
        return (double)(m)*n / N;
    }

    if (omega == 0.) {
        if (n > N - m) reportError("Not enough items with nonzero weight in WalleniusNCHypergeometric::mean");
        return 0.;
    }

    if (xmin == xmax) return xmin;

    // calculate Cornfield mean of Fisher noncentral hypergeometric distribution as first guess
    a = (m + n) * omega + (N - m - n);
    b = a * a - 4. * omega * (omega - 1.) * m * n;
    b = b > 0. ? sqrt(b) : 0.;
    mean = (a - b) / (2. * (omega - 1.));
    if (mean < xmin) mean = xmin;
    if (mean > xmax) mean = xmax;

    m1r = 1. / m;  m2r = 1. / (N - m);
    iter = 0;

    if (omega > 1.) {
        do { // Newton Raphson iteration
            mean1 = mean;
            e1 = 1. - (n - mean) * m2r;
            if (e1 < 1E-14) {
                e2 = 0.;     // avoid underflow
            }
            else {
                e2 = pow(e1, omega - 1.);
            }
            g = e2 * e1 + (mean - m) * m1r;
            gd = e2 * omega * m2r + m1r;
            mean -= g / gd;
            if (mean < xmin) mean = xmin;
            if (mean > xmax) mean = xmax;
            if (++iter > 40) {
                reportError("Search for mean failed in function WalleniusNCHypergeometric::mean");
            }
        } while (fabs(mean1 - mean) > 2E-6);
    }
    else { // omega < 1
        omegar = 1. / omega;
        do { // Newton Raphson iteration
            mean1 = mean;
            e1 = 1. - mean * m1r;
            if (e1 < 1E-14) {
                e2 = 0.;     // avoid underflow
            }
            else {
                e2 = pow(e1, omegar - 1.);
            }
            g = 1. - (n - mean) * m2r - e2 * e1;
            gd = e2 * omegar * m1r + m2r;
            mean -= g / gd;
            if (mean < xmin) mean = xmin;
            if (mean > xmax) mean = xmax;
            if (++iter > 40) {
                reportError("Search for mean failed in function WalleniusNCHypergeometric::mean");
            }
        } while (fabs(mean1 - mean) > 2E-6);
    }
    return mean;
}

double WalleniusNCHypergeometric::variance(void) {
    // find approximate variance (poor approximation)    
    double my = mean(); // approximate mean
    // find approximate variance from Fisher's noncentral hypergeometric approximation
    double r1 = my * (m - my); double r2 = (n - my) * (my + N - n - m);
    if (r1 <= 0. || r2 <= 0.) return 0.;
    double var = N * r1 * r2 / ((N - 1) * (m * r2 + (N - m) * r1));
    if (var < 0.) var = 0.;
    return var;
}

int32_t WalleniusNCHypergeometric::mode(void) {
    // find mode
    int32_t mode;

    if (omega == 1.) {
        // simple hypergeometric
        int32_t L = m + n - N;
        int32_t m1 = m + 1, n1 = n + 1;
        mode = int32_t((double)m1 * n1 * omega / ((m1 + n1) * omega - L));
    }
    else {
        // find mode
        double f, f2 = 0.; // f2 = -1.; 
        int32_t xi, x2;
        int32_t xmin = m + n - N;  if (xmin < 0) xmin = 0; // calculate xmin
        int32_t xmax = n;  if (xmax > m) xmax = m;         // calculate xmax

        mode = (int32_t)mean();                            // floor(mean)
        if (omega < 1.) {
            if (mode < xmax) mode++;                       // ceil(mean)
            x2 = xmin;                                     // lower limit
            if (omega > 0.294 && N <= 10000000) {
                x2 = mode - 1;
            }                                              // search for mode can be limited
            for (xi = mode; xi >= x2; xi--) {
                f = probability(xi);
                if (f <= f2) break;
                mode = xi;  f2 = f;
            }
        }
        else {
            if (mode < xmin) mode++;
            x2 = xmax;                                     // upper limit
            if (omega < 3.4 && N <= 10000000) {
                x2 = mode + 1;
            }                                              // search for mode can be limited
            for (xi = mode; xi <= x2; xi++) {
                f = probability(xi);
                if (f <= f2) break;
                mode = xi; f2 = f;
            }
        }
    }
    return mode;
}

double WalleniusNCHypergeometric::probability(int32_t x_) {
    // calculate probability function. choosing best method
    x = x_;
    if (x < xmin || x > xmax) return 0.;
    if (xmin == xmax) return 1.;

    if (omega == 1.) {       // hypergeometric
        return exp(lnbico() + lnFac(n) + lnFac(N - n) - lnFac(N));
    }

    if (omega == 0.) {
        if (n > N - m) reportError("Not enough items with nonzero weight in WalleniusNCHypergeometric::probability");
        return x == 0;
    }

    int32_t x2 = n - x;
    int32_t x0 = x < x2 ? x : x2;
    int em = (x == m || x2 == N - m);

    if (x0 == 0 && n > 500) {
        return binoexpand();
    }

    if (double(n) * x0 < 1000 || (double(n) * x0 < 10000 && (N > 1000. * n || em))) {
        return recursive();
    }

    if (x0 <= 1 && N - n <= 1) {
        return binoexpand();
    }

    findpars();

    if (w < 0.04 && E < 10 && (!em || w > 0.004)) {
        return laplace();
    }

    return integrate();
}

void WalleniusNCHypergeometric::findpars() {
    // calculate d, E, r, w
    if (x == xLastFindpars) {
        return;              // all values are unchanged since last call
    }

    // find r to center peak of integrand at 0.5
    double dd, d1, z, zd, rr, lastr, rrc, rt, r2, r21, a, b, dummy;
    double oo[2];
    double xx[2] = { double(x), double(n - x) };
    int i, j = 0;
    if (omega > 1.) {        // make both omegas <= 1 to avoid overflow
        oo[0] = 1.;  oo[1] = 1. / omega;
    }
    else {
        oo[0] = omega;  oo[1] = 1.;
    }
    dd = oo[0] * (m - x) + oo[1] * (N - m - xx[1]);
    d1 = 1. / dd;
    E = (oo[0] * m + oo[1] * (N - m)) * d1;
    rr = r;
    if (rr <= d1) rr = 1.2 * d1;                 // initial guess
    // Newton-Raphson iteration to find r
    do {
        lastr = rr;
        rrc = 1. / rr;
        z = dd - rrc;
        zd = rrc * rrc;
        for (i = 0; i < 2; i++) {
            rt = rr * oo[i];
            if (rt < 100.) {                     // avoid overflow if rt big
                r21 = pow2_1(rt, &r2);           // r2=2^r, r21=1.-2^r
                a = oo[i] / r21;                 // omegai/(1.-2^r)
                b = xx[i] * a;                   // x*omegai/(1.-2^r)
                z += b;
                zd += b * a * ln2 * r2;
            }
        }
        if (zd == 0) reportError("can't find r in function WalleniusNCHypergeometric::findpars");
        rr -= z / zd;
        if (rr <= d1) rr = lastr * 0.125 + d1 * 0.875;
        if (++j == 70) reportError("convergence problem searching for r in function WalleniusNCHypergeometric::findpars");
    } while (fabs(rr - lastr) > rr * 1.E-6);
    if (omega > 1) {
        dd *= omega;  rr *= oo[1];
    }
    r = rr;  rd = rr * dd;

    // find peak width
    double ro, k1, k2;
    ro = r * omega;
    if (ro < 300) {                              // avoid overflow
        k1 = pow2_1(ro, &dummy);
        k1 = -1. / k1;
        k1 = omega * omega * (k1 + k1 * k1);
    }
    else k1 = 0.;
    if (r < 300) {                               // avoid overflow
        k2 = pow2_1(r, &dummy);
        k2 = -1. / k2;
        k2 = (k2 + k2 * k2);
    }
    else k2 = 0.;
    phi2d = -4. * r * r * (x * k1 + (n - x) * k2);
    if (phi2d >= 0.) {
        reportError("peak width undefined in function WalleniusNCHypergeometric::findpars");
        /* wr = r = 0.; */
    }
    else {
        wr = sqrt(-phi2d); w = 1. / wr;
    }
    xLastFindpars = x;
}

double WalleniusNCHypergeometric::binoexpand() {
    // calculate by binomial expansion of integrand
    // only for x < 2 or n-x < 2 (not implemented for higher x because of loss of precision)
    int32_t x1, m1, m2;
    double o;
    if (x > n / 2) { // invert
        x1 = n - x; m1 = N - m; m2 = m; o = 1. / omega;
    }
    else {
        x1 = x; m1 = m; m2 = N - m; o = omega;
    }
    if (x1 == 0) {
        return exp(fallingFactorial(m2, n) - fallingFactorial(m2 + o * m1, n));
    }
    if (x1 == 1) {
        double d, e, q, q0, q1;
        q = fallingFactorial(m2, n - 1);
        e = o * m1 + m2;
        q1 = q - fallingFactorial(e, n);
        e -= o;
        q0 = q - fallingFactorial(e, n);
        d = e - (n - 1);
        return m1 * d * (exp(q0) - exp(q1));
    }
    reportError("x > 1 not supported by function WalleniusNCHypergeometric::binoexpand");
    return 0;
}

double WalleniusNCHypergeometric::recursive() {
    // recursive calculation
    // Wallenius noncentral hypergeometric distribution by recursion formula
    // Approximate by ignoring probabilities < accuracy and minimize storage requirement
    const int BUFSIZE = 512;           // buffer size
    double p[BUFSIZE + 2];             // probabilities
    double * p1, * p2;                 // offset into p
    double mxo;                        // (m-x)*omega
    double Nmnx;                       // N-m-nu+x
    double y, y1;                      // save old p[x] before it is overwritten
    double d1, d2, dcom;               // divisors in probability formula
    double accuracya;                  // absolute accuracy
    int32_t xi, nu;                    // xi, nu = recursion values of x, n
    int32_t x1, x2;                    // xi_min, xi_max

    accuracya = 0.005f * accuracy;     // absolute accuracy
    p1 = p2 = p + 1;                   // make space for p1[-1]
    p1[-1] = 0.;  p1[0] = 1.;          // initialize for recursion
    x1 = x2 = 0;
    for (nu = 1; nu <= n; nu++) {
        if (n - nu < x - x1 || p1[x1] < accuracya) {
            x1++;                      // increase lower limit when breakpoint passed or probability negligible
            p2--;                      // compensate buffer offset in order to reduce storage space
        }
        if (x2 < x && p1[x2] >= accuracya) {
            x2++;  y1 = 0.;            // increase upper limit until x has been reached
        }
        else {
            y1 = p1[x2];
        }
        if (x1 > x2) return 0.;
        if (p2 + x2 - p > BUFSIZE) reportError("buffer overrun in function WalleniusNCHypergeometric::recursive");

        mxo = (m - x2) * omega;
        Nmnx = N - m - nu + x2 + 1;
        for (xi = x2; xi >= x1; xi--) {// backwards loop
            d2 = mxo + Nmnx;
            mxo += omega; Nmnx--;
            d1 = mxo + Nmnx;
            dcom = 1. / (d1 * d2);     // save a division by making common divisor
            y = p1[xi - 1] * mxo * d2 * dcom + y1 * (Nmnx + 1) * d1 * dcom;
            y1 = p1[xi - 1];           // (warning: pointer alias, can't swap instruction order)
            p2[xi] = y;
        }
        p1 = p2;
    }
    if (x < x1 || x > x2) return 0.;
    return p1[x];
}

double WalleniusNCHypergeometric::laplace() {
    // Laplace's method with narrow integration interval, 
    // using error function residues table, defined below
    // Note that this function can only be used when the integrand peak is narrow.
    // findpars() must be called before this function.

    const int COLORS = 2;              // number of colors
    const int MAXDEG = 40;             // arraysize, maximum expansion degree
    int degree;                        // max expansion degree
    double accur;                      // stop expansion when terms below this threshold
    double omegai[COLORS] = { omega, 1.0 }; // weights for each color
    double xi[COLORS] = { double(x), double(n - x) }; // number of each color sampled
    double f0;                         // factor outside integral
    double rho[COLORS];                // r*omegai
    double qi;                         // 2^(-rho)
    double qi1;                        // 1-qi
    double qq[COLORS];                 // qi / qi1
    double eta[COLORS + 1][MAXDEG + 1];// eta coefficients
    double phideri[MAXDEG + 1];        // derivatives of phi
    double PSIderi[MAXDEG + 1];        // derivatives of PSI
    double * erfresp;                  // pointer to table of error function residues

    // variables in asymptotic summation
    static const double sqrt8 = 2.828427124746190098; // sqrt(8)
    double qqpow;                      // qq^j
    double pow2k;                      // 2^k
    double bino;                       // binomial coefficient  
    double vr;                         // 1/v, v = integration interval
    double v2m2;                       // (2*v)^(-2)
    double v2mk1;                      // (2*v)^(-k-1)
    double s;                          // summation term
    double sum;                        // Taylor sum

    int i;                             // loop counter for color
    int j;                             // loop counter for derivative
    int k;                             // loop counter for expansion degree
    int ll;                            // k/2
    int converg = 0;                   // number of consequtive terms below accuracy
    int precisionIndex;                // index into erfRes table according to desired precision

    // initialize
    for (k = 0; k <= 2; k++)  phideri[k] = PSIderi[k] = 0;

    // find rho[i], qq[i], first eta coefficients, and zero'th derivative of phi
    for (i = 0; i < COLORS; i++) {
        rho[i] = r * omegai[i];
        if (rho[i] > 40.) {
            qi = 0.;  qi1 = 1.;
        }                    // avoid underflow
        else {
            qi1 = pow2_1(-rho[i], &qi);
        }                   // qi=2^(-rho), qi1=1.-2^(-rho)
        qq[i] = qi / qi1;              // 2^(-r*omegai)/(1.-2^(-r*omegai))
        // peak = zero'th derivative
        phideri[0] += xi[i] * log1mx(qi, qi1);
        // eta coefficients
        eta[i][0] = 0.;
        eta[i][1] = eta[i][2] = rho[i] * rho[i];
    }

    // r, rd, and w must be calculated by findpars()
    // zero'th derivative
    phideri[0] -= (rd - 1.) * ln2;
    // scaled factor outside integral
    f0 = rd * exp(phideri[0] + lnbico());

    vr = sqrt8 * w;
    phideri[2] = phi2d;

    // get table according to desired precision
    precisionIndex = (-floorLog2((float)accuracy) - ERFRES_B + ERFRES_S - 1) / ERFRES_S;
    if (precisionIndex < 0) precisionIndex = 0;
    if (precisionIndex > ERFRES_N - 1) precisionIndex = ERFRES_N - 1;
    while (w * numSDev[precisionIndex] > 0.3) {
        // check if integration interval is too wide
        if (precisionIndex == 0) {
            reportError("Laplace method failed. Peak width too high in function WalleniusNCHypergeometric::laplace");
            break;
        }
        precisionIndex--;                        // reduce precision to keep integration interval narrow
    }
    erfresp = erfRes[precisionIndex];            // choose desired table

    degree = MAXDEG;                             // max expansion degree
    if (degree >= ERFRES_L * 2) degree = ERFRES_L * 2 - 2;

    // set up for starting loop at k=3
    v2m2 = 0.25 * vr * vr;                       // (2*v)^(-2)
    PSIderi[0] = 1.;
    pow2k = 8.;
    sum = 0.5 * vr * erfresp[0];
    v2mk1 = 0.5 * vr * v2m2 * v2m2;
    accur = accuracy * sum;

    // summation loop
    for (k = 3; k <= degree; k++) {
        phideri[k] = 0.;

        // loop for all (2) colors
        for (i = 0; i < COLORS; i++) {
            eta[i][k] = 0.;
            // backward loop for all powers
            for (j = k; j > 0; j--) {
                // find coefficients recursively from previous coefficients
                eta[i][j] = eta[i][j] * (j * rho[i] - (k - 2)) + eta[i][j - 1] * rho[i] * (j - 1);
            }
            qqpow = 1.;
            // forward loop for all powers
            for (j = 1; j <= k; j++) {
                qqpow *= qq[i];                  // qq^j
                // contribution to derivative
                phideri[k] += xi[i] * eta[i][j] * qqpow;
            }
        }

        // finish calculation of derivatives
        phideri[k] = -pow2k * phideri[k] + 2 * (1 - k) * phideri[k - 1];

        pow2k *= 2.;    // 2^k

        // loop to calculate derivatives of PSI from derivatives of psi.
        // terms # 0, 1, 2, k-2, and k-1 are zero and not included in loop.
        // The j'th derivatives of psi are identical to the derivatives of phi for j>2, and
        // zero for j=1,2. Hence we are using phideri[j] for j>2 here.
        PSIderi[k] = phideri[k];                 // this is term # k
        bino = 0.5 * (k - 1) * (k - 2);          // binomial coefficient for term # 3
        for (j = 3; j < k - 2; j++) {            // loop for remaining nonzero terms (if k>5)
            PSIderi[k] += PSIderi[k - j] * phideri[j] * bino;
            bino *= double(k - j) / double(j);
        }
        if ((k & 1) == 0) {                      // only for even k
            ll = k / 2;
            s = PSIderi[k] * v2mk1 * erfresp[ll];
            sum += s;

            // check for convergence of Taylor expansion
            if (fabs(s) < accur) converg++; else converg = 0;
            if (converg > 1) break;

            // update recursive expressions
            v2mk1 *= v2m2;
        }
    }
    // multiply by terms outside integral  
    return f0 * sum;
}

double WalleniusNCHypergeometric::integrate() {
    // Wallenius non-central hypergeometric distribution function
    // calculation by numerical integration with variable-length steps
    // Note: findpars() must be called before this function.
    double s;                                    // result of integration step
    double sum;                                  // integral
    double ta, tb;                               // subinterval for integration step

    lnbico();                                    // compute log of binomial coefficients

    // choose method:
    if (w < 0.02 || (w < 0.1 && (x == m || n - x == N - m) && accuracy > 1E-6)) {
        // normal method. Step length determined by peak width w
        double delta, s1;
        s1 = accuracy < 1E-9 ? 0.5 : 1.;
        delta = s1 * w;                          // integration steplength
        ta = 0.5 + 0.5 * delta;
        sum = integrate_step(1. - ta, ta);       // first integration step around center peak
        do {
            tb = ta + delta;
            if (tb > 1.) tb = 1.;
            s = integrate_step(ta, tb);          // integration step to the right of peak
            s += integrate_step(1. - tb, 1. - ta); // integration step to the left of peak
            sum += s;
            if (s < accuracy * sum) break;       // stop before interval finished if accuracy reached
            ta = tb;
            if (tb > 0.5 + w) delta *= 2.;       // increase step length far from peak
        } while (tb < 1.);
    }
    else {
        // difficult situation. Step length determined by inflection points
        double t1, t2, tinf, delta, delta1;
        sum = 0.;
        // do left and right half of integration interval separately:
        for (t1 = 0., t2 = 0.5; t1 < 1.; t1 += 0.5, t2 += 0.5) {
            // integrate from 0 to 0.5 or from 0.5 to 1
            tinf = search_inflect(t1, t2);       // find inflection point
            delta = tinf - t1; if (delta > t2 - tinf) delta = t2 - tinf; // distance to nearest endpoint
            delta *= 1. / 7.;                    // 1/7 will give 3 steps to nearest endpoint
            if (delta < 1E-4) delta = 1E-4;
            delta1 = delta;
            // integrate from tinf forwards to t2
            ta = tinf;
            do {
                tb = ta + delta1;
                if (tb > t2 - 0.25 * delta1) tb = t2;      // last step of this subinterval
                s = integrate_step(ta, tb);                // integration step
                sum += s;
                delta1 *= 2;                               // double steplength
                if (s < sum * 1E-4) delta1 *= 8.;          // large step when s small
                ta = tb;
            } while (tb < t2);
            if (tinf) {
                // integrate from tinf backwards to t1
                tb = tinf;
                do {
                    ta = tb - delta;
                    if (ta < t1 + 0.25 * delta) ta = t1;   // last step of this subinterval
                    s = integrate_step(ta, tb);            // integration step
                    sum += s;
                    delta *= 2;                            // double steplength
                    if (s < sum * 1E-4) delta *= 8.;       // large step when s small
                    tb = ta;
                } while (ta > t1);
            }
        }
    }
    return sum * rd;
}

double WalleniusNCHypergeometric::integrate_step(double ta, double tb) {
    // integration subprocedure used by integrate()
    // makes one integration step from ta to tb using Gauss-Legendre method.
    // result is scaled by multiplication with exp(bico)
    double ab, delta, tau, ltau, y, sum, taur, rdm1;
    int i;

    // define constants for Gauss-Legendre integration with IPOINTS points
#define IPOINTS  8  // number of points in each integration step

#if   IPOINTS == 3
    const double xval[3] = { -.774596669241,0,0.774596668241 };
    const double weights[3] = { .5555555555555555,.88888888888888888,.55555555555555 };
#elif IPOINTS == 4
    const double xval[4] = { -0.861136311594,-0.339981043585,0.339981043585,0.861136311594 },
        const double weights[4] = { 0.347854845137,0.652145154863,0.652145154863,0.347854845137 };
#elif IPOINTS == 5
    const double xval[5] = { -0.906179845939,-0.538469310106,0,0.538469310106,0.906179845939 };
    const double weights[5] = { 0.236926885056,0.478628670499,0.568888888889,0.478628670499,0.236926885056 };
#elif IPOINTS == 6
    const double xval[6] = { -0.932469514203,-0.661209386466,-0.238619186083,0.238619186083,0.661209386466,0.932469514203 };
    const double weights[6] = { 0.171324492379,0.360761573048,0.467913934573,0.467913934573,0.360761573048,0.171324492379 };
#elif IPOINTS == 8
    const double xval[8] = { -0.960289856498,-0.796666477414,-0.525532409916,-0.183434642496,0.183434642496,0.525532409916,0.796666477414,0.960289856498 };
    const double weights[8] = { 0.10122853629,0.222381034453,0.313706645878,0.362683783378,0.362683783378,0.313706645878,0.222381034453,0.10122853629 };
#elif IPOINTS == 12
    const double xval[12] = { -0.981560634247,-0.90411725637,-0.769902674194,-0.587317954287,-0.367831498998,-0.125233408511,0.125233408511,0.367831498998,0.587317954287,0.769902674194,0.90411725637,0.981560634247 };
    const double weights[12] = { 0.0471753363866,0.106939325995,0.160078328543,0.203167426723,0.233492536538,0.249147045813,0.249147045813,0.233492536538,0.203167426723,0.160078328543,0.106939325995,0.0471753363866 };
#elif IPOINTS == 16
    const double xval[16] = { -0.989400934992,-0.944575023073,-0.865631202388,-0.755404408355,-0.617876244403,-0.458016777657,-0.281603550779,-0.0950125098376,0.0950125098376,0.281603550779,0.458016777657,0.617876244403,0.755404408355,0.865631202388,0.944575023073,0.989400934992 };
    const double weights[16] = { 0.027152459411,0.0622535239372,0.0951585116838,0.124628971256,0.149595988817,0.169156519395,0.182603415045,0.189450610455,0.189450610455,0.182603415045,0.169156519395,0.149595988817,0.124628971256,0.0951585116838,0.0622535239372,0.027152459411 };
#else
#error // IPOINTS must be a value for which the tables are defined
#endif

    delta = 0.5 * (tb - ta);
    ab = 0.5 * (ta + tb);
    rdm1 = rd - 1.;
    sum = 0;

    for (i = 0; i < IPOINTS; i++) {
        tau = ab + delta * xval[i];
        ltau = log(tau);
        taur = r * ltau;
        // possible loss of precision due to subtraction here:
        y = log1pow(taur * omega, x) + log1pow(taur, n - x) + rdm1 * ltau + bico;
        if (y > -50.) sum += weights[i] * exp(y);
    }
    return delta * sum;
}

double WalleniusNCHypergeometric::search_inflect(double t_from, double t_to) {
    // search for an inflection point of the integrand PHI(t) in the interval
    // t_from < t < t_to
    const int COLORS = 2;              // number of colors
    double t, t1;                      // independent variable
    double rho[COLORS];                // r*omega[i]
    double q;                          // t^rho[i] / (1-t^rho[i])
    double q1;                         // 1-t^rho[i]
    double xx[COLORS];                 // x[i]
    double zeta[COLORS][4][4];         // zeta[i,j,k] coefficients
    double phi[4];                     // derivatives of phi(t) = log PHI(t)
    double z2;                         // PHI''(t)/PHI(t)
    double zd;                         // derivative in Newton Raphson iteration
    double rdm1;                       // r * d - 1
    double tr;                         // 1/t
    double log2t;                      // log2(t)
    double method;                     // 0 for z2'(t) method, 1 for z3(t) method
    int i;                             // color
    int iter;                          // count iterations

    rdm1 = rd - 1.;
    if (t_from == 0 && rdm1 <= 1.) return 0.;    // no inflection point
    rho[0] = r * omega;  rho[1] = r;
    xx[0] = x;  xx[1] = n - x;
    t = 0.5 * (t_from + t_to);
    for (i = 0; i < COLORS; i++) {               // calculate zeta coefficients
        zeta[i][1][1] = rho[i];
        zeta[i][1][2] = rho[i] * (rho[i] - 1.);
        zeta[i][2][2] = rho[i] * rho[i];
        zeta[i][1][3] = zeta[i][1][2] * (rho[i] - 2.);
        zeta[i][2][3] = zeta[i][1][2] * rho[i] * 3.;
        zeta[i][3][3] = zeta[i][2][2] * rho[i] * 2.;
    }
    iter = 0;

    do {
        t1 = t;
        tr = 1. / t;
        log2t = log(t) * (1. / ln2);
        phi[1] = phi[2] = phi[3] = 0.;
        for (i = 0; i < COLORS; i++) {           // calculate first 3 derivatives of phi(t)
            q1 = pow2_1(rho[i] * log2t, &q);
            q /= q1;
            phi[1] -= xx[i] * zeta[i][1][1] * q;
            phi[2] -= xx[i] * q * (zeta[i][1][2] + q * zeta[i][2][2]);
            phi[3] -= xx[i] * q * (zeta[i][1][3] + q * (zeta[i][2][3] + q * zeta[i][3][3]));
        }
        phi[1] += rdm1;
        phi[2] -= rdm1;
        phi[3] += 2. * rdm1;
        phi[1] *= tr;
        phi[2] *= tr * tr;
        phi[3] *= tr * tr * tr;
        method = (iter & 2) >> 1;                // alternate between the two methods
        z2 = phi[1] * phi[1] + phi[2];
        zd = method * phi[1] * phi[1] * phi[1] + (2. + method) * phi[1] * phi[2] + phi[3];

        if (t < 0.5) {
            if (z2 > 0) {
                t_from = t;
            }
            else {
                t_to = t;
            }
            if (zd >= 0) {
                // use binary search if Newton-Raphson iteration makes problems
                t = (t_from ? 0.5 : 0.2) * (t_from + t_to);
            }
            else {
                // Newton-Raphson iteration
                t -= z2 / zd;
            }
        }
        else {
            if (z2 < 0) {
                t_from = t;
            }
            else {
                t_to = t;
            }
            if (zd <= 0) {
                // use binary search if Newton-Raphson iteration makes problems
                t = 0.5 * (t_from + t_to);
            }
            else {
                // Newton-Raphson iteration
                t -= z2 / zd;
            }
        }
        if (t >= t_to) t = (t1 + t_to) * 0.5;
        if (t <= t_from) t = (t1 + t_from) * 0.5;
        if (++iter > 20) reportError("Search for inflection point failed in function WalleniusNCHypergeometric::search_inflect");
    } while (fabs(t - t1) > 1E-5);
    return t;
}

int WalleniusNCHypergeometric::bernouilliH(int32_t x_, double h, double rh, RandomVariates *sto) {
    // This function generates a Bernouilli variate with probability proportional
    // to the univariate Wallenius' noncentral hypergeometric distribution.
    // The return value will be 1 with probability f(x_)/h and 0 with probability
    // 1-f(x_)/h.
    // This is equivalent to calling sto->bernouilli(probability(x_)/h),
    // but this method is faster. The method used here avoids calculating the
    // Wallenius probability by sampling in the t-domain.
    // rh is a uniform random number in the interval 0 <= rh < h. The function
    // uses additional random numbers generated from sto.
    // This function is intended for use in rejection methods for sampling from
    // the Wallenius distribution. It is called from walleniusNCHypRatioOfUnifoms
    double f0;                         // Lambda*Phi(0.5)
    double phideri0;                   // phi(0.5)/rd
    double qi;                         // 2^(-r*omega[i])
    double qi1;                        // 1-qi
    double omegai[2] = { omega,1. };     // weights for each color
    double romegi;                     // r*omega[i]
    double xi[2] = { double(x_), double(n - x_) };    // number of each color sampled
    double k;                          // adjusted width for majorizing function Ypsilon(t)
    double erfk;                       // erf correction
    double rdm1;                       // rd - 1
    double G_integral;                 // integral of majorizing function Ypsilon(t)
    double ts;                         // t sampled from Ypsilon(t) distribution
    double logts;                      // log(ts)
    double rlogts;                     // r*log(ts)
    double fts;                        // phi(ts)/rd
    double rgts;                       // 1/(Ypsilon(ts)/rd)
    double t2;                         // temporary in calculation of Ypsilon(ts)
    int i, j;                          // loop counters
    const double rsqrt8 = 0.3535533905932737622; // 1/sqrt(8)
    const double sqrt2pi = 2.506628274631000454; // sqrt(2*pi)

    x = x_;                            // save x in class object
    lnbico();                          // calculate bico = log(Lambda)
    findpars();                        // calculate r, d, rd, w, E
    if (E > 0.) {
        k = log(E);                    // correction for majorizing function
        k = 1. + 0.0271 * (k * sqrt(k));
    }
    else k = 1.;
    k *= w;                            // w * k   
    rdm1 = rd - 1.;

    // calculate phi(0.5)/rd
    phideri0 = -ln2 * rdm1;
    for (i = 0; i < 2; i++) {
        romegi = r * omegai[i];
        if (romegi > 40.) {
            qi = 0.;  qi1 = 1.;        // avoid underflow
        }
        else {
            qi1 = pow2_1(-romegi, &qi);
        }
        phideri0 += xi[i] * log1mx(qi, qi1);
    }

    erfk = erf(rsqrt8 / k);
    f0 = rd * exp(phideri0 + bico);
    G_integral = f0 * sqrt2pi * k * erfk;

    if (G_integral <= h) {             // G fits under h-hat
        do {
            ts = sto->normal(0, k);    // sample ts from normal distribution
        } while (fabs(ts) >= 0.5);     // reject values outside interval, and avoid ts = 0
        ts += 0.5;                     // ts = normal distributed in interval (0,1)

        for (fts = 0., j = 0; j < 2; j++) { // calculate (Phi(ts)+Phi(1-ts))/2
            logts = log(ts);  rlogts = r * logts; // (ts = 0 avoided above)
            fts += exp(log1pow(rlogts * omega, xi[0]) + log1pow(rlogts, xi[1]) + rdm1 * logts + bico);
            ts = 1. - ts;
        }
        fts *= 0.5;

        t2 = (ts - 0.5) / k;           // calculate 1/Ypsilon(ts)
        rgts = exp(-(phideri0 + bico - 0.5 * t2 * t2));
        return rh < G_integral * fts * rgts;   // Bernouilli variate
    }

    else { // G > h: can't use sampling in t-domain
        return rh < probability(x);
    }
}

double WalleniusNCHypergeometric::lnbico() {
    // natural log of binomial coefficients.
    // returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!)
    int32_t x2 = n - x, m2 = N - m;
    if (xLastBico < 0) { // m, n, N have changed
        mFac = lnFac(m) + lnFac(m2);
    }
    if (m < FAK_LEN && m2 < FAK_LEN)  goto DEFLT;
    switch (x - xLastBico) {
    case 0: // x unchanged
        break;
    case 1: // x incremented. calculate from previous value
        xFac += log(double(x) * (m2 - x2) / (double(x2 + 1) * (m - x + 1)));
        break;
    case -1: // x decremented. calculate from previous value
        xFac += log(double(x2) * (m - x) / (double(x + 1) * (m2 - x2 + 1)));
        break;
    default: DEFLT: // calculate all
        xFac = lnFac(x) + lnFac(x2) + lnFac(m - x) + lnFac(m2 - x2);
    }
    xLastBico = x;
    return bico = mFac - xFac;
}

void WalleniusNCHypergeometric::reportError(const char* text) {
    // report an error message via global object
    errors.reportError(text);
}


/***********************************************************************
 Various common functions:
***********************************************************************/

double pow2_1(double q, double * y0) {
    // calculate 2^q and (1-2^q) without loss of precision.
    // return value is (1-2^q). 2^q is returned in *y0
    double y, y1;
    q *= ln2;
    if (fabs(q) > 0.1) {
        y = exp(q);                    // 2^q
        y1 = 1. - y;                   // 1-2^q
    }
    else { // Use expm1
        y1 = expm1(q);                 // 2^q-1
        y = y1 + 1;                    // 2^q
        y1 = -y1;                      // 1-2^q
    }
    if (y0) *y0 = y;                   // Return y if not void pointer
    return y1;                         // Return y1
}

double log1mx(double x, double x1) {
    // Calculate log(1-x) without loss of precision when x is small.
    // Parameter x1 must be = 1-x.
    if (fabs(x) > 0.03) {
        return log(x1);
    }
    else { // use log1p(x) = log(1+x)
        return log1p(-x);
    }
}

double log1pow(double q, double x) {
    // calculate log((1-e^q)^x) without loss of precision.
    // Combines the methods of the above two functions.
    double y, y1;

    if (fabs(q) > 0.1) {
        y = exp(q);                    // e^q
        y1 = 1. - y;                   // 1-e^q
    }
    else { // Use expm1
        y1 = expm1(q);                 // e^q-1
        y = y1 + 1;                    // e^q
        y1 = -y1;                      // 1-e^q
    }

    if (y > 0.1) { // (1-y)^x calculated without problem
        return x * log(y1);
    }
    else { // Use log1p
        return x * log1p(-y);
    }
}

double fallingFactorial(double a, double b) {
    // calculates ln(a*(a-1)*(a-2)* ... * (a-b+1))

    if (b < 30 && int(b) == b && a < 1E10) {
        // direct calculation
        double f = 1.;
        for (int i = 0; i < b; i++) f *= a--;
        return log(f);
    }

    if (a > 100. * b && b > 1.) {
        // combine Stirling formulas for a and (a-b) to avoid loss of precision
        double ar = 1. / a;
        double cr = 1. / (a - b);
        // calculate -log(1-b/a) by Taylor expansion
        double s = 0., lasts, n = 1., ba = b * ar, f = ba;
        do {
            lasts = s;
            s += f / n;
            f *= ba;
            n++;
        } while (s != lasts);
        return (a + 0.5) * s + b * log(a - b) - b + (1. / 12.) * (ar - cr)
            /* - (1./360.)*(ar*ar*ar-cr*cr*cr) */;
    }
    // use lnFacr function
    return lnFacr(a) - lnFacr(a - b);
}

double lnFacr(double x) {
    // log factorial of non-integer x
    int32_t ix = (int32_t)(x);
    if (x == ix) return lnFac(ix);     // x is integer
    double r, r2, d = 1., f;
    static const double
        c0 = 0.918938533204672722,     // ln(sqrt(2*pi))
        c1 = 1. / 12.,
        c3 = -1. / 360.,
        c5 = 1. / 1260.,
        c7 = -1. / 1680.;
    if (x < 6.) {
        if (x == 0 || x == 1) return 0;
        while (x < 6) d *= ++x;
    }
    r = 1. / x;  r2 = r * r;
    f = (x + 0.5) * log(x) - x + c0 + r * (c1 + r2 * (c3 + r2 * (c5 + r2 * c7)));
    if (d != 1.) f -= log(d);
    return f;
}

int32_t floorLog2(float x) {
    // This function calculates floor(log2(x)) for positive x.
    // The return value is <= -127 for x <= 0.
    union UfloatInt {  // Union for extracting bits from a float
        float   f;
        int32_t i;
        UfloatInt(float ff) { f = ff; }    // constructor
    };
    // Running on a platform known to use IEEE-754 floating point format
    //int32_t n = *(int32_t*)&x;
    int32_t n = UfloatInt(x).i;
    return (n >> 23) - 0x7F;
}


/***********************************************************************
Multivariate Wallenius noncentral hypergeometric distribution
***********************************************************************/

void RandomVariates::multiWalleniusNCHyp(int32_t * destination,
    int32_t * source, double * weights, int32_t n, int colors) {
    /*
    This function generates a vector of random variables with the
    multivariate Wallenius noncentral hypergeometric distribution.

    The multivariate Wallenius noncentral hypergeometric distribution is
    the distribution you get when drawing colored balls from an urn
    with any number of colors, without replacement, and with bias.

    The weights are defined so that the probability of taking a particular
    ball is proportional to its weight.

    Parameters:
    destination:   An output array to receive the number of balls of each color.
                   Must have space for at least 'colors' elements.
    source:        An input array containing the number of balls of each color in the urn.
                   Must have 'colors' elements. All elements must be non-negative.
    weights:       The odds of each color. Must have 'colors' elements.
                   All elements must be non-negative.
    n:             The number of balls to draw from the urn.
                   Cannot exceed the total number of balls with nonzero weight in source.
    colors:        The number of possible colors.

    MAXCOLORS  (defined in random.h): You may adjust MAXCOLORS to the maximum
    number of colors you need.

    The function will reduce the number of colors, if possible, by eliminating
    colors with zero weight or zero number and pooling together colors with the
    same weight. The problem thus reduced is handled in the arrays osource,
    urn, oweights and osample of size colors2.

    The sampling proceeds by either of two methods: simulating urn experiment,
    or conditional method followed by Metropolis-Hastings sampling.

    Simulating the urn experiment is simply taking one ball at a time, requiring
    n uniform random variates. The problem is reduced whenever a color has been
    exhausted.

    The conditional method divides the colors into groups where the number of
    balls in each group is determined by sampling from the marginal distribution
    which is approximated by the univariate Wallenius distribution. Each group
    is then subdivided by sampling one color at a time until all colors have
    been sampled.

    The sample from the conditional method does not have the exact distribution,
    but it is used as a starting point for the Metropolis-Hastings sampling,
    which proceeds as follows: colors c1 and c2 are re-sampled using the
    univariate Wallenius distribution, keeping the samples of all other colors
    constant. The new sample is accepted or the old sample retained, according
    to the Metropolis formula which corrects for the slight error introduced
    by not using the true conditional distribution. c1 and c2 are rotated in
    an order determined by the variance of each color. This rotation (scan) is
    repeated nHastings times.
    */

    // variables 
    int order1[MAXCOLORS];             // sort order, index into source and destination
    int order2[MAXCOLORS];             // corresponding index into arrays when equal weights pooled together
    int order3[MAXCOLORS];             // secondary index for sorting by variance
    int32_t osource[MAXCOLORS];        // contents of source, sorted by weight with equal weights pooled together
    int32_t urn[MAXCOLORS];            // balls from osource not taken yet
    int32_t osample[MAXCOLORS];        // balls sampled
    double oweights[MAXCOLORS];        // sorted list of weights
    double wcum[MAXCOLORS];            // list of accumulated probabilities
    double var[MAXCOLORS];             // sorted list of variance
    double w = 0.;                     // weight of balls of one color
    double w1, w2;                     // odds within group; mean weight in group
    double wsum;                       // total weight of all balls of several or all colors
    double p;                          // probability
    double f0, f1;                     // multivariate probability function
    double g0, g1;                     // conditional probability function
    double r1, r2;                     // temporaries in calculation of variance
    int32_t nn;                        // number of balls left to sample
    int32_t m;                         // number of balls of one color
    int32_t msum;                      // total number of balls of several or all colors
    int32_t N;                         // total number of balls with nonzero weight
    int32_t x0, x = 0;                 // sample of one color
    int32_t n1, n2, ng;                // size of weight group sample or partial sample
    int32_t m1, m2;                    // size of weight group
    int i, j, k;                       // loop counters
    int c, c1, c2;                     // color index
    int colors2;                       // reduced number of colors
    int a, b;                          // color index delimiting weight group
    int nHastings;                     // number of scans in Metropolis-Hastings sampling

    // check validity of parameters
    if (n < 0 || colors < 0 || colors > MAXCOLORS) reportError("Parameter out of range in function MultiWalleniusNCHyp");
    if (colors == 0) return;
    if (n == 0) {
        for (i = 0; i < colors; i++) destination[i] = 0; return;
    }

    // check validity of array parameters
    for (i = 0, msum = 0; i < colors; i++) {
        m = source[i];  w = weights[i];
        if (m < 0 || w < 0) reportError("Parameter negative in function MultiWalleniusNCHyp");
        if (w) msum += m;
    }
    N = msum;

    // sort colors by weight, heaviest first
    for (i = 0; i < colors; i++) order1[i] = order3[i] = i;
    for (i = 0; i < colors - 1; i++) {
        c = order1[i];  k = i;
        w = weights[c];
        if (source[c] == 0) w = 0;     // zero number treated as zero weight
        for (j = i + 1; j < colors; j++) {
            c2 = order1[j];
            if (weights[c2] > w && source[c2]) {
                w = weights[c2];  k = j;
            }
        }
        order1[i] = order1[k];  order1[k] = c;
    }

    // skip any colors with zero weight or zero number.
    // this solves all problems with zero weights
    while (colors > 0 && (weights[c = order1[colors - 1]] == 0 || source[c] == 0)) {
        colors--;  destination[c] = 0;
    }

    // check if there are more than n balls with nonzero weight
    if (n >= N) {
        if (n > N) reportError("Taking more items than there are in function MultiWalleniusNCHyp");
        for (i = 0; i < colors; i++) {
            c = order1[i];  destination[c] = source[c];
        }
        return;
    }

    // copy source and weights into ordered lists 
    // and pool together colors with same weight
    for (i = 0, c2 = -1; i < colors; i++) {
        c = order1[i];
        if (i == 0 || weights[c] != w) {
            c2++;
            x = source[c];
            oweights[c2] = w = weights[c];
        }
        else {
            x += source[c];            // join colors with same weight
        }
        urn[c2] = osource[c2] = x;
        order2[i] = c2;
        osample[c2] = 0;
    }
    colors2 = c2 + 1;

    // check number of colors left
    if (colors2 < 3) {
        // simple cases
        if (colors2 == 1) osample[0] = n;
        if (colors2 == 2) {
            x = walleniusNCHyp(n, osource[0], N, oweights[0] / oweights[1]);
            osample[0] = x;  osample[1] = n - x;
        }
    }
    else {
        // more than 2 colors
        nn = n;

        // decide which method to use
        if (nn < 5000 * colors2) {

            // Simulate urn experiment

            // Make list of accumulated probabilities of each color
            for (i = 0, wsum = 0; i < colors2; i++) {
                wsum += urn[i] * oweights[i];
                wcum[i] = wsum;
            }

            // take one item nn times
            j = colors2 - 1;
            do {
                // get random color according to probability distribution wcum
                p = random() * wcum[colors2 - 1];
                // get color from search in probability distribution wcum
                for (i = 0; i < j; i++) {
                    if (p < wcum[i]) break;
                }

                // sample one ball of color i
                osample[i]++;  urn[i]--;  nn--;

                // check if this color has been exhausted
                if (urn[i] == 0) {
                    if (i != j) {
                        // put exhausted color at the end of lists so that colors2 can be reduced
                        m = osource[i]; osource[i] = osource[j]; osource[j] = m;
                        m = urn[i]; urn[i] = urn[j]; urn[j] = m;
                        m = osample[i]; osample[i] = osample[j]; osample[j] = m;
                        w = oweights[i]; oweights[i] = oweights[j]; oweights[j] = w;
                        // update order2 list (no longer sorted by weight)
                        for (k = 0; k < colors; k++) {
                            if (order2[k] == i) order2[k] = j; else
                                if (order2[k] == j) order2[k] = i;
                        }
                    }
                    colors2--;  j = colors2 - 1;  // decrement number of colors left in urn

                    if (colors2 == 2 && nn > 50) {
                        // two colors left. use univariate distribution for the rest
                        x = walleniusNCHyp(nn, urn[0], urn[0] + urn[1], oweights[0] / oweights[1]);
                        osample[0] += x;
                        osample[1] += nn - x;
                        break;
                    }

                    if (colors2 == 1) {
                        // only one color left. The rest is deterministic
                        osample[0] += nn;
                        break;
                    }

                    // make sure wcum is re-calculated from beginning
                    i = 0;
                }

                // update list of accumulated probabilities
                wsum = i > 0 ? wcum[i - 1] : 0.;
                for (k = i; k < colors2; k++) {
                    wsum += urn[k] * oweights[k];
                    wcum[k] = wsum;
                }
            } while (nn);
        }

        else {
            // use conditional method to make starting point for
            // Metropolis-Hastings sampling

            // divide weights into two groups, heavy and light
            a = 0;  b = colors2 - 1;
            w = sqrt(oweights[0] * oweights[colors2 - 1]);
            do {
                c = (a + b) / 2;
                if (oweights[c] > w) a = c; else b = c;
            } while (b > a + 1);
            // heavy group goes from 0 to b-1, light group goes from b to colors2-1

            // calculate mean weight for heavy color group
            for (i = 0, m1 = 0, wsum = 0; i < b; i++) {
                m1 += urn[i];  wsum += oweights[i] * urn[i];
            }
            w1 = wsum / m1;

            // calculate mean weight for light color group
            for (i = b, m2 = 0, wsum = 0; i < colors2; i++) {
                m2 += urn[i];  wsum += oweights[i] * urn[i];
            }
            w2 = wsum / m2;

            // split partial sample n into heavy (n1) and light (n2)
            n1 = walleniusNCHyp(n, m1, m1 + m2, w1 / w2);
            n2 = n - n1;

            // set parameters for first group (heavy)
            a = 0;  ng = n1;

            // loop twice, for the two groops
            for (k = 0; k < 2; k++) {

                // split group into single colors by calling univariate distribution b-a-1 times
                for (i = a; i < b - 1; i++) {
                    m = urn[i];  w = oweights[i];

                    // calculate mean weight of remaining colors
                    for (j = i + 1, msum = 0, wsum = 0; j < b; j++) {
                        m1 = urn[j];  w1 = oweights[j];
                        msum += m1;  wsum += m1 * w1;
                    }

                    // sample color i in group
                    x = wsum ? walleniusNCHyp(ng, m, msum + m, w * msum / wsum) : ng;

                    osample[i] = x;
                    ng -= x;
                }

                // get the last one in the group
                osample[i] = ng;

                // set parameters for second group (light)
                a = b;  b = colors2;  ng = n2;
            }

            // finished with conditional method. 
            // osample contains starting point for Metropolis-Hastings sampling

            // make object for calculating probabilities and mean
            MultiWalleniusNCHypergeometric wmnc(n, osource, oweights, colors2);

            wmnc.mean(var); // calculate mean
            // calculate approximate variance from mean
            for (i = 0; i < colors; i++) {
                r1 = var[i] * (osource[i] - var[i]);
                r2 = (n - var[i]) * (var[i] + N - n - osource[i]);
                if (r1 <= 0. || r2 <= 0.) {
                    var[i] = 0.;
                }
                else {
                    var[i] = N * r1 * r2 / ((N - 1) * (osource[i] * r2 + (N - osource[i]) * r1));
                }
            }

            // sort again, this time by variance
            for (i = 0; i < colors2 - 1; i++) {
                c = order3[i];  k = i;
                w = var[c];
                for (j = i + 1; j < colors2; j++) {
                    c2 = order3[j];
                    if (var[c2] > w) {
                        w = var[c2];  k = j;
                    }
                }
                order3[i] = order3[k];  order3[k] = c;
            }

            // number of scans (this value of nHastings has not been fine-tuned)
            nHastings = 4;
            if (accuracy < 1E-6) nHastings = 6;
            if (colors2 > 5) nHastings++;

            // Metropolis-Hastings sampler
            f0 = -1.;
            for (k = 0; k < nHastings; k++) {
                for (i = 0; i < colors2; i++) {
                    j = i + 1;
                    if (j >= colors2) j = 0;
                    c1 = order3[i];  c2 = order3[j];
                    w = oweights[c1] / oweights[c2];
                    n1 = osample[c1] + osample[c2];
                    x0 = osample[c1];
                    x = walleniusNCHyp(n1, osource[c1], osource[c1] + osource[c2], w);
                    if (x == x0) continue; // accepted
                    if (f0 < 0.) f0 = wmnc.probability(osample);
                    WalleniusNCHypergeometric nc(n1, osource[c1], osource[c1] + osource[c2], w, accuracy);
                    g0 = nc.probability(x0);
                    g1 = nc.probability(x);
                    osample[c1] = x;
                    osample[c2] = n1 - x;
                    f1 = wmnc.probability(osample);
                    g0 = f1 * g0;  g1 = f0 * g1;
                    if (g0 >= g1 || g0 > g1 * random()) {
                        // new state accepted
                        f0 = -1.;
                    }
                    else {
                        // rejected. restore old sample
                        osample[c1] = x0;
                        osample[c2] = n1 - x0;
                    }
                }
            }
        }
    }

    // finished sampling by either method
    // un-sort sample into destination and untangle re-orderings
    for (i = 0; i < colors; i++) {
        c1 = order1[i];  c2 = order2[i];
        if (source[c1] == osource[c2]) {
            destination[c1] = osample[c2];
        }
        else {
            // split colors with same weight that have been treated as one
            x = hypergeometric(osample[c2], source[c1], osource[c2]);
            destination[c1] = x;
            osample[c2] -= x;
            osource[c2] -= source[c1];
        }
    }
}

void RandomVariates::complMultiWalleniusNCHyp(int32_t * destination, int32_t * source, double * weights, int32_t n, int colors) {
    // Complementary Multivariate Wallenius noncentral hypergeometric distribution
    // This function models negative selection with more than two phenotypes or genotypes.
    // The number of individuals killed follows the Wallenius distribution,
    // The number of individuals surviving follows the complementary Wallenius distribution.
    int32_t killed[MAXCOLORS];                   // number killed
    int32_t source1[MAXCOLORS];                  // same as source, excluding colors with zero weight
    double reciprocal_weights[MAXCOLORS];        // 1./weights
    int32_t N = 0;                               // total individuals = sum of source[i]
    int i;                                       // loop counter
    for (i = 0; i < colors; i++) {
        source1[i] = source[i];
        if (weights[i] != 0.) {                  // avoid dividing by zero
            reciprocal_weights[i] = 1. / weights[i];
        }
        else {
            reciprocal_weights[i] = 0;
            source1[i] = 0;
        }
        N += source1[i];                         // sum of individuals with nonzero weight
    }
    multiWalleniusNCHyp(killed, source1, reciprocal_weights, N - n, colors); // find number killed
    for (i = 0; i < colors; i++) {
        destination[i] = source1[i] - killed[i]; // find number surviving
    }
}


// calculation methods in class MultiWalleniusNCHypergeometric:

MultiWalleniusNCHypergeometric::MultiWalleniusNCHypergeometric(int32_t n_, int32_t * m_, double * odds_, int colors_, double accuracy_) {
    // constructor
    accuracy = accuracy_;
    setParameters(n_, m_, odds_, colors_);
}

void MultiWalleniusNCHypergeometric::setParameters(int32_t n_, int32_t * m_, double * odds_, int colors_) {
    // change parameters
    int32_t N1;
    int i;
    n = n_;  m = m_;  omega = odds_;  colors = colors_;
    r = 1.;
    for (N = N1 = 0, i = 0; i < colors; i++) {
        if (m[i] < 0 || omega[i] < 0) reportError("Parameter negative in constructor for CMultiWalleniusNCHypergeometric");
        N += m[i];
        if (omega[i]) N1 += m[i];
    }
    if (N < n) reportError("Not enough items in constructor for MultiWalleniusNCHypergeometric");
    if (N1 < n) reportError("Not enough items with nonzero weight in constructor for MultiWalleniusNCHypergeometric");
}

void MultiWalleniusNCHypergeometric::mean(double * mu) {
    // calculate approximate mean of multivariate Wallenius noncentral hypergeometric 
    // distribution. Result is returned in mu[0..colors-1]
    double omeg[MAXCOLORS];            // scaled weights
    double omr;                        // reciprocal mean weight
    double t, t1;                      // independent variable in iteration
    double to, to1;                    // exp(t*omega[i]), 1-exp(t*omega[i])
    double h;                          // function to find root of
    double hd;                         // derivative of h
    double dummy;                      // unused return
    int i;                             // color index
    int iter;                          // number of iterations

    if (n == 0) {
        // needs special case
        for (i = 0; i < colors; i++) {
            mu[i] = 0.;
        }
        return;
    }

    // calculate mean weight
    for (omr = 0., i = 0; i < colors; i++) omr += omega[i] * m[i];
    omr = N / omr;
    // scale weights to make mean = 1
    for (i = 0; i < colors; i++) omeg[i] = omega[i] * omr;
    // Newton Raphson iteration
    iter = 0;  t = -1.;                // first guess
    do {
        t1 = t;
        h = hd = 0.;
        // calculate H and HD
        for (i = 0; i < colors; i++) {
            if (omeg[i] != 0.) {
                to1 = pow2_1(t * (1. / ln2) * omeg[i], &to);
                h += m[i] * to1;
                hd -= m[i] * omeg[i] * to;
            }
        }
        t -= (h - n) / hd;
        if (t >= 0) t = 0.5 * t1;
        if (++iter > 20) {
            reportError("Search for mean failed in function MultiWalleniusNCHypergeometric::mean");
        }
    } while (fabs(h - n) > 1E-3);
    // finished iteration. Get all mu[i]
    for (i = 0; i < colors; i++) {
        if (omeg[i] != 0.) {
            to1 = pow2_1(t * (1. / ln2) * omeg[i], &dummy);
            mu[i] = m[i] * to1;
        }
        else {
            mu[i] = 0.;
        }
    }
}

double MultiWalleniusNCHypergeometric::probability(int32_t * x_) {
    // calculate probability function. choosing best method
    int i, j, em;
    bool central;
    int32_t xsum;
    x = x_;

    for (xsum = i = 0; i < colors; i++)  xsum += x[i];
    if (xsum != n) {
        reportError("sum of x values not equal to n in function MultiWalleniusNCHypergeometric::probability");
    }

    if (colors < 3) {
        if (colors <= 0) return 1.;
        if (colors == 1) return x[0] == m[0];
        // colors = 2
        if (omega[1] == 0.) return x[0] == m[0];
        return WalleniusNCHypergeometric(n, m[0], N, omega[0] / omega[1], accuracy).probability(x[0]);
    }

    // check if weights are equal
    central = true;
    for (i = j = em = 0; i < colors; i++) {
        if (x[i] > m[i] || x[i] < 0 || x[i] < n - N + m[i]) return 0.;
        if (x[i] > 0) j++;
        if (omega[i] == 0. && x[i]) return 0.;
        if (x[i] == m[i] || omega[i] == 0.) em++;
        if (i > 0 && omega[i] != omega[i - 1]) central = false;
    }

    if (n == 0 || em == colors) return 1.;

    if (central) {
        // All omega's are equal. 
        // This is multivariate central hypergeometric distribution
        int32_t sx = n, sm = N;
        double p = 1.;
        for (i = 0; i < colors - 1; i++) {
            // Use univariate hypergeometric (usedcolors-1) times
            p *= WalleniusNCHypergeometric(sx, m[i], sm, 1.).probability(x[i]);
            sx -= x[i];  sm -= m[i];
        }
        return p;
    }

    if (j == 1) {
        return binoexpand();
    }

    findpars();
    if (w < 0.04 && E < 10 && (!em || w > 0.004)) {
        return laplace();
    }

    return integrate();
}

void MultiWalleniusNCHypergeometric::findpars(void) {
    // calculate d, E, r, w

    // find r to center peak of integrand at 0.5
    double dd;                         // scaled d
    double dr;                         // 1/d

    double z, zd, rr, lastr, rrc, rt, r2, r21, a, b, ro, k1, dummy;
    double omax;                       // highest omega
    double omaxr;                      // 1/omax
    double omeg[MAXCOLORS];            // scaled weights
    int i, j = 0;

    // find highest omega
    for (omax = 0., i = 0; i < colors; i++) {
        if (omega[i] > omax) omax = omega[i];
    }
    omaxr = 1. / omax;
    dd = E = 0.;
    for (i = 0; i < colors; i++) {
        // scale weights to make max = 1
        omeg[i] = omega[i] * omaxr;
        // calculate d and E
        dd += omeg[i] * (m[i] - x[i]);
        E += omeg[i] * m[i];
    }
    dr = 1. / dd;
    E *= dr;
    rr = r * omax;
    if (rr <= dr) rr = 1.2 * dr;       // initial guess
    // Newton-Raphson iteration to find r
    do {
        lastr = rr;
        rrc = 1. / rr;
        z = dd - rrc;                  // z(r)
        zd = rrc * rrc;                // z'(r)
        for (i = 0; i < colors; i++) {
            rt = rr * omeg[i];
            if (rt < 100. && rt > 0) { // avoid overflow and division by 0
                r21 = pow2_1(rt, &r2); // r2=2^r, r21=1.-2^r
                a = omeg[i] / r21;     // omegai/(1.-2^r)
                b = x[i] * a;          // x*omegai/(1.-2^r)
                z += b;
                zd += b * a * r2 * ln2;
            }
        }
        if (zd == 0) reportError("can't find r in function MultiWalleniusNCHypergeometric::findpars");
        rr -= z / zd;                  // next r
        if (rr <= dr) rr = lastr * 0.125 + dr * 0.875;
        if (++j == 70) reportError("convergence problem searching for r in function MultiWalleniusNCHypergeometric::findpars");
    } while (fabs(rr - lastr) > rr * 1.E-5);
    rd = rr * dd;
    r = rr * omaxr;

    // find peak width
    phi2d = 0.;
    for (i = 0; i < colors; i++) {
        ro = rr * omeg[i];
        if (ro < 300 && ro > 0.) {     // avoid overflow and division by 0
            k1 = pow2_1(ro, &dummy);
            k1 = -1. / k1;
            k1 = omeg[i] * omeg[i] * (k1 + k1 * k1);
        }
        else k1 = 0.;
        phi2d += x[i] * k1;
    }
    phi2d *= -4. * rr * rr;
    if (phi2d > 0.) reportError("peak width undefined in function MultiWalleniusNCHypergeometric::findpars");
    wr = sqrt(-phi2d);  w = 1. / wr;
}

double MultiWalleniusNCHypergeometric::lnbico(void) {
    // natural log of binomial coefficients
    bico = 0.;
    int i;
    for (i = 0; i < colors; i++) {
        if (x[i] < m[i] && omega[i]) {
            bico += lnFac(m[i]) - lnFac(x[i]) - lnFac(m[i] - x[i]);
        }
    }
    return bico;
}

// implementations of different calculation methods:

double MultiWalleniusNCHypergeometric::binoexpand(void) {
    // binomial expansion of integrand
    // only implemented for x[i] = 0 for all but one i
    int i, j, k;
    double W = 0.;                     // total weight
    for (i = j = k = 0; i < colors; i++) {
        W += omega[i] * m[i];
        if (x[i]) {
            j = i; k++;                // find the nonzero x[i]
        }
    }
    if (k > 1) reportError("More than one x[i] nonzero in MultiWalleniusNCHypergeometric::binoexpand");
    return exp(fallingFactorial(m[j], n) - fallingFactorial(W / omega[j], n));
}

double MultiWalleniusNCHypergeometric::laplace(void) {
    // Laplace's method with narrow integration interval, 
    // using error function residues table, defined below
    // Note that this function can only be used when the integrand peak is narrow.
    // findpars() must be called before this function.

    const int MAXDEG = 40;             // arraysize
    int degree;                        // max expansion degree
    double accur;                      // stop expansion when terms below this threshold
    double f0;                         // factor outside integral
    double rho[MAXCOLORS];             // r*omegai
    double qi;                         // 2^(-rho)
    double qi1;                        // 1-qi
    double qq[MAXCOLORS];              // qi / qi1
    double eta[MAXCOLORS + 1][MAXDEG + 1]; // eta coefficients
    double phideri[MAXDEG + 1];        // derivatives of phi
    double PSIderi[MAXDEG + 1];        // derivatives of PSI
    double * erfresp;                  // pointer to table of error function residues

    // variables in asymptotic summation
    static const double sqrt8 = 2.828427124746190098; // sqrt(8)
    double qqpow;                      // qq^j
    double pow2k;                      // 2^k
    double bino;                       // binomial coefficient  
    double vr;                         // 1/v, v = integration interval
    double v2m2;                       // (2*v)^(-2)
    double v2mk1;                      // (2*v)^(-k-1)
    double s;                          // summation term
    double sum;                        // Taylor sum

    int i;                             // loop counter for color
    int j;                             // loop counter for derivative
    int k;                             // loop counter for expansion degree
    int ll;                            // k/2
    int converg = 0;                   // number of consequtive terms below accuracy
    int precisionIndex;                // index into ErfRes table according to desired precision

    // initialize
    for (k = 0; k <= 2; k++)  phideri[k] = PSIderi[k] = 0;

    // find rho[i], qq[i], first eta coefficients, and zero'th derivative of phi
    for (i = 0; i < colors; i++) {
        rho[i] = r * omega[i];
        if (rho[i] == 0.) continue;
        if (rho[i] > 40.) {
            qi = 0.;  qi1 = 1.;        // avoid underflow
        }
        else {
            qi1 = pow2_1(-rho[i], &qi);// qi=2^(-rho), qi1=1.-2^(-rho)
        }
        qq[i] = qi / qi1;              // 2^(-r*omegai)/(1.-2^(-r*omegai))
        // peak = zero'th derivative
        phideri[0] += x[i] * log1mx(qi, qi1);
        // eta coefficients
        eta[i][0] = 0.;
        eta[i][1] = eta[i][2] = rho[i] * rho[i];
    }

    // d, r, and w must be calculated by findpars()
    // zero'th derivative
    phideri[0] -= (rd - 1.) * ln2;
    // scaled factor outside integral
    f0 = rd * exp(phideri[0] + lnbico());
    // calculate narrowed integration interval
    vr = sqrt8 * w;
    phideri[2] = phi2d;

    // get table according to desired precision
    precisionIndex = (-floorLog2((float)accuracy) - ERFRES_B + ERFRES_S - 1) / ERFRES_S;
    if (precisionIndex < 0) precisionIndex = 0;
    if (precisionIndex > ERFRES_N - 1) precisionIndex = ERFRES_N - 1;
    while (w * numSDev[precisionIndex] > 0.3) {
        // check if integration interval is too wide
        if (precisionIndex == 0) {
            reportError("Laplace method failed. Peak width too high in function MultiWalleniusNCHypergeometric::laplace");
            break;
        }
        precisionIndex--;              // reduce precision to keep integration interval narrow
    }
    erfresp = erfRes[precisionIndex];  // choose desired table

    degree = MAXDEG;                   // max expansion degree
    if (degree >= ERFRES_L * 2) degree = ERFRES_L * 2 - 2;

    // set up for starting loop at k=3
    v2m2 = 0.25 * vr * vr;             // (2*v)^(-2)
    PSIderi[0] = 1.;
    pow2k = 8.;
    sum = 0.5 * vr * erfresp[0];
    v2mk1 = 0.5 * vr * v2m2 * v2m2;
    accur = accuracy * sum;

    // summation loop
    for (k = 3; k <= degree; k++) {
        phideri[k] = 0.;

        // loop for all colors
        for (i = 0; i < colors; i++) {
            if (rho[i] == 0.) continue;
            eta[i][k] = 0.;
            // backward loop for all powers
            for (j = k; j > 0; j--) {
                // find coefficients recursively from previous coefficients
                eta[i][j] = eta[i][j] * (j * rho[i] - (k - 2)) + eta[i][j - 1] * rho[i] * (j - 1);
            }
            qqpow = 1.;
            // forward loop for all powers
            for (j = 1; j <= k; j++) {
                qqpow *= qq[i];        // qq^j
                // contribution to derivative
                phideri[k] += x[i] * eta[i][j] * qqpow;
            }
        }

        // finish calculation of derivatives
        phideri[k] = -pow2k * phideri[k] + 2 * (1 - k) * phideri[k - 1];

        pow2k *= 2.;                   // 2^k

        // loop to calculate derivatives of PSI from derivatives of psi.
        // terms # 0, 1, 2, k-2, and k-1 are zero and not included in loop.
        // The j'th derivatives of psi are identical to the derivatives of phi for j>2, and
        // zero for j=1,2. Hence we are using phideri[j] for j>2 here.
        PSIderi[k] = phideri[k];       // this is term # k
        bino = 0.5 * (k - 1) * (k - 2);// binomial coefficient for term # 3
        for (j = 3; j < k - 2; j++) {  // loop for remaining nonzero terms (if k>5)
            PSIderi[k] += PSIderi[k - j] * phideri[j] * bino;
            bino *= double(k - j) / double(j);
        }

        if ((k & 1) == 0) { // only for even k
            ll = k / 2;
            s = PSIderi[k] * v2mk1 * erfresp[ll];
            sum += s;

            // check for convergence of Taylor expansion
            if (fabs(s) < accur) converg++; else converg = 0;
            if (converg > 1) break;

            // update recursive expressions
            v2mk1 *= v2m2;
        }
    }

    // multiply by terms outside integral  
    return f0 * sum;
}

double MultiWalleniusNCHypergeometric::integrate(void) {
    // Wallenius non-central hypergeometric distribution function
    // calculation by numerical integration with variable-length steps
    // NOTE: findpars() must be called before this function.
    double s;                          // result of integration step
    double sum;                        // integral
    double ta, tb;                     // subinterval for integration step

    lnbico();                          // compute log of binomial coefficients

    // choose method:
    if (w < 0.02) {
        // normal method. Step length determined by peak width w
        double delta, s1;
        s1 = accuracy < 1E-9 ? 0.5 : 1.;
        delta = s1 * w;                          // integration steplength
        ta = 0.5 + 0.5 * delta;
        sum = integrate_step(1. - ta, ta);       // first integration step around center peak
        do {
            tb = ta + delta;
            if (tb > 1.) tb = 1.;
            s = integrate_step(ta, tb);          // integration step to the right of peak
            s += integrate_step(1. - tb, 1. - ta);   // integration step to the left of peak
            sum += s;
            if (s < accuracy * sum) break;       // stop before interval finished if accuracy reached
            ta = tb;
            if (tb > 0.5 + w) delta *= 2.;       // increase step length far from peak
        } while (tb < 1.);
    }

    else {
        // difficult situation. Step length determined by inflection points
        double t1, t2, tinf, delta, delta1;
        sum = 0.;
        // do left and right half of integration interval separately:
        for (t1 = 0., t2 = 0.5; t1 < 1.; t1 += 0.5, t2 += 0.5) {
            // integrate from 0 to 0.5 or from 0.5 to 1
            tinf = search_inflect(t1, t2);                 // find inflection point
            delta = tinf - t1; if (delta > t2 - tinf) delta = t2 - tinf; // distance to nearest endpoint
            delta *= 1. / 7.;                              // 1/7 will give 3 steps to nearest endpoint
            if (delta < 1E-4) delta = 1E-4;
            delta1 = delta;
            // integrate from tinf forwards to t2
            ta = tinf;
            do {
                tb = ta + delta1;
                if (tb > t2 - 0.25 * delta1) tb = t2;      // last step of this subinterval
                s = integrate_step(ta, tb);                // integration step
                sum += s;
                delta1 *= 2;                               // double steplength
                if (s < sum * 1E-4) delta1 *= 8.;          // large step when s small
                ta = tb;
            } while (tb < t2);
            if (tinf) {
                // integrate from tinf backwards to t1
                tb = tinf;
                do {
                    ta = tb - delta;
                    if (ta < t1 + 0.25 * delta) ta = t1;   // last step of this subinterval
                    s = integrate_step(ta, tb);            // integration step
                    sum += s;
                    delta *= 2;                            // double steplength
                    if (s < sum * 1E-4) delta *= 8.;       // large step when s small
                    tb = ta;
                } while (ta > t1);
            }
        }
    }
    return sum * rd;
}

double MultiWalleniusNCHypergeometric::integrate_step(double ta, double tb) {
    // integration subprocedure used by MultiWalleniusNCHypergeometric::integrate()
    // makes one integration step from ta to tb using Gauss-Legendre method.
    // result is scaled by multiplication with exp(bico)
    double ab, delta, tau, ltau, y, sum, taur, rdm1;
    int i, j;

    // define constants for Gauss-Legendre integration with IPOINTS points
#define IPOINTS  8  // number of points in each integration step

#if   IPOINTS == 3
    static const double xval[3] = { -.774596669241,0,0.774596668241 };
    static const double weights[3] = { .5555555555555555,.88888888888888888,.55555555555555 };
#elif IPOINTS == 4
    static const double xval[4] = { -0.861136311594,-0.339981043585,0.339981043585,0.861136311594 },
        static const double weights[4] = { 0.347854845137,0.652145154863,0.652145154863,0.347854845137 };
#elif IPOINTS == 5
    static const double xval[5] = { -0.906179845939,-0.538469310106,0,0.538469310106,0.906179845939 };
    static const double weights[5] = { 0.236926885056,0.478628670499,0.568888888889,0.478628670499,0.236926885056 };
#elif IPOINTS == 6
    static const double xval[6] = { -0.932469514203,-0.661209386466,-0.238619186083,0.238619186083,0.661209386466,0.932469514203 };
    static const double weights[6] = { 0.171324492379,0.360761573048,0.467913934573,0.467913934573,0.360761573048,0.171324492379 };
#elif IPOINTS == 8
    static const double xval[8] = { -0.960289856498,-0.796666477414,-0.525532409916,-0.183434642496,0.183434642496,0.525532409916,0.796666477414,0.960289856498 };
    static const double weights[8] = { 0.10122853629,0.222381034453,0.313706645878,0.362683783378,0.362683783378,0.313706645878,0.222381034453,0.10122853629 };
#elif IPOINTS == 12
    static const double xval[12] = { -0.981560634247,-0.90411725637,-0.769902674194,-0.587317954287,-0.367831498998,-0.125233408511,0.125233408511,0.367831498998,0.587317954287,0.769902674194,0.90411725637,0.981560634247 };
    static const double weights[12] = { 0.0471753363866,0.106939325995,0.160078328543,0.203167426723,0.233492536538,0.249147045813,0.249147045813,0.233492536538,0.203167426723,0.160078328543,0.106939325995,0.0471753363866 };
#elif IPOINTS == 16
    static const double xval[16] = { -0.989400934992,-0.944575023073,-0.865631202388,-0.755404408355,-0.617876244403,-0.458016777657,-0.281603550779,-0.0950125098376,0.0950125098376,0.281603550779,0.458016777657,0.617876244403,0.755404408355,0.865631202388,0.944575023073,0.989400934992 };
    static const double weights[16] = { 0.027152459411,0.0622535239372,0.0951585116838,0.124628971256,0.149595988817,0.169156519395,0.182603415045,0.189450610455,0.189450610455,0.182603415045,0.169156519395,0.149595988817,0.124628971256,0.0951585116838,0.0622535239372,0.027152459411 };
#else
#error // IPOINTS must be a value for which the tables are defined
#endif

    delta = 0.5 * (tb - ta);
    ab = 0.5 * (ta + tb);
    rdm1 = rd - 1.;
    sum = 0;

    for (j = 0; j < IPOINTS; j++) {
        tau = ab + delta * xval[j];
        ltau = log(tau);
        taur = r * ltau;
        y = 0.;
        for (i = 0; i < colors; i++) {
            // possible loss of precision due to subtraction here:
            if (omega[i]) {
                y += log1pow(taur * omega[i], x[i]);       // ln((1-e^taur*omegai)^xi)
            }
        }
        y += rdm1 * ltau + bico;
        if (y > -50.) sum += weights[j] * exp(y);
    }
    return delta * sum;
}

double MultiWalleniusNCHypergeometric::search_inflect(double t_from, double t_to) {
    // search for an inflection point of the integrand PHI(t) in the interval
    // t_from < t < t_to
    double t, t1;                      // independent variable
    double rho[MAXCOLORS];             // r*omega[i]
    double q;                          // t^rho[i] / (1-t^rho[i])
    double q1;                         // 1-t^rho[i]
    double zeta[MAXCOLORS][4][4];      // zeta[i,j,k] coefficients
    double phi[4];                     // derivatives of phi(t) = log PHI(t)
    double z2;                         // PHI''(t)/PHI(t)
    double zd;                         // derivative in Newton Raphson iteration
    double rdm1;                       // r * d - 1
    double tr;                         // 1/t
    double log2t;                      // log2(t)
    double method;                     // 0 for z2'(t) method, 1 for z3(t) method
    int i;                             // color
    int iter;                          // count iterations

    rdm1 = rd - 1.;
    if (t_from == 0 && rdm1 <= 1.) return 0.;    //no inflection point
    t = 0.5 * (t_from + t_to);
    for (i = 0; i < colors; i++) {               // calculate zeta coefficients
        rho[i] = r * omega[i];
        zeta[i][1][1] = rho[i];
        zeta[i][1][2] = rho[i] * (rho[i] - 1.);
        zeta[i][2][2] = rho[i] * rho[i];
        zeta[i][1][3] = zeta[i][1][2] * (rho[i] - 2.);
        zeta[i][2][3] = zeta[i][1][2] * rho[i] * 3.;
        zeta[i][3][3] = zeta[i][2][2] * rho[i] * 2.;
    }
    iter = 0;

    do {
        t1 = t;
        tr = 1. / t;
        log2t = log(t) * (1. / ln2);
        phi[1] = phi[2] = phi[3] = 0.;
        for (i = 0; i < colors; i++) {           // calculate first 3 derivatives of phi(t)
            if (rho[i] == 0.) continue;
            q1 = pow2_1(rho[i] * log2t, &q);
            q /= q1;
            phi[1] -= x[i] * zeta[i][1][1] * q;
            phi[2] -= x[i] * q * (zeta[i][1][2] + q * zeta[i][2][2]);
            phi[3] -= x[i] * q * (zeta[i][1][3] + q * (zeta[i][2][3] + q * zeta[i][3][3]));
        }
        phi[1] += rdm1;
        phi[2] -= rdm1;
        phi[3] += 2. * rdm1;
        phi[1] *= tr;
        phi[2] *= tr * tr;
        phi[3] *= tr * tr * tr;
        method = (iter & 2) >> 1;                // alternate between the two methods
        z2 = phi[1] * phi[1] + phi[2];
        zd = method * phi[1] * phi[1] * phi[1] + (2. + method) * phi[1] * phi[2] + phi[3];

        if (t < 0.5) {
            if (z2 > 0) {
                t_from = t;
            }
            else {
                t_to = t;
            }
            if (zd >= 0) {
                // use binary search if Newton-Raphson iteration makes problems
                t = (t_from ? 0.5 : 0.2) * (t_from + t_to);
            }
            else {
                // Newton-Raphson iteration
                t -= z2 / zd;
            }
        }
        else {
            if (z2 < 0) {
                t_from = t;
            }
            else {
                t_to = t;
            }
            if (zd <= 0) {
                // use binary search if Newton-Raphson iteration makes problems
                t = 0.5 * (t_from + t_to);
            }
            else {
                // Newton-Raphson iteration
                t -= z2 / zd;
            }
        }
        if (t >= t_to) t = (t1 + t_to) * 0.5;
        if (t <= t_from) t = (t1 + t_from) * 0.5;
        if (++iter > 20) reportError("Search for inflection point failed in function MultiWalleniusNCHypergeometric::search_inflect");
    } while (fabs(t - t1) > 1E-5);
    return t;
}

// Members of class MultiWalleniusNCHypergeometricMoments:

// calculation of mean and variance by complete enumeration of all x-combinations
double MultiWalleniusNCHypergeometricMoments::moments(double * mu, double * variance, int32_t * combinations) {
    // calculates mean and variance of multivariate Wallenius noncentral 
    // hypergeometric distribution by calculating all combinations of x-values.
    // Return value = sum of all probabilities. The deviation of this value 
    // from 1 is a measure of the accuracy.
    // Returns the mean to mean[0...colors-1]
    // Returns the variance to variance[0...colors-1]
    double sumf;                       // sum of all f(x) values
    int32_t msum;                      // temporary sum
    int i;                             // loop counter

    // get approximate mean
    mean(sx);
    // round mean to integers
    for (i = 0; i < colors; i++) {
        xm[i] = (int32_t)(sx[i] + 0.4999999);
    }

    // set up for recursive loops
    for (i = colors - 1, msum = 0; i >= 0; i--) {
        remaining[i] = msum;  msum += m[i];
    }
    for (i = 0; i < colors; i++) sx[i] = sxx[i] = 0.;
    sn = 0;

    // recursive loops to calculate sums  
    sumf = loop(n, 0);

    // calculate mean and variance
    for (i = 0; i < colors; i++) {
        mu[i] = sx[i] / sumf;
        variance[i] = sxx[i] / sumf - sx[i] * sx[i] / (sumf * sumf);
        //variance[i] = sxx[i] - mu[i] * mu[i];
        if (variance[i] < 0.) variance[i] = 0.;
    }

    // return combinations and sum
    if (combinations) *combinations = sn;
    return sumf;
}

double MultiWalleniusNCHypergeometricMoments::loop(int32_t n, int c) {
    // recursive function to loop through all combinations of x-values.
    // used by moments()
    int32_t x, x0;                     // x of color c
    int32_t xmin, xmax;                // min and max of x[c]
    double s1, s2, sum = 0.;           // sum of f(x) values
    int i;                             // loop counter

    if (c < colors - 1) {
        // not the last color
        // calculate min and max of x[c] for given x[0]..x[c-1]
        xmin = n - remaining[c];  if (xmin < 0) xmin = 0;
        xmax = m[c];  if (xmax > n) xmax = n;
        x0 = xm[c];  if (x0 < xmin) x0 = xmin;  if (x0 > xmax) x0 = xmax;
        // loop for all x[c] from mean and up
        for (x = x0, s2 = 0.; x <= xmax; x++) {
            xi[c] = x;
            sum += s1 = loop(n - x, c + 1);      // recursive loop for remaining colors
            if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
            s2 = s1;
        }
        // loop for all x[c] from mean and down
        for (x = x0 - 1; x >= xmin; x--) {
            xi[c] = x;
            sum += s1 = loop(n - x, c + 1);      // recursive loop for remaining colors
            if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
            s2 = s1;
        }
    }
    else {
        // last color
        xi[c] = n;
        s1 = probability(xi);
        for (i = 0; i < colors; i++) {
            sx[i] += s1 * xi[i];
            sxx[i] += s1 * xi[i] * xi[i];
        }
        sn++;
        sum = s1;
    }
    return sum;
}

void MultiWalleniusNCHypergeometric::reportError(const char* text) {
    // report an error message via global object
    errors.reportError(text);
}



/***********************************************************************
Univariate Fishers noncentral hypergeometric distribution
***********************************************************************/
int32_t RandomVariates::fishersNCHyp(int32_t n, int32_t m, int32_t N, double odds) {
    /*
    This function generates a random variate with Fisher's noncentral
    hypergeometric distribution.

    This distribution simulates selection when animal fates are independent,
    where the survival rate is fixed.

    This function uses inversion by chop-down search from zero when parameters
    are small, and the ratio-of-uniforms rejection method when the former
    method would be too slow or would give overflow.
    */
    int32_t fak, addd;                 // used for undoing transformations
    int32_t x;                         // result

    // check if parameters are valid
    if (n > N || m > N || n < 0 || m < 0 || odds <= 0.) {
        if (odds == 0.) {
            if (n > N - m) reportError("Not enough items with nonzero weight in function fishersNCHyp");
            return 0;
        }
        reportError("Parameter out of range in function fishersNCHyp");
    }

    if (odds == 1.) {
        // use hypergeometric function if odds == 1
        return hypergeometric(n, m, N);
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
    if (n == 0 || odds == 0.) return addd;

    if (fak == -1) {
        // reciprocal odds if inverting
        odds = 1. / odds;
    }

    // choose method
    if (n < 30 && N < 1024 && odds > 1.E-5 && odds < 1.E5) {
        // use inversion by chop down method
        x = fishersNCHypInversion(n, m, N, odds);
    }
    else {
        // use ratio-of-uniforms method
        x = fishersNCHypRatioOfUnifoms(n, m, N, odds);
    }

    // undo symmetry transformations  
    return x * fak + addd;
}

// Subfunctions used by FishersNCHyp

int32_t RandomVariates::fishersNCHypInversion(int32_t n, int32_t m, int32_t N, double odds) {
    /*
    Subfunction for fishersNCHyp distribution.
    Implements Fisher's noncentral hypergeometric distribution by inversion
    method, using chop-down search starting at zero.

    Valid only for 0 <= n <= m <= N/2.
    Without overflow check the parameters must be limited to n < 30, N < 1024,
    and 1.E-5 < odds < 1.E5. This limitation is acceptable because this method
    is slow for higher n.

    The execution time of this function grows with n.
    */
    int32_t x;                         // x value
    int32_t ll;                        // derived parameter
    double f;                          // scaled function value 
    double sum;                        // scaled sum of function values
    double a1, a2, b1, b2, f1, f2;     // factors in recursive calculation
    double u;                          // uniform random variate

    ll = N - m - n;

    if (n != fnc_n_last || m != fnc_m_last || N != fnc_N_last || odds != fnc_o_last) {
        // parameters have changed. set-up
        fnc_n_last = n; fnc_m_last = m; fnc_N_last = N; fnc_o_last = odds;

        // f(0) is set to an arbitrary value because it cancels out.
        // A low value is chosen to avoid overflow.
        fnc_f0 = 1.E-100;

        // calculate summation of e(x), using the formula:
        // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
        // All divisions are avoided by scaling the parameters
        sum = f = fnc_f0;  fnc_scale = 1.;
        a1 = m;  a2 = n;  b1 = 1;  b2 = ll + 1;
        for (x = 1; x <= n; x++) {
            f1 = a1 * a2 * odds;
            f2 = b1 * b2;
            a1--;  a2--;  b1++;  b2++;
            f *= f1;
            sum *= f2;
            fnc_scale *= f2;
            sum += f;
            // overflow check. not needed if parameters are limited:
            // if (sum > 1E100) {sum *= 1E-100; f *= 1E-100; fnc_scale *= 1E-100;}
        }
        fnc_f0 *= fnc_scale;
        fnc_scale = sum;
        // now f(0) = fnc_f0 / fnc_scale.
        // We are still avoiding all divisions by saving the scale factor
    }

    // uniform random
    u = random() * fnc_scale;

    // recursive calculation:
    // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
    f = fnc_f0;  x = 0;  a1 = m;  a2 = n;  b1 = 0;  b2 = ll;
    do {
        u -= f;
        if (u <= 0) break;
        x++;  b1++;  b2++;
        f *= a1 * a2 * odds;
        u *= b1 * b2;
        // overflow check. not needed if parameters are limited:
        // if (u > 1.E100) {u *= 1E-100;  f *= 1E-100;}
        a1--;  a2--;
    } while (x < n);
    return x;
}

int32_t RandomVariates::fishersNCHypRatioOfUnifoms(int32_t n, int32_t m, int32_t N, double odds) {
    /*
    Subfunction for FishersNCHyp distribution.
    Valid for 0 <= n <= m <= N/2, odds != 1

    Fisher's noncentral hypergeometric distribution by ratio-of-uniforms
    rejection method.

    The execution time of this function is almost independent of the parameters.
    */
    int32_t ll;                        // N-m-n
    int32_t mode;                      // mode
    double mean;                       // mean
    double variance;                   // variance
    double x;                          // real sample
    int32_t k;                         // integer sample
    double u;                          // uniform random
    double lf;                         // ln(f(x))
    double aa, bb, g1, g2;             // temporary

    ll = N - m - n;

    if (n != fnc_n_last || m != fnc_m_last || N != fnc_N_last || odds != fnc_o_last) {
        // parameters have changed. set-up
        fnc_n_last = n;  fnc_m_last = m;  fnc_N_last = N;  fnc_o_last = odds;

        // find approximate mean
        aa = (m + n) * odds + ll; bb = sqrt(aa * aa - 4 * odds * (odds - 1) * m * n);
        mean = (aa - bb) / (2 * (odds - 1));

        // find approximate variance
        aa = mean * (m - mean); bb = (n - mean) * (mean + ll);
        variance = N * aa * bb / ((N - 1) * (m * bb + (n + ll) * aa));

        // compute log(odds)
        fnc_logb = log(odds);

        // find center and width of hat function
        fnc_a = mean + 0.5;
        fnc_h = 1.028 + 1.717 * sqrt(variance + 0.5) + 0.032 * fabs(fnc_logb);

        // find safety bound
        fnc_bound = (int32_t)(mean + 4.0 * fnc_h);
        if (fnc_bound > n) fnc_bound = n;

        // find mode
        mode = (int32_t)(mean);
        g1 = (double)(m - mode) * (n - mode) * odds;
        g2 = (double)(mode + 1) * (ll + mode + 1);
        if (g1 > g2 && mode < n) mode++;

        // value at mode to scale with:
        fnc_lfm = mode * fnc_logb - fc_lnpk(mode, ll, m, n);
    }

    while (true) {
        u = random();
        if (u == 0) continue;                      // avoid divide by 0
        x = fnc_a + fnc_h * (random() - 0.5) / u;
        if (x < 0. || x > 2E9) continue;           // reject, avoid overflow
        k = (int32_t)(x);                          // truncate
        if (k > fnc_bound) continue;               // reject if outside safety bound
        lf = k * fnc_logb - fc_lnpk(k, ll, m, n) - fnc_lfm;// compute function value
        if (u * (4.0 - u) - 3.0 <= lf) break;      // lower squeeze accept
        if (u * (u - lf) > 1.0) continue;            // upper squeeze reject
        if (2.0 * log(u) <= lf) break;             // final acceptance
    }
    return k;
}

/***********************************************************************
Members of class FishersNCHypergeometric, for calculating probability etc.
***********************************************************************/

FishersNCHypergeometric::FishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy) {
    // constructor
    // set parameters
    this->n = n;  this->m = m;  this->N = N;
    this->odds = odds;  this->accuracy = accuracy;

    // check validity of parameters
    if (n < 0 || m < 0 || N < 0 || odds < 0. || n > N || m > N) {
        reportError("Parameter out of range in class FishersNCHypergeometric");
    }
    if (accuracy < 0) accuracy = 0;
    if (accuracy > 1) accuracy = 1;
    // initialize
    logodds = log(odds);  scale = rsum = 0.;
    parametersChanged = 1;
    // calculate xmin and xmax
    xmin = m + n - N;  if (xmin < 0) xmin = 0;
    xmax = n;  if (xmax > m) xmax = m;
}

double FishersNCHypergeometric::probability(int32_t x) {
    // calculate probability function
    const double accur = accuracy * 0.1;// accuracy of calculation

    if (x < xmin || x > xmax) return 0;

    if (n == 0) return 1.;

    if (odds == 1.) {
        // central hypergeometric
        return exp(
            lnFac(m) - lnFac(x) - lnFac(m - x) +
            lnFac(N - m) - lnFac(n - x) - lnFac((N - m) - (n - x)) -
            (lnFac(N) - lnFac(n) - lnFac(N - n)));
    }

    if (odds == 0.) {
        if (n > N - m) reportError("Not enough items with nonzero weight in FishersNCHypergeometric::probability");
        return x == 0;
    }

    if (!rsum) {
        // first time. calculate rsum = reciprocal of sum of proportional 
        // function over all probable x values
        int32_t x1, x2;                // x loop
        double y;                      // value of proportional function
        x1 = (int32_t)mean();          // start at mean
        if (x1 < xmin) x1 = xmin;
        x2 = x1 + 1;
        scale = 0.; scale = lng(x1);   // calculate scale to avoid overflow
        rsum = 1.;                     // = exp(lng(x1)) with this scale
        for (x1--; x1 >= xmin; x1--) {
            rsum += y = exp(lng(x1));  // sum from x1 and down 
            if (y < accur) break;      // until value becomes negligible
        }
        for (; x2 <= xmax; x2++) {     // sum from x2 and up
            rsum += y = exp(lng(x2));
            if (y < accur) break;      // until value becomes negligible
        }
        rsum = 1. / rsum;              // save reciprocal sum
    }
    return exp(lng(x)) * rsum;         // function value
}

double FishersNCHypergeometric::mean(void) {
    // Find approximate mean
    // Calculation analogous with mode
    double a, b;                       // temporaries in calculation
    double mean;                       // mean

    if (odds == 1.) {                  // simple hypergeometric
        return double(m) * n / N;
    }
    // calculate Cornfield mean
    a = (m + n) * odds + (N - m - n);
    b = a * a - 4. * odds * (odds - 1.) * m * n;
    b = b > 0. ? sqrt(b) : 0.;
    mean = (a - b) / (2. * (odds - 1.));
    return mean;
}

int32_t FishersNCHypergeometric::mode(void) {
    // Find mode (exact)
    // Uses the method of Liao and Rosen, The American Statistician, vol 55,
    // no 4, 2001, p. 366-369.
    // Note that there is an error in Liao and Rosen's formula. 
    // Replace sgn(b) with -1 in Liao and Rosen's formula. 

    double a, b, c, d;                 // coefficients for quadratic equation
    double x;                          // mode
    int32_t L = m + n - N;
    int32_t m1 = m + 1, n1 = n + 1;

    if (odds == 1.) {
        // simple hypergeometric
        x = (m + 1.) * (n + 1.) / (N + 2.);
    }
    else {
        // calculate analogously to Cornfield's mean
        a = 1. - odds;
        b = (m1 + n1) * odds - L;
        c = -(double)m1 * n1 * odds;
        d = b * b - 4 * a * c;
        d = d > 0. ? sqrt(d) : 0.;
        x = (d - b) / (a + a);
    }
    return (int32_t)x;
}

double FishersNCHypergeometric::moments(double * mean_, double * var_) {
    // calculate exact mean and variance
    // return value = sum of f(x), expected = 1.
    double y, sy = 0, sxy = 0, sxxy = 0, me1;
    int32_t x, xm, x1;
    const double accur = 0.1 * accuracy;         // accuracy of calculation
    xm = (int32_t)mean();                        // approximation to mean
    for (x = xm; x <= xmax; x++) {
        y = probability(x);
        x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
        sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
        if (y < accur && x != xm) break;
    }
    for (x = xm - 1; x >= xmin; x--) {
        y = probability(x);
        x1 = x - xm;  // subtract approximate mean to avoid loss of precision in sums
        sy += y; sxy += x1 * y; sxxy += x1 * x1 * y;
        if (y < accur) break;
    }
    me1 = sxy / sy;
    *mean_ = me1 + xm;
    y = sxxy / sy - me1 * me1;
    if (y < 0) y = 0;
    *var_ = y;
    return sy;
}

double FishersNCHypergeometric::lng(int32_t x) {
    // natural log of proportional function
    // returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!*odds^x)
    int32_t x2 = n - x, m2 = N - m;
    if (parametersChanged) {
        mFac = lnFac(m) + lnFac(m2);
        xLast = -99; parametersChanged = false;
    }
    if (m < FAK_LEN && m2 < FAK_LEN)  goto DEFLT;
    switch (x - xLast) {
    case 0:   // x unchanged
        break;
    case 1:   // x incremented. calculate from previous value
        xFac += log(double(x) * (m2 - x2) / (double(x2 + 1) * (m - x + 1)));
        break;
    case -1:  // x decremented. calculate from previous value
        xFac += log(double(x2) * (m - x) / (double(x + 1) * (m2 - x2 + 1)));
        break;
    default: DEFLT: // calculate all
        xFac = lnFac(x) + lnFac(x2) + lnFac(m - x) + lnFac(m2 - x2);
    }
    xLast = x;
    return mFac - xFac + x * logodds - scale;
}

double FishersNCHypergeometric::variance(void) {
    // find approximate variance (poor approximation)    
    double my = mean(); // approximate mean
    // find approximate variance from Fisher's noncentral hypergeometric approximation
    double r1 = my * (m - my); double r2 = (n - my) * (my + N - n - m);
    if (r1 <= 0. || r2 <= 0.) return 0.;
    double var = N * r1 * r2 / ((N - 1) * (m * r2 + (N - m) * r1));
    if (var < 0.) var = 0.;
    return var;
}

void FishersNCHypergeometric::reportError(const char* text) {
    // report an error message via global object
    errors.reportError(text);
}



/***********************************************************************
Multivariate Fishers noncentral hypergeometric distribution
***********************************************************************/
void RandomVariates::multiFishersNCHyp(int32_t * destination,
    int32_t * source, double * weights, int32_t n, int colors) {
    /*
    This function generates a vector of random variates with the
    multivariate Fisher's noncentral hypergeometric distribution.

    Parameters:
    destination:    An output array to receive the number of balls of each
    color. Must have space for at least 'colors' elements.
    source:         An input array containing the number of balls of each
    color in the urn. Must have 'colors' elements.
    All elements must be non-negative.
    weights:        The odds of each color. Must have 'colors' elements.
    All elements must be non-negative.
    n:              The number of balls drawn from the urn.
    Can't exceed the total number of balls with nonzero weight
    in the urn.
    colors:         The number of possible colors.

    Method: The conditional method is used for generating a sample with the
    approximate distribution. This sample is used as a starting point for
    a Gibbs sampler. The accuracy depends on the number of scans with the
    Gibbs sampler.

    The function will reduce the number of colors, if possible, by eliminating
    colors with zero weight or zero number and pooling together colors with the
    same weight. A symmetry transformation is used if more than half the balls
    are taken. The problem thus reduced is handled in the arrays osource,
    oweights and osample of dimension colors2.
    */
    int order1[MAXCOLORS];             // sort order, index into source and destination
    int order2[MAXCOLORS];             // corresponding index into osource when equal weights pooled together
    int order3[MAXCOLORS];             // secondary index for sorting by variance
    int32_t osource[MAXCOLORS];        // contents of source, sorted by weight with equal weights pooled together
    int32_t osample[MAXCOLORS];        // balls sampled, sorted by weight
    double oweights[MAXCOLORS];        // sorted list of weights
    double var[MAXCOLORS];             // sorted list of variance
    int32_t x = 0;                     // univariate sample
    int32_t m;                         // number of items of one color
    int32_t m1, m2;                    // number of items in each weight group
    int32_t msum;                      // total number of items of several or all colors
    int32_t n0;                        // remaining balls to sample
    int32_t n1, n2;                    // sample size for each weight group
    double w = 0.;                     // weight or variance of items of one color
    double w1, w2;                     // mean weight of each weight group  
    double wsum;                       // total weight of all items of several or all colors
    double odds;                       // weight ratio
    int i, j, k;                       // loop counters
    int a, b;                          // limits for weight group
    int c, c1, c2;                     // color index
    int colors2;                       // reduced number of colors, number of entries in osource
    int ngibbs;                        // number of scans in Gibbs sampler
    int invert = 0;                    // 1 if symmetry transformation used

    // check validity of parameters
    if (n < 0 || colors < 0 || colors > MAXCOLORS) reportError("Parameter out of range in function multiFishersNCHyp");
    if (colors == 0) return;
    if (n == 0) {
        for (i = 0; i < colors; i++) destination[i] = 0; return;
    }

    // check validity of array parameters
    for (i = 0, msum = 0; i < colors; i++) {
        m = source[i];  w = weights[i];
        if (m < 0 || w < 0) reportError("Parameter negative in function multiFishersNCHyp");
        if (w) msum += m;
    }

    // sort by weight, heaviest first
    for (i = 0; i < colors; i++) order1[i] = order3[i] = i;
    for (i = 0; i < colors - 1; i++) {
        c = order1[i];  k = i;
        w = weights[c];  if (source[c] == 0) w = 0;
        for (j = i + 1; j < colors; j++) {
            c2 = order1[j];
            if (weights[c2] > w && source[c2]) {
                w = weights[c2];  k = j;
            }
        }
        order1[i] = order1[k];  order1[k] = c;
    }

    // Skip any items with zero weight
    // this solves all problems with zero weights
    while (colors && (weights[c = order1[colors - 1]] == 0 || source[c] == 0)) {
        colors--;  destination[c] = 0;
    }

    // check if we are taking all, or too many, balls
    if (n >= msum) {
        if (n > msum) reportError("Taking more items than there are in function multiFishersNCHyp");
        for (i = 0; i < colors; i++) { c = order1[i];  destination[c] = source[c]; }
        return;
    }

    if (n > msum / 2) {
        // improve accuracy by symmetry transformation
        for (i = 0, j = colors - 1; i < j; i++, j--) { // reverse order list
            c = order1[i];  order1[i] = order1[j];  order1[j] = c;
        }
        n = msum - n;  invert = 1;
    }

    // copy source and weights into ordered lists and pool together colors with same weight
    for (i = 0, c2 = -1; i < colors; i++) {
        c = order1[i];
        if (i == 0 || weights[c] != w) {
            c2++;
            x = source[c];
            oweights[c2] = w = invert ? 1. / weights[c] : weights[c];
        }
        else {
            x += source[c];
        }
        osource[c2] = x;
        order2[i] = c2;
        osample[c2] = 0;
    }
    colors2 = c2 + 1;

    // simple cases  
    if (colors2 == 1) osample[0] = n;
    if (colors2 == 2) {
        x = fishersNCHyp(n, osource[0], msum, oweights[0] / oweights[1]);
        osample[0] = x;  osample[1] = n - x;
    }

    if (colors2 > 2) {
        // divide weights into two groups, heavy and light
        a = 0;  b = colors2 - 1;
        w = sqrt(oweights[0] * oweights[colors2 - 1]);
        do {
            c = (a + b) / 2;
            if (oweights[c] > w) a = c; else b = c;
        } while (b > a + 1);
        a = 0; // heavy group goes from a to b-1, light group goes from b to colors2-1

        // calculate mean weight for heavy group
        for (i = a, m1 = 0, wsum = 0; i < b; i++) {
            m1 += osource[i];  wsum += oweights[i] * osource[i];
        }
        w1 = wsum / m1;

        // calculate mean weight for light group
        for (i = b, m2 = 0, wsum = 0; i < colors2; i++) {
            m2 += osource[i];  wsum += oweights[i] * osource[i];
        }
        w2 = wsum / m2;

        // split sample n into heavy (n1) and light (n2) groups
        n1 = fishersNCHyp(n, m1, m1 + m2, w1 / w2);
        n2 = n - n1;
        n0 = n1;

        // loop twice, for the two groops
        for (k = 0; k < 2; k++) {

            // split group into single colors by calling FishersNCHyp b-a-1 times
            for (i = a; i < b - 1; i++) {
                m = osource[i];  w = oweights[i];

                // calculate mean weight of remaining colors
                for (j = i + 1, msum = 0, wsum = 0; j < b; j++) {
                    m1 = osource[j];  w1 = oweights[j];
                    msum += m1;  wsum += m1 * w1;
                }

                // split out color i
                if (w == w1) {
                    x = hypergeometric(n0, m, msum + m);
                }
                else {
                    if (wsum == 0) {
                        x = n0;
                    }
                    else {
                        odds = w * msum / wsum;
                        x = fishersNCHyp(n0, m, msum + m, odds);
                    }
                }
                osample[i] += x;
                n0 -= x;
            }

            // get the last color in the group
            osample[i] += n0;

            // set parameters for second group
            a = b;  b = colors2;  n0 = n2;
        }

        // calculate variance
        MultiFishersNCHypergeometric(n, osource, oweights, colors2).variance(var);

        // sort again, this time by variance
        for (i = 0; i < colors2 - 1; i++) {
            c = order3[i];  k = i;
            w = var[c];
            for (j = i + 1; j < colors2; j++) {
                c2 = order3[j];
                if (var[c2] > w) {
                    w = var[c2];  k = j;
                }
            }
            order3[i] = order3[k];  order3[k] = c;
        }

        // determine number of scans (not fine-tuned):
        ngibbs = 4;  if (accuracy < 1E-6) ngibbs = 6;  if (colors2 > 5) ngibbs++;

        // Gibbs sampler
        for (k = 0; k < ngibbs; k++) {
            for (i = 0; i < colors2; i++) {
                c1 = order3[i];
                j = i + 1;  if (j == colors2) j = 0;
                c2 = order3[j];
                n1 = osample[c1] + osample[c2];
                x = fishersNCHyp(n1, osource[c1], osource[c1] + osource[c2], oweights[c1] / oweights[c2]);
                osample[c1] = x;
                osample[c2] = n1 - x;
            }
        }
    }

    if (invert) {
        // reverse symmetry transformation on result
        for (i = 0; i < colors2; i++) {
            osample[i] = osource[i] - osample[i];
        }
    }

    // un-sort sample into destination
    for (i = 0; i < colors; i++) {
        c1 = order1[i];  c2 = order2[i];
        if (source[c1] == osource[c2]) {
            destination[c1] = osample[c2];
        }
        else {
            x = hypergeometric(osample[c2], source[c1], osource[c2]);
            destination[c1] = x;
            osample[c2] -= x;
            osource[c2] -= source[c1];
        }
    }
}

// calculation methods in class CMultiFishersNCHypergeometric

MultiFishersNCHypergeometric::MultiFishersNCHypergeometric(int32_t n_, int32_t * m_, double * odds_, int colors_, double accuracy_) {
    // constructor
    int i;                             // loop counter

    // copy parameters
    n = n_;  colors = colors_;  accuracy = accuracy_;

    // check if parameters are valid
    reduced = 2;  N = Nu = 0;  usedcolors = 0;
    for (i = 0; i < colors; i++) {
        nonzero[i] = 1;                // remember if color i has m > 0 and odds > 0
        m[usedcolors] = m_[i];         // copy m
        N += m_[i];                    // sum of m
        if (m_[i] <= 0) {
            nonzero[i] = 0;            // color i unused
            reduced |= 1;
            if (m_[i] < 0) reportError("Parameter m negative in constructor for MultiFishersNCHypergeometric");
        }
        odds[usedcolors] = odds_[i];   // copy odds
        if (odds_[i] <= 0) {
            nonzero[i] = 0;            // color i unused
            reduced |= 1;
            if (odds_[i] < 0) reportError("Parameter odds negative in constructor for MultiFishersNCHypergeometric");
        }
        if (usedcolors > 0 && nonzero[i] && odds[usedcolors] != odds[usedcolors - 1]) {
            reduced &= ~2;             // odds are not all equal
        }
        if (nonzero[i]) {
            Nu += m[usedcolors];       // sum of m for used colors
            usedcolors++;              // skip color i if zero
        }
    }
    if (N < n)  reportError("Taking more items than there are in constructor for MultiFishersNCHypergeometric");
    if (Nu < n) reportError("Not enough items with nonzero weight in constructor for MultiFishersNCHypergeometric");

    // calculate mFac and logodds
    for (i = 0, mFac = 0.; i < usedcolors; i++) {
        mFac += lnFac(m[i]);
        logodds[i] = log(odds[i]);
    }
    // initialize
    sn = 0;
}

void MultiFishersNCHypergeometric::mean(double * mu) {
    // calculates approximate mean of multivariate Fisher's noncentral
    // hypergeometric distribution. Result is returned in mu[0..colors-1].
    // The calculation is reasonably fast.
    int i, j;                           // color index
    double mur[MAXCOLORS];              // mean for used colors

    // get mean of used colors
    mean1(mur);

    // resolve unused colors
    for (i = j = 0; i < colors; i++) {
        if (nonzero[i]) {
            mu[i] = mur[j++];
        }
        else {
            mu[i] = 0.;
        }
    }
}

void MultiFishersNCHypergeometric::mean1(double * mu) {
    // calculates approximate mean of multivariate Fisher's noncentral
    // hypergeometric distribution, except for unused colors
    double r, r1;                       // iteration variable
    double q;                           // mean of color i
    double W;                           // total weight
    int i;                              // color index
    int iter = 0;                       // iteration counter

    if (usedcolors < 3) {
        // simple cases
        if (usedcolors == 1) mu[0] = n;
        if (usedcolors == 2) {
            mu[0] = FishersNCHypergeometric(n, m[0], Nu, odds[0] / odds[1]).mean();
            mu[1] = n - mu[0];
        }
    }
    else if (n == Nu) {
        // Taking all balls
        for (i = 0; i < usedcolors; i++) mu[i] = m[i];
    }
    else {
        // not a special case

        // initial guess for r
        for (i = 0, W = 0.; i < usedcolors; i++) W += m[i] * odds[i];
        r = (double)n * Nu / ((Nu - n) * W);

        if (r > 0.) {
            // iteration loop to find r
            do {
                r1 = r;
                for (i = 0, q = 0.; i < usedcolors; i++) {
                    q += m[i] * r * odds[i] / (r * odds[i] + 1.);
                }
                r *= n * (Nu - q) / (q * (Nu - n));
                if (++iter > 100) reportError("convergence problem in function MultiFishersNCHypergeometric::mean");
            } while (fabs(r - r1) > 1E-5);
        }

        // get result
        for (i = 0; i < usedcolors; i++) {
            mu[i] = m[i] * r * odds[i] / (r * odds[i] + 1.);
        }
    }
}

void MultiFishersNCHypergeometric::variance(double * var, double * mean_) {
    // calculates approximate variance of multivariate Fisher's noncentral
    // hypergeometric distribution (accuracy is not too good).
    // Variance is returned in variance[0..colors-1].
    // Mean is returned in mean_[0..colors-1] if not NULL.
    // The calculation is reasonably fast.
    double r1, r2;
    double mu[MAXCOLORS];
    int i, j;

    mean1(mu);         // Mean of used colors

    for (i = j = 0; i < colors; i++) {
        if (nonzero[i]) {
            r1 = mu[j] * (m[j] - mu[j]);
            r2 = (n - mu[j]) * (mu[j] + Nu - n - m[j]);
            if (r1 <= 0. || r2 <= 0.) {
                var[i] = 0.;
            }
            else {
                var[i] = Nu * r1 * r2 / ((Nu - 1) * (m[j] * r2 + (Nu - m[j]) * r1));
            }
            j++;
        }
        else {  // unused color
            var[i] = 0.;
        }
    }

    // Store mean if mean_ is not NULL
    if (mean_) {
        // resolve unused colors
        for (i = j = 0; i < colors; i++) {
            if (nonzero[i]) {
                mean_[i] = mu[j++];
            }
            else {
                mean_[i] = 0.;
            }
        }
    }
}

double MultiFishersNCHypergeometric::probability(int32_t * x) {
    // Calculate probability function.
    // Note: The first-time call takes very long time because it requires
    // a calculation of all possible x combinations with probability >
    // accuracy, which may be extreme.
    // The calculation uses logarithms to avoid overflow. 
    // (Recursive calculation may be faster, but this has not been implemented)
    int i, j;                             // color index
    int32_t xsum = 0;                     // sum of x
    int32_t Xu[MAXCOLORS];                // x for used colors

    // resolve unused colors
    for (i = j = 0; i < colors; i++) {
        if (nonzero[i]) {
            Xu[j++] = x[i];               // copy x to array of used colors
            xsum += x[i];                 // sum of x
        }
        else {
            if (x[i]) return 0.;          // taking balls with zero weight
        }
    }

    if (xsum != n) {
        reportError("sum of x values not equal to n in function MultiFishersNCHypergeometric::probability");
    }

    for (i = 0; i < usedcolors; i++) {
        if (Xu[i] > m[i] || Xu[i] < 0 || Xu[i] < n - Nu + m[i]) return 0.;  // Outside bounds for x
    }

    if (n == 0 || n == Nu) return 1.;   // deterministic cases

    if (usedcolors < 3) {               // cases with < 3 colors
        if (usedcolors < 2) return 1.;
        // Univariate probability
        return FishersNCHypergeometric(n, m[0], Nu, odds[0] / odds[1], accuracy).probability(Xu[0]);
    }

    if (reduced & 2) {
        // All odds are equal. This is multivariate central hypergeometric distribution
        int32_t sx = n, sm = N;
        double p = 1.;
        for (i = 0; i < usedcolors - 1; i++) {
            // Use univariate hypergeometric (usedcolors-1) times
            p *= FishersNCHypergeometric(sx, m[i], sm, 1.).probability(x[i]);
            sx -= x[i];  sm -= m[i];
        }
        return p;
    }

    // all special cases eliminated. Calculate sum of all function values
    if (sn == 0) sumOfAll();            // first time initialize

    return exp(lng(Xu)) * rsum;         // function value
}

double MultiFishersNCHypergeometric::moments(double * mean, double * variance, int32_t * combinations) {
    // calculates mean and variance of the Fisher's noncentral hypergeometric 
    // distribution by calculating all combinations of x-values with
    // probability > accuracy.
    // Return value = 1.
    // Returns the mean in mean[0...colors-1]
    // Returns the variance in variance[0...colors-1]

    int i, j;                           // color index
    if (sn == 0) {
        // first time initialization includes calculation of mean and variance
        sumOfAll();
    }
    // copy results and resolve unused colors
    for (i = j = 0; i < colors; i++) {
        if (nonzero[i]) {
            mean[i] = sx[j];
            variance[i] = sxx[j];
            j++;
        }
        else {
            mean[i] = variance[i] = 0.;
        }
    }
    if (combinations) *combinations = sn;
    return 1.;
}

void MultiFishersNCHypergeometric::sumOfAll() {
    // this function does the very time consuming job of calculating the sum
    // of the proportional function g(x) over all possible combinations of
    // the x[i] values with probability > accuracy. These combinations are 
    // generated by the recursive function loop().
    // The mean and variance are generated as by-products.

    int i;                             // color index
    int32_t msum;                      // sum of m[i]

    // get approximate mean
    mean1(sx);

    // round mean to integers
    for (i = 0, msum = 0; i < usedcolors; i++) {
        msum += xm[i] = (int32_t)(sx[i] + 0.4999999);
    }
    // adjust truncated x values to make the sum = n
    msum -= n;
    for (i = 0; msum < 0; i++) {
        if (xm[i] < m[i]) {
            xm[i]++; msum++;
        }
    }
    for (i = 0; msum > 0; i++) {
        if (xm[i] > 0) {
            xm[i]--; msum--;
        }
    }

    // adjust scale factor to g(mean) to avoid overflow
    scale = 0.; scale = lng(xm);

    // initialize for recursive loops
    sn = 0;
    for (i = usedcolors - 1, msum = 0; i >= 0; i--) {
        remaining[i] = msum;  msum += m[i];
    }
    for (i = 0; i < usedcolors; i++) {
        sx[i] = 0;  sxx[i] = 0;
    }

    // recursive loops to calculate sums of g(x) over all x combinations
    rsum = 1. / loop(n, 0);

    // calculate mean and variance
    for (i = 0; i < usedcolors; i++) {
        sxx[i] = sxx[i] * rsum - sx[i] * sx[i] * rsum * rsum;
        sx[i] = sx[i] * rsum;
    }
}

double MultiFishersNCHypergeometric::loop(int32_t n, int c) {
    // recursive function to loop through all combinations of x-values.
    // used by SumOfAll
    int32_t x, x0;                     // x of color c
    int32_t xmin, xmax;                // min and max of x[c]
    double s1, s2, sum = 0.;           // sum of g(x) values
    int i;                             // loop counter

    if (c < usedcolors - 1) {
        // not the last color
        // calculate min and max of x[c] for given x[0]..x[c-1]
        xmin = n - remaining[c];  if (xmin < 0) xmin = 0;
        xmax = m[c]; if (xmax > n) xmax = n;
        x0 = xm[c];  if (x0 < xmin) x0 = xmin;  if (x0 > xmax) x0 = xmax;
        // loop for all x[c] from mean and up
        for (x = x0, s2 = 0.; x <= xmax; x++) {
            xi[c] = x;
            sum += s1 = loop(n - x, c + 1);      // recursive loop for remaining colors
            if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
            s2 = s1;
        }
        // loop for all x[c] from mean and down
        for (x = x0 - 1; x >= xmin; x--) {
            xi[c] = x;
            sum += s1 = loop(n - x, c + 1);      // recursive loop for remaining colors
            if (s1 < accuracy && s1 < s2) break; // stop when values become negligible
            s2 = s1;
        }
    }
    else {
        // last color
        xi[c] = n;
        // sums and squaresums    
        s1 = exp(lng(xi));                       // proportional function g(x)
        for (i = 0; i < usedcolors; i++) {       // update sums
            sx[i] += s1 * xi[i];
            sxx[i] += s1 * xi[i] * xi[i];
        }
        sn++;
        sum += s1;
    }
    return sum;
}

double MultiFishersNCHypergeometric::lng(int32_t * x) {
    // natural log of proportional function g(x)
    double y = 0.;
    int i;
    for (i = 0; i < usedcolors; i++) {
        y += x[i] * logodds[i] - lnFac(x[i]) - lnFac(m[i] - x[i]);
    }
    return mFac + y - scale;
}

void MultiFishersNCHypergeometric::reportError(const char* text) {
    // report an error message via global object
    errors.reportError(text);
}


/*****************************************************************************
Tables of residues of a certain expansion of the error function.
These tables are used in the Laplace method for calculating Wallenius noncentral
hypergeometric distribution. Used in WalleniusNCHypergeometric::laplace() and
MultiWalleniusNCHypergeometric::laplace().

These tables are generated by ERFRESMK.CPP from https://www.agner.org/random/stocc.zip
Please see the file ERFRESMK.CPP for a detailed description.
You must re-run ERFRESMK.CPP if the following constants are changed.

The following constants have been used for making the tables below:
ERFRES_B =   16    (-log2 of lowest precision)
ERFRES_E =   40    (-log2 of highest precision)
ERFRES_S =    2    (step size from begin to end)
ERFRES_N =   13    (number of tables)
ERFRES_L =   48    (length of each table)
*****************************************************************************/

//number of standard deviations to integrate
double numSDev[ERFRES_N] = {
    4.324919041, 4.621231001, 4.900964208, 5.16657812, 5.419983175, 5.662697617, 5.895951217, 6.120756286, 6.337957755, 6.548269368, 6.752300431, 6.950575948, 7.143552034 };

//tables of error function residues
double erfRes[ERFRES_N][ERFRES_L] = {
    // 0: precision 1.53E-05
    {1.77242680540608204400E+00, 4.42974050453076994800E-01, 5.52683719287987914000E-02, 4.57346771067359261300E-03,
    2.80459064155823224600E-04, 1.34636065677244878500E-05, 5.21352785817798300800E-07, 1.65832271688171705300E-08,
    4.38865717471213472100E-10, 9.76518286165874680600E-12, 1.84433013221606645200E-13, 2.98319658966723379900E-15,
    4.16751049288581722800E-17, 5.06844293411881381200E-19, 5.40629927341885830200E-21, 5.09268600245963099700E-23,
    4.26365286677037947600E-25, 3.19120961809492396300E-27, 2.14691825888024309100E-29, 1.30473994083903636000E-31,
    7.19567933922698314600E-34, 3.61655672748362805300E-36, 1.66299275803871018000E-38, 7.02143932105206679000E-41,
    2.73122271211734530800E-43, 9.81824938600123102500E-46, 3.27125155121613401700E-48, 1.01290491600297417870E-50,
    2.92208589554240568800E-53, 7.87247562929246970200E-56, 1.98510836143160618600E-58, 4.69476368999432417500E-61,
    1.04339442450396263710E-63, 2.18317315734482557700E-66, 4.30811606197931495800E-69, 8.03081062303437395000E-72,
    1.41637813978528824300E-74, 2.36693694351427741600E-77, 3.75309000199992425400E-80, 5.65409397708564003600E-83,
    8.10322084538751956300E-86, 1.10610328893385430400E-88, 1.43971150303803736000E-91, 1.78884532267880002700E-94,
    2.12393968173898899400E-97, 2.41222807417272408400E-100, 2.62311608532487946600E-103, 2.73362126618952541200E-106},
    // 1: precision 3.81E-06
    {1.77244708953065753100E+00, 4.43074113723358004800E-01, 5.53507546366094128100E-02, 4.60063583541917741200E-03,
    2.85265530531727983900E-04, 1.39934570721569428400E-05, 5.61234181715130108200E-07, 1.87635216633109792000E-08,
    5.29386567604284238200E-10, 1.27170893476994027400E-11, 2.62062404027629145800E-13, 4.66479837413316034000E-15,
    7.22069968938298529400E-17, 9.78297384753513147400E-19, 1.16744590415498861200E-20, 1.23448081765041655900E-22,
    1.16327347874717650400E-24, 9.82084801488552519700E-27, 7.46543820883360082800E-29, 5.13361419796185362400E-31,
    3.20726459674397306300E-33, 1.82784782995019591600E-35, 9.53819678596992509200E-38, 4.57327699736894183000E-40,
    2.02131302843758583500E-42, 8.26035836048709995200E-45, 3.13004443753993537100E-47, 1.10264466279388735400E-49,
    3.62016356599029098800E-52, 1.11028768672354227000E-54, 3.18789098809699663200E-57, 8.58660896411902915800E-60,
    2.17384332055877431800E-62, 5.18219413865915035000E-65, 1.16526530012222654600E-67, 2.47552943408735877700E-70,
    4.97637013794934320200E-73, 9.47966949394160838200E-76, 1.71361124212171341900E-78, 2.94335699587741039100E-81,
    4.80983789654609513600E-84, 7.48676877660738410200E-87, 1.11129798477201315100E-89, 1.57475145101473103400E-92,
    2.13251069867015016100E-95, 2.76249093386952224300E-98, 3.42653604413897348900E-101, 4.07334940102519697800E-104},
    // 2: precision 9.54E-07
    {1.77245216056180140300E+00, 4.43102496776356791100E-01, 5.53772601883593673800E-02, 4.61054749828262358400E-03,
    2.87253302758514987700E-04, 1.42417784632842086400E-05, 5.82408831964509309600E-07, 2.00745450404117050700E-08,
    5.91011604093749423400E-10, 1.49916022838813094600E-11, 3.29741365965300606900E-13, 6.32307780683001018100E-15,
    1.06252674842175897800E-16, 1.57257431560311360800E-18, 2.06034642322747725700E-20, 2.40159615347654528000E-22,
    2.50271435589313449400E-24, 2.34271631492982176000E-26, 1.97869636045309031700E-28, 1.51440731538936707000E-30,
    1.05452976534458622500E-32, 6.70612854853490875900E-35, 3.90863249061728208500E-37, 2.09490406980039604000E-39,
    1.03572639732910843160E-41, 4.73737271771599553200E-44, 2.01016799853191990700E-46, 7.93316727009805559200E-49,
    2.91896910080597410900E-51, 1.00361556207253403120E-53, 3.23138481735358914000E-56, 9.76266225260763484100E-59,
    2.77288342251948021500E-61, 7.41751660051554639600E-64, 1.87191699537047863600E-66, 4.46389809367038823800E-69,
    1.00740435367143552990E-71, 2.15468537440631290200E-74, 4.37372804933525238000E-77, 8.43676369508201162800E-80,
    1.54845094802349484100E-82, 2.70727577941653793200E-85, 4.51412388960109772800E-88, 7.18605932463221426200E-91,
    1.09328719452457957600E-93, 1.59123500193816486400E-96, 2.21770259794482485600E-99, 2.96235081914900644200E-102},
    // 3: precision 2.38E-07
    {1.77245342831958737100E+00, 4.43110438095780200600E-01, 5.53855581791170228000E-02, 4.61401880234106439000E-03,
    2.88031928895194049600E-04, 1.43505456256023050800E-05, 5.92777558091362167400E-07, 2.07920891418090254000E-08,
    6.28701715960960909000E-10, 1.65457546101845217200E-11, 3.81394501062348919800E-13, 7.73640169798996619200E-15,
    1.38648618664047143200E-16, 2.20377376795474051600E-18, 3.11871105901085320300E-20, 3.94509797765438339700E-22,
    4.47871054279593642800E-24, 4.58134444141001287500E-26, 4.23915369932833545200E-28, 3.56174643985755223000E-30,
    2.72729562179570597400E-32, 1.90986605998546816600E-34, 1.22720072734085613700E-36, 7.25829034260272865500E-39,
    3.96321699645874596800E-41, 2.00342049456074966200E-43, 9.40055798441764717800E-46, 4.10462275003981738400E-48,
    1.67166813346582579800E-50, 6.36422340874443565900E-53, 2.26969100679582421400E-55, 7.59750937838053600600E-58,
    2.39149482673471882600E-60, 7.09134153544718378800E-63, 1.98415128824311335000E-65, 5.24683837588056156800E-68,
    1.31326161465641387500E-70, 3.11571024962460536800E-73, 7.01627137211411880000E-76, 1.50162731270605666400E-78,
    3.05816530510335364700E-81, 5.93355048535012188600E-84, 1.09802441010335521600E-86, 1.94008240128183308800E-89,
    3.27631821921541675800E-92, 5.29343480369738200400E-95, 8.19001419434114020600E-98, 1.21456436757992622700E-100},
    // 4: precision 5.96E-08
    {1.77245374525903386300E+00, 4.43112635580628681700E-01, 5.53880993417431935600E-02, 4.61519508177347361400E-03,
    2.88323830371235781500E-04, 1.43956506488931199600E-05, 5.97533121516696046900E-07, 2.11560073234896927000E-08,
    6.49836113541376862800E-10, 1.75091216044688314800E-11, 4.16782737060155846600E-13, 8.80643257335436424800E-15,
    1.65748420791207225100E-16, 2.78707349086274968000E-18, 4.19899868515935354900E-20, 5.68498078698629510200E-22,
    6.93816222596422139400E-24, 7.65747618996655475200E-26, 7.66779861336649418200E-28, 6.98905143723583695400E-30,
    5.81737537190421990800E-32, 4.43568540037466870600E-34, 3.10768227888207447300E-36, 2.00640852664381818400E-38,
    1.19706367104711013300E-40, 6.61729939738396217600E-43, 3.39784063694262711800E-45, 1.62450416252839296200E-47,
    7.24798161653719932800E-50, 3.02428684730111423300E-52, 1.18255348374176440700E-54, 4.34156802253088795200E-57,
    1.49931575039307549400E-59, 4.87879082698754128200E-62, 1.49836511723882777600E-64, 4.34998243416684050900E-67,
    1.19554618884894856000E-69, 3.11506828608539767000E-72, 7.70504604851319512900E-75, 1.81153231245726529100E-77,
    4.05332288179748454100E-80, 8.64127160751002389800E-83, 1.75723563299790750600E-85, 3.41217779987510142000E-88,
    6.33324341504830543600E-91, 1.12470466360665277900E-93, 1.91282818505057981800E-96, 3.11838272111119088500E-99},
    // 5: precision 1.49E-08
    {1.77245382449389548700E+00, 4.43113238150016054000E-01, 5.53888635367372804600E-02, 4.61558298326459057200E-03,
    2.88429374592283566800E-04, 1.44135302457832808700E-05, 5.99599530816354110000E-07, 2.13293263207088596800E-08,
    6.60866899904610148200E-10, 1.80600922150303605400E-11, 4.38957621672449876700E-13, 9.54096365498724593600E-15,
    1.86125270560486321400E-16, 3.26743200260750243300E-18, 5.17322947745786073000E-20, 7.40303709577309752000E-22,
    9.59703297362487960100E-24, 1.12979041959758568400E-25, 1.21090586780714120800E-27, 1.18477600671972569200E-29,
    1.06110784945102789800E-31, 8.72301430014194580800E-34, 6.59978694597213862400E-36, 4.60782503988683505400E-38,
    2.97629996764696360400E-40, 1.78296967476668997800E-42, 9.92947813649120231300E-45, 5.15238281451496107200E-47,
    2.49648080941516617600E-49, 1.13183145876711695200E-51, 4.81083885812771760200E-54, 1.92068525483444959800E-56,
    7.21538203720691761200E-59, 2.55484244329461795400E-61, 8.54021947322263940200E-64, 2.69922457940407460300E-66,
    8.07806757099831088400E-69, 2.29233505413233278200E-71, 6.17627451352383776600E-74, 1.58198519435517862400E-76,
    3.85682833066898009900E-79, 8.96007783937447061800E-82, 1.98575880907873828900E-84, 4.20275001914011054200E-87,
    8.50301055680340658200E-90, 1.64613519849643900900E-92, 3.05222294684008316300E-95, 5.42516704506242119200E-98},
    // 6: precision 3.73E-09
    {1.77245384430261089200E+00, 4.43113402125597019200E-01, 5.53890898808651020700E-02, 4.61570802060252211600E-03,
    2.88466397094702578100E-04, 1.44203545983349722400E-05, 6.00457657669759309400E-07, 2.14076280553580130200E-08,
    6.66287908992827087900E-10, 1.83546080772263722600E-11, 4.51849203153760888400E-13, 1.00053478654150626250E-14,
    2.00133542358651377800E-16, 3.62647881190865840300E-18, 5.96489800325831839200E-20, 8.92069144951359438200E-22,
    1.21499978844978062400E-23, 1.50969159775091919100E-25, 1.71458470816131592700E-27, 1.78354149193378771000E-29,
    1.70298947555869630200E-31, 1.49600537831395400600E-33, 1.21186208172570666700E-35, 9.07362642179266008600E-38,
    6.29382543478586469600E-40, 4.05352760000606626000E-42, 2.42933889358226154400E-44, 1.35768914148821438100E-46,
    7.09017160688256911600E-49, 3.46664168532600651800E-51, 1.58991153690202909500E-53, 6.85218984466549798200E-56,
    2.77986852228382907500E-58, 1.06333492956411188200E-60, 3.84102521375678317000E-63, 1.31221496031384552800E-65,
    4.24584095965170648000E-68, 1.30291378525223696900E-70, 3.79687911940099574200E-73, 1.05205378465263412500E-75,
    2.77502269989758744900E-78, 6.97601832816401403200E-81, 1.67315109709482392200E-83, 3.83268665565667928900E-86,
    8.39358376033290752000E-89, 1.75907817494562062400E-91, 3.53115954806899335200E-94, 6.79562013989671425000E-97},
    // 7: precision 9.31E-10
    {1.77245384925478974400E+00, 4.43113446460012284000E-01, 5.53891560601252504200E-02, 4.61574755288994634700E-03,
    2.88479053368568788400E-04, 1.44228769021976818600E-05, 6.00800544645992949800E-07, 2.14414502554089331400E-08,
    6.68819005926294320800E-10, 1.85032367193584636900E-11, 4.58880445172944815400E-13, 1.02790650461108873560E-14,
    2.09055796622121955200E-16, 3.87357904265687446300E-18, 6.55355746022352119400E-20, 1.01398465283490267200E-21,
    1.43654532753298842400E-23, 1.86580454392148962200E-25, 2.22454554378132065200E-27, 2.43828788210971585600E-29,
    2.46099438567553070000E-31, 2.29136593939231572900E-33, 1.97178483051357608300E-35, 1.57129911859150760300E-37,
    1.16187715309016251400E-39, 7.98791034830625946600E-42, 5.11610271388176540200E-44, 3.05861085454619325800E-46,
    1.71006575230074253400E-48, 8.95787473757552059200E-51, 4.40426750636187741200E-53, 2.03593329808165663200E-55,
    8.86319619094250260800E-58, 3.63949556302483252000E-60, 1.41180525527432472100E-62, 5.18110448656726197600E-65,
    1.80130976146235507900E-67, 5.94089489436009998000E-70, 1.86108901096460881000E-72, 5.54453617603266634800E-75,
    1.57273231131712670500E-77, 4.25229555550383344000E-80, 1.09708064410784368000E-82, 2.70363777400980301400E-85,
    6.37064773173804957600E-88, 1.43666982549400138800E-90, 3.10359876850474266200E-93, 6.42822304267944541900E-96},
    // 8: precision 2.33E-10
    {1.77245385049283445600E+00, 4.43113458380306853400E-01, 5.53891751960330686200E-02, 4.61575984524613369300E-03,
    2.88483285115404915700E-04, 1.44237837119469849000E-05, 6.00933085215778545800E-07, 2.14555059613473259000E-08,
    6.69949807134525424700E-10, 1.85746173246056176400E-11, 4.62510251141501895600E-13, 1.04309449728125451550E-14,
    2.14376794695367282400E-16, 4.03195345507914206800E-18, 6.95901230873262760600E-20, 1.10422005968960415700E-21,
    1.61274044622451622200E-23, 2.17010646570190394600E-25, 2.69272585719737993500E-27, 3.08406442023150341400E-29,
    3.26412756902204044100E-31, 3.19659762892894327800E-33, 2.90079234489442113000E-35, 2.44307440922101839900E-37,
    1.91280099578638699700E-39, 1.39463784147443818800E-41, 9.48568383329895892700E-44, 6.02906080392955580400E-46,
    3.58720420688290561300E-48, 2.00136767763554841800E-50, 1.04877885428425423540E-52, 5.17045929753308956200E-55,
    2.40183088534749939500E-57, 1.05288434613857573000E-59, 4.36191374659545444200E-62, 1.71017740178796946700E-64,
    6.35417287308090154000E-67, 2.24023617204667066100E-69, 7.50388817892399787300E-72, 2.39087016939309798700E-74,
    7.25439736654156264700E-77, 2.09846227207024494800E-79, 5.79315651373498761100E-82, 1.52786617607871741100E-84,
    3.85332605389629328300E-87, 9.30196261538477647000E-90, 2.15126632809118648300E-92, 4.77058936290696223500E-95},
    // 9: precision 5.82E-11
    {1.77245385080234563500E+00, 4.43113461569894215700E-01, 5.53891806760746538300E-02, 4.61576361260268991600E-03,
    2.88484673044866409200E-04, 1.44241019771415521500E-05, 6.00982861902849871600E-07, 2.14611541966231908200E-08,
    6.70435999307504633400E-10, 1.86074527008731886600E-11, 4.64296589104966284700E-13, 1.05109058078120195880E-14,
    2.17373506425627932200E-16, 4.12736258800510237200E-18, 7.22027572389545573000E-20, 1.16641031427122158000E-21,
    1.74261574594878846800E-23, 2.40999131874158664000E-25, 3.08741471404781296800E-27, 3.66622899027160893300E-29,
    4.03832398444680182100E-31, 4.12964092806000764200E-33, 3.92459969957984993300E-35, 3.47023698321199047400E-37,
    2.85870037656881575800E-39, 2.19701222983622897200E-41, 1.57757442199878062800E-43, 1.05998290283581317870E-45,
    6.67461794578944750100E-48, 3.94493775265477963400E-50, 2.19180590286711897200E-52, 1.14647284342367091100E-54,
    5.65409064942635909000E-57, 2.63281413190197920300E-59, 1.15914855705146421000E-61, 4.83173813806023163900E-64,
    1.90931412007029721900E-66, 7.16152712238209948300E-69, 2.55277823724126351900E-71, 8.65775632882397637500E-74,
    2.79685049229469435800E-76, 8.61535752145576873700E-79, 2.53319381071928112300E-81, 7.11686161831786026200E-84,
    1.91227899461300469000E-86, 4.91879425560043181900E-89, 1.21226578717106016000E-91, 2.86511260628508142200E-94},
    // 10: precision 1.46E-11
    {1.77245385087972342800E+00, 4.43113462419744630200E-01, 5.53891822321947835700E-02, 4.61576475266972634100E-03,
    2.88485120632836570100E-04, 1.44242113476668549100E-05, 6.01001089101483108200E-07, 2.14633579957941871400E-08,
    6.70638121912630560800E-10, 1.86219965341716152100E-11, 4.65139560168398521100E-13, 1.05511053035457485150E-14,
    2.18978467579008781700E-16, 4.18179627467181890600E-18, 7.37905600609363562400E-20, 1.20666925770415139000E-21,
    1.83216676939141016100E-23, 2.58616160243870388400E-25, 3.39612594393133643000E-27, 4.15117456105401982300E-29,
    4.72512355800254106200E-31, 5.01108411105699264300E-33, 4.95452692086540934200E-35, 4.57052259669118191500E-37,
    3.93757613394119041600E-39, 3.17143225730425447800E-41, 2.39087136989889684400E-43, 1.68918677399352864600E-45,
    1.11992962513487784300E-47, 6.97720003652956407000E-50, 4.09017183052803247800E-52, 2.25925194899934230000E-54,
    1.17743902383784437300E-56, 5.79751618317805258800E-59, 2.70049127204827368400E-61, 1.19150157862632851000E-63,
    4.98581510751975724600E-66, 1.98102566456273457700E-68, 7.48277410614888503600E-71, 2.68994458637406843000E-73,
    9.21308680313745922900E-76, 3.00957175301701607000E-78, 9.38604174484261857600E-81, 2.79745691952436047200E-83,
    7.97548757616816228000E-86, 2.17700350714256603000E-88, 5.69442820814374326200E-91, 1.42855756885812751800E-93},
    // 11: precision 3.64E-12
    {1.77245385089906787700E+00, 4.43113462645337308000E-01, 5.53891826707801996000E-02, 4.61576509382801447000E-03,
    2.88485262834342722100E-04, 1.44242482379506758200E-05, 6.01007615943023924400E-07, 2.14641957411498484200E-08,
    6.70719685646245707700E-10, 1.86282265411023575000E-11, 4.65522856702499667400E-13, 1.05705070352080171380E-14,
    2.19800647930093079100E-16, 4.21139261151871749000E-18, 7.47068213693802656400E-20, 1.23132525686457329000E-21,
    1.89037080673535316000E-23, 2.70767450402634975900E-25, 3.62208731605653583200E-27, 4.52783644780645903400E-29,
    5.29116794891083221600E-31, 5.78191926529856774600E-33, 5.91019131357709915300E-35, 5.65375339320520942200E-37,
    5.06448494950527399600E-39, 4.25125004489814020300E-41, 3.34702040997479327500E-43, 2.47392597585772167100E-45,
    1.71856809642179370600E-47, 1.12329116466680264100E-49, 6.91635006957699099400E-52, 4.01648185933072044700E-54,
    2.20256743728563483200E-56, 1.14197705850825122000E-58, 5.60474946818590333800E-61, 2.60701847612354797700E-63,
    1.15061401831998511400E-65, 4.82402847794291118400E-68, 1.92339714685666953300E-70, 7.30092195189691915600E-73,
    2.64114863236683700200E-75, 9.11500639536260716600E-78, 3.00399043312000082200E-80, 9.46306767642663343000E-83,
    2.85205432245625504600E-85, 8.23120145271503093200E-88, 2.27678649791096140000E-90, 6.04082678746563674000E-93},
    // 12: precision 9.09E-13
    {1.77245385090390399000E+00, 4.43113462705021723200E-01, 5.53891827935733966800E-02, 4.61576519490408572200E-03,
    2.88485307416075940900E-04, 1.44242604760223605000E-05, 6.01009907022372119900E-07, 2.14645068933581115800E-08,
    6.70751738699247757000E-10, 1.86308168994678478700E-11, 4.65691470353760117700E-13, 1.05795367138350319200E-14,
    2.20205466324054638500E-16, 4.22680889851439179400E-18, 7.52117118137557251000E-20, 1.24569747014608843200E-21,
    1.92626007811754286900E-23, 2.78693040917777943300E-25, 3.77798094465194860200E-27, 4.80270052176922369800E-29,
    5.72806202403284098500E-31, 6.41118455649104110000E-33, 6.73530071235990996000E-35, 6.64287180769401900600E-37,
    6.15272463485746774200E-39, 5.35401292372264035500E-41, 4.37964050507321407500E-43, 3.37013878900376065400E-45,
    2.44151902553507999600E-47, 1.66674472552984171500E-49, 1.07324838386391679300E-51, 6.52532932562465070600E-54,
    3.75007759408864456600E-56, 2.03933010598440151000E-58, 1.05056269424470639500E-60, 5.13240427502016103000E-63,
    2.38044205354512290600E-65, 1.04929890842558070320E-67, 4.40052237815903136000E-70, 1.75760526644875492000E-72,
    6.69249991110777975200E-75, 2.43182093294000139800E-77, 8.44044451319186471300E-80, 2.80086205952805676200E-82,
    8.89407469263960473600E-85, 2.70501913533005623200E-87, 7.88617413146613817400E-90, 2.20568290007963387700E-92}
};
