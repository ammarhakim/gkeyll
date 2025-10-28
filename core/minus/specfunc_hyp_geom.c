#include <complex.h>
#include <math.h>
//#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <stdlib.h>

#include <specfunc_hyp_geom.h>

/*
 * to compile as a library (that can be used by specfunc.py):
 *
 *      cc -fPIC -c -O2 specfunc.c -o specfunc.o
 *      cc --shared specfunc.o -o libspecfunc.so 
 *
 * if the gamma-function stuff an dthe parabolic cylinder function is commented
 * in, you need
 *
 *      cc --shared specfunc.o -o libspecfunc.so -lgslcblas -lgsl
 *
 * instead
 *
 * hyp2f1 and hyp1f1 seem to work well, the parabolic cylinder functions pcfd 
 * are problematic for anything but small arguments
 *
 * felix, june 2014
 */

#define PREC_WARN_LIMIT 1e99 // warn if individual terms in the series are larger than this
// XXX 4. 12. 14: set to high value because contrary to what I thought I had 
// understood, there seem to be many perfectly fine cases where individual
// terms in the sum are higher than 10^15


// hypergeometric function for |z| < 1
double complex hyp2f1(const double complex a, const double complex b,
        const double complex c, const double complex z, sf_prms_and_info_t* p)
{
    double frac = p->tol;
    double complex n = 0;
    double complex summand = a*b/c * z;
    double complex sum = 1. + summand;

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (frac >= p->tol); n++)
    {
        summand *= (a+n-1)*(b+n-1)/(c+n-1) * z/n;
        p->prec_warning = p->prec_warning | (creal(summand) > PREC_WARN_LIMIT) 
                                          | (cimag(summand) > PREC_WARN_LIMIT);
        sum += summand;
        frac = fabs(creal(summand/sum)) + fabs(cimag((summand/sum)));
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = frac;
    return sum;
}

// confluent hypergeometric function
double complex hyp1f1(const double complex a, const double complex b, 
        const double complex z, sf_prms_and_info_t* p)
{
    double frac = p->tol;
    double complex n = 0;
    double complex summand = a/b * z;
    double complex sum = 1. + summand;
    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (frac >= p->tol) ; n++)
    {
        summand = summand/(b+n-1) * (a+n-1) * (z/n);
        p->prec_warning = p->prec_warning | (creal(summand) > PREC_WARN_LIMIT) 
                                          | (cimag(summand) > PREC_WARN_LIMIT);
        sum += summand;
        frac = fabs(creal(summand/sum)) + fabs(cimag((summand/sum)));
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = frac;
    return sum;
}

void hyp1f1_a_arr(const double complex* a, const double complex b, 
        const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a[i]/b * z;
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = 1./(b+n-1) * z / n; // same for all a
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f * (a[i]+n-1);
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}

void hyp1f1_b_arr(const double complex a, const double complex* b, 
        const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a/b[i] * z;
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = (a+n-1) * z / n; // same for all b
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f / (b[i]+n-1);
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}

void hyp1f1_z_arr(const double complex a, const double complex b, 
        const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a/b * z[i];
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = (a+n-1)/(b+n-1) / n; // same for all z
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f * z[i];
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}

void hyp1f1_all_arr(const double complex* a, const double complex* b, 
        const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    int i = 0;
    for (i = 0; i < n_arr; i++)
        out[i] = hyp1f1(a[i], b[i], z[i], p);
}

/*
double complex cgamma(const double complex z)
{
    gsl_sf_result lnabs = {0, 0};
    gsl_sf_result arg = {0, 0};
    gsl_sf_lngamma_complex_e(creal(z), cimag(z), &lnabs, &arg);
    return cexp(lnabs.val+1j*arg.val);
}

// parabolic cylinder function -- this is only ok for small parameters
double complex pcfd(const double complex nu, const double complex z,
        sf_prms_and_info_t* p)
{
    return cpow(2, nu/2)*cexp(-z*z/4) * csqrt(M_PI) * 
        (1./cgamma((1-nu)/2) * hyp1f1(-nu/2, 1./2, z*z/2, p)
            - csqrt(2.)*z/cgamma(-nu/2) 
            * hyp1f1((1-nu)/2, 3./2, z*z/2, p)); 
}

void pcfd_nu_arr(const double complex* nu, const double complex z, 
        double complex* out, size_t n_arr, sf_prms_and_info_t* p) 
{
    int i = 0;
    for (i = 0; i < n_arr; i++)
    {
        out[i] = pcfd(nu[i], z, p);
    }
}

void pcfd_z_arr(const double complex nu, const double complex* z, 
        double complex* out, size_t n_arr, sf_prms_and_info_t* p) 
{
    int i = 0;
    for (i = 0; i < n_arr; i++)
    {
        out[i] = pcfd(nu, z[i], p);
    }
}
*/

void hyp2f1_a_arr(const double complex* a, const double complex b, 
        const double complex c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a[i]*b/c * z;
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = (b+n-1)/(c+n-1) * z / n; // same for all a
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f * (a[i]+n-1);
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}

void hyp2f1_b_arr(const double complex a, const double complex* b, 
        const double complex c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    hyp2f1_a_arr(b, a, c, z, out, n_arr, p);
}

void hyp2f1_c_arr(const double complex a, const double complex b, 
        const double complex* c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a*b/c[i] * z;
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = (a+n-1)*(b+n-1) * z / n; // same for all a
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f / (c[i]+n-1);
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}


// hyp2d1 for an array of z values
void hyp2f1_z_arr(const double complex a, const double complex b, 
        const double complex c, const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    double complex f = 0.;
    double frac = 0.;
    double maxfrac = p->tol;
    double complex n = 0;
    double complex* summand = 
        (double complex*)malloc(sizeof(double complex) * n_arr);
    int i = 0;
    
    for (i = 0; i < n_arr; i++)
    {
        summand[i] = a*b/c * z[i];
        out[i] = 1. + summand[i];
    }

    p->prec_warning = 0;
    for (n = 2; (creal(n) < p->max_iter) && (maxfrac >= p->tol) ; n++)
    {
        maxfrac = 0;
        f = (a+n-1)*(b+n-1)/(c+n-1) / n; // same for all zs
        for (i = 0; i < n_arr; i++)
        {
            summand[i] *= f * z[i];
            p->prec_warning = p->prec_warning | (creal(summand[i]) > PREC_WARN_LIMIT) 
                                              | (cimag(summand[i]) > PREC_WARN_LIMIT);
            out[i] += summand[i];
            frac = fabs(creal(summand[i]/out[i])) 
                + fabs(cimag((summand[i]/out[i])));
            maxfrac = frac > maxfrac ? frac : maxfrac;            
        }
    }
    p->iters_needed = (int)creal(n);
    p->tol_achieved = maxfrac;
    free(summand);    
}

void hyp2f1_all_arr(const double complex* a, const double complex* b, 
        const double complex* c, const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p)
{
    int i = 0;
    for (i = 0; i < n_arr; i++)
        out[i] = hyp2f1(a[i], b[i], c[i], z[i], p);
}

