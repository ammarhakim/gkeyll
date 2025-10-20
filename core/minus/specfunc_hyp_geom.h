#include <complex.h>
#include <stddef.h>

typedef struct sf_prms_and_info
{
    int max_iter;
    double tol;
    int iters_needed;
    double tol_achieved;
    int prec_warning;
} sf_prms_and_info_t;

double complex hyp2f1(const double complex a, const double complex b,
        const double complex c, const double complex z, sf_prms_and_info_t* p);

double complex hyp1f1(const double complex a, const double complex b, 
        const double complex z, sf_prms_and_info_t* p);

void hyp1f1_a_arr(const double complex* a, const double complex b, 
        const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp1f1_b_arr(const double complex a, const double complex* b, 
        const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp1f1_z_arr(const double complex a, const double complex b, 
        const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp1f1_all_arr(const double complex* a, const double complex* b, 
        const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

double complex cgamma(const double complex z);

double complex pcfd(const double complex nu, const double complex z,
        sf_prms_and_info_t* p);

void hyp2f1_a_arr(const double complex* a, const double complex b, 
        const double complex c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp2f1_b_arr(const double complex a, const double complex* b, 
        const double complex c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp2f1_c_arr(const double complex a, const double complex b, 
        const double complex* c, const double complex z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp2f1_z_arr(const double complex a, const double complex b, 
        const double complex c, const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void hyp2f1_all_arr(const double complex* a, const double complex* b, 
        const double complex* c, const double complex* z, double complex* out,
        size_t n_arr, sf_prms_and_info_t* p);

void pcfd_nu_arr(const double complex* nu, const double complex z, 
        double complex* out, size_t n_arr, sf_prms_and_info_t* p); 

void pcfd_z_arr(const double complex nu, const double complex* z, 
        double complex* out, size_t n_arr, sf_prms_and_info_t* p); 


