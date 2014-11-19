//SystemData.h
#ifndef _SYSTEMDATA_H_
#define _SYSTEMDATA_H_

#include "Structure.h"	

//SystemData.cpp
/**Returns true if ca[0] is the same letter as cb[0] regardless of case.**/
long int lsame_(const char *ca, const char *cb);

double pow_di (double *ap, long *bp);

static int dlamc1_ (long *beta, long *t, long *rnd, long *ieee1);

static int dlamc2_ (long *beta, long *t, long *rnd, double *eps, long *emin, double *rmin, long *emax, double *rmax);

static double dlamc3_ (double *a, double *b);

static int dlamc4_ (long *emin, double *start, long *base);

static int dlamc5_ (long *beta, long *p, long *emin, long *ieee, long *emax, double *rmax);

double getSystemData(const char *cmach);
/*  Purpose
    =======
    NUMblas_dlamch determines double machine parameters.
    Arguments
    =========
    CMACH   (input) char*
            Specifies the value to be returned by DLAMCH:
            = 'E' or 'e',   DLAMCH := eps
            = 'S' or 's ,   DLAMCH := sfmin
            = 'B' or 'b',   DLAMCH := base
            = 'P' or 'p',   DLAMCH := eps*base
            = 'N' or 'n',   DLAMCH := t
            = 'R' or 'r',   DLAMCH := rnd
            = 'M' or 'm',   DLAMCH := emin
            = 'U' or 'u',   DLAMCH := rmin
            = 'L' or 'l',   DLAMCH := emax
            = 'O' or 'o',   DLAMCH := rmax

            where

            eps   = relative machine precision
            sfmin = safe minimum, such that 1/sfmin does not overflow
            base  = base of the machine
            prec  = eps*base
            t     = number of (base) digits in the mantissa
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
            emin  = minimum exponent before (gradual) underflow
            rmin  = underflow threshold - base**(emin-1)
            emax  = largest exponent before overflow
            rmax  = overflow threshold  - (base**emax)*(1-eps)
   =====================================================================
*/
#endif
