#include "SystemData.h" 
#include "math.h"
#include <iostream>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))

long int lsame_ (const char *ca, const char *cb) {
	int a = *(unsigned char *)ca;
	int b = *(unsigned char *)cb;
	return tolower(a) == tolower(b);
}

double pow_di (double *ap, long *bp) {
	double pow, x;
	long n;
	unsigned long u;

	pow = 1;
	x = *ap;
	n = *bp;

	if (n != 0) {
		if (n < 0) {
			n = -n;
			x = 1 / x;
		}
		for (u = n; ;) {
			if (u & 01) {  ///(u AND 01) == 1? 
				pow *= x;
			}
			if (u >>= 1) {
				x *= x;
			} else {
				break;
			}
		}
	}
	return (pow);
}


static int dlamc1_(long *beta, long *t, long *rnd, long *ieee1) {
	/* -- LAPACK auxiliary routine (version 3.0) -- Univ. of Tennessee, Univ.
	   of California Berkeley, NAG Ltd., Courant Institute, Argonne National
	   Lab, and Rice University October 31, 1992

	   Purpose =======
	   DLAMC1 determines the machine parameters given by BETA, T, RND, and
	   IEEE1.

	   Arguments =========

	   BETA (output) INTEGER The base of the machine.

	   T (output) INTEGER The number of ( BETA ) digits in the mantissa.

	   RND (output) LOGICAL Specifies whether proper rounding ( RND = .TRUE.
	   ) or chopping ( RND = .FALSE. ) occurs in addition. This may not

	   be a reliable guide to the way in which the machine performs

	   its arithmetic.

	   IEEE1 (output) LOGICAL Specifies whether rounding appears to be done
	   in the IEEE 'round to nearest' style.

	   Further Details ===============

	   The routine is based on the routine ENVRON by Malcolm and incorporates
	   suggestions by Gentleman and Marovich. See

	   Malcolm M. A. (1972) Algorithms to reveal properties of floating-point
	   arithmetic. Comms. of the ACM, 15, 949-951.

	   Gentleman W. M. and Marovich S. B. (1974) More on algorithms that
	   reveal properties of floating point arithmetic units. Comms. of the
	   ACM, 17, 276-277.

	   ===================================================================== */
	/* Initialized data */
	static long first = TRUE;

	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	static long lrnd;
	static double a, b, c, f;
	static long lbeta;
	static double savec;
	static long lieee1;
	static double t1, t2;
	static long lt;
	static double one, qtr;

	if (first) {
		first = FALSE;
		one = 1.;

		/* LBETA, LIEEE1, LT and LRND are the local values of BETA, IEEE1, T
		   and RND.

		   Throughout this routine we use the function DLAMC3 to ensure that
		   relevant values are stored and not held in registers, or are not
		   affected by optimizers. Compute a = 2.0**m with the smallest
		   positive integer m such that fl( a + 1.0 ) = a. */

		a = 1.;
		c = 1.;

		/* + WHILE( C.EQ.ONE )LOOP */
  L10:
		if (c == one) {
			a *= 2;
			c = dlamc3_(&a, &one);
			d__1 = -a;
			c = dlamc3_(&c, &d__1);
			goto L10;
		}
		/* + END WHILE

		   Now compute b = 2.0**m with the smallest positive integer m such
		   that fl( a + b ) .gt. a. */

		b = 1.;
		c = dlamc3_(&a, &b);

		/* + WHILE( C.EQ.A )LOOP */
  L20:
		if (c == a) {
			b *= 2;
			c = dlamc3_(&a, &b);
			goto L20;
		}
		/* + END WHILE

		   Now compute the base.  a and c are neighbouring floating point
		   numbers in the interval ( beta**t, beta**( t + 1 ) ) and so their
		   difference is beta. Adding 0.25 to c is to ensure that it is
		   truncated to beta and not ( beta - 1 ). */

		qtr = one / 4;
		savec = c;
		d__1 = -a;
		c = dlamc3_(&c, &d__1);
		lbeta = (long)(c + qtr);

		/* Now determine whether rounding or chopping occurs, by adding a bit
		   less than beta/2 and a bit more than beta/2 to a. */

		b = (double) lbeta;
		d__1 = b / 2;
		d__2 = -b / 100;
		f = dlamc3_(&d__1, &d__2);
		c = dlamc3_(&f, &a);
		if (c == a) 
			lrnd = TRUE;
	    else 
			lrnd = FALSE;
			
		d__1 = b / 2;
		d__2 = b / 100;
		f = dlamc3_(&d__1, &d__2);
		c = dlamc3_(&f, &a);
		if (lrnd && c == a) {
			lrnd = FALSE;
		}

		/* Try and decide whether rounding is done in the IEEE 'round to
		   nearest' style. B/2 is half a unit in the last place of the two
		   numbers A and SAVEC. Furthermore, A is even, i.e. has last bit
		   zero, and SAVEC is odd. Thus adding B/2 to A should not change A,
		   but adding B/2 to SAVEC should change SAVEC. */

		d__1 = b / 2;
		t1 = dlamc3_(&d__1, &a);
		d__1 = b / 2;
		t2 = dlamc3_(&d__1, &savec);
		lieee1 = t1 == a && t2 > savec && lrnd;

		/* Now find the mantissa, t. It should be the integer part of log to
		   the base beta of a, however it is safer to determine t by
		   powering.  So we find t as the smallest positive integer for which
		   fl( beta**t + 1.0 ) = 1.0. */

		lt = 0;
		a = 1.;
		c = 1.;

		/* + WHILE( C.EQ.ONE )LOOP */
  L30:
		if (c == one) {
			++lt;
			a *= lbeta;
			c = dlamc3_(&a, &one);
			d__1 = -a;
			c = dlamc3_(&c, &d__1);
			goto L30;
		}
		/* + END WHILE */
	}

	*beta = lbeta;
	*t = lt;
	*rnd = lrnd;
	*ieee1 = lieee1;
	return 0;
}								/* dlamc1_ */

static int dlamc2_(long *beta, long *t, long *rnd, double *eps, long *emin, double *rmin, long *emax, double *rmax) {
	/* -- LAPACK auxiliary routine (version 3.0) -- Univ. of Tennessee, Univ.
	   of California Berkeley, NAG Ltd., Courant Institute, Argonne National
	   Lab, and Rice University October 31, 1992
	   Purpose =======
	   DLAMC2 determines the machine parameters specified in its argument
	   list.
	   Arguments =========

	   BETA (output) INTEGER The base of the machine.

	   T (output) INTEGER The number of ( BETA ) digits in the mantissa.

	   RND (output) LOGICAL Specifies whether proper rounding ( RND = .TRUE.
	   ) or chopping ( RND = .FALSE. ) occurs in addition. This may not
	   be a reliable guide to the way in which the machine performs
	   its arithmetic.

	   EPS (output) DOUBLE PRECISION The smallest positive number such that
	   fl( 1.0 - EPS ) .LT. 1.0, where fl denotes the computed value.

	   EMIN (output) INTEGER The minimum exponent before (gradual) underflow
	   occurs.

	   RMIN (output) DOUBLE PRECISION The smallest normalized number for the
	   machine, given by BASE**( EMIN - 1 ), where BASE is the floating point
	   value of BETA.

	   EMAX (output) INTEGER The maximum exponent before overflow occurs.

	   RMAX (output) DOUBLE PRECISION The largest positive number for the
	   machine, given by BASE**EMAX * ( 1 - EPS ), where BASE is the floating
	   point value of BETA.

	   Further Details ===============

	   The computation of EPS is based on a routine PARANOIA by W. Kahan of
	   the University of California at Berkeley.

	   ===================================================================== */
	/* Table of constant values */
	/* Initialized data */
	static long first = TRUE;
	static long iwarn = FALSE;

	/* System generated locals */
	long i__1;
	double d__1, d__2, d__3, d__4, d__5;

	/* Builtin functions */
	/* Local variables */
	static long ieee;
	static double half;
	static long lrnd;
	static double leps, zero, a, b, c;
	static long i, lbeta;
	static double rbase;
	static long lemin, lemax, gnmin;
	static double smal;
	static long gpmin;
	static double third, lrmin, lrmax, sixth;
	static long lieee1;
	static long lt, ngnmin, ngpmin;
	static double one, two;

	if (first) {
		first = FALSE;
		zero = 0.;
		one = 1.;
		two = 2.;

		/* LBETA, LT, LRND, LEPS, LEMIN and LRMIN are the local values of
		   BETA, T, RND, EPS, EMIN and RMIN. Throughout this routine we use
		   the function DLAMC3 to ensure that relevant values are stored and
		   not held in registers, or are not affected by optimizers. DLAMC1
		   returns the parameters LBETA, LT, LRND and LIEEE1. */

		dlamc1_(&lbeta, &lt, &lrnd, &lieee1);          ///lbeta = -20.0, lt = 1, lrnd = FALSE,  lieee1 = 0(FALSE)

		/* Start to find EPS. */

		b = (double)lbeta;
		i__1 = -lt;
		a = pow_di(&b, &i__1);
		leps = a;

		/* Try some tricks to see whether or not this is the correct EPS. */
		b = two / 3;
		half = one / 2;
		d__1 = -half;
		sixth = dlamc3_(&b, &d__1);
		third = dlamc3_(&sixth, &sixth);
		d__1 = -half;
		b = dlamc3_(&third, &d__1);
		b = dlamc3_(&b, &sixth);
		b = fabs(b);
		if (b < leps) {
			b = leps;
		}
		leps = 1.;

		/* + WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
  L10:
		if (leps > b && b > zero) {
			leps = b;
			d__1 = half * leps;
			/* Computing 5th power */
			d__3 = two, d__4 = d__3, d__3 *= d__3;
			/* Computing 2nd power */
			d__5 = leps;
			d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
			c = dlamc3_(&d__1, &d__2);
			d__1 = -c;
			c = dlamc3_(&half, &d__1);
			b = dlamc3_(&half, &c);
			d__1 = -b;
			c = dlamc3_(&half, &d__1);
			b = dlamc3_(&half, &c);
			goto L10;
		}
		/* + END WHILE */

		if (a < leps) {
			leps = a;
		}

		/* Computation of EPS complete. Now find EMIN.  Let A = + or - 1, and
		   + or - (1 + BASE**(-3)). Keep dividing A by BETA until (gradual)
		   underflow occurs. This is detected when we cannot recover the
		   previous A. */

		rbase = one / lbeta;
		smal = one;
		for (i = 1; i <= 3; ++i) {
			d__1 = smal * rbase;
			smal = dlamc3_(&d__1, &zero);
			/* L20: */
		}
		a = dlamc3_(&one, &smal);
		dlamc4_(&ngpmin, &one, &lbeta);
		d__1 = -one;
		dlamc4_(&ngnmin, &d__1, &lbeta);
		dlamc4_(&gpmin, &a, &lbeta);
		d__1 = -a;
		dlamc4_(&gnmin, &d__1, &lbeta);
		ieee = FALSE;

		if (ngpmin == ngnmin && gpmin == gnmin) {
			if (ngpmin == gpmin) {
				lemin = ngpmin;
				/* ( Non twos-complement machines, no gradual underflow;
				   e.g., VAX ) */
			} else if (gpmin - ngpmin == 3) {
				lemin = ngpmin - 1 + lt;
				ieee = TRUE;
				/* ( Non twos-complement machines, with gradual underflow;
				   e.g., IEEE standard followers ) */
			} else {
				lemin = MIN (ngpmin, gpmin);
				/* ( A guess; no known machine ) */
				iwarn = TRUE;
			}

		} else if (ngpmin == gpmin && ngnmin == gnmin) {
			if ( (i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
				lemin = MAX (ngpmin, ngnmin);
				/* ( Twos-complement machines, no gradual underflow; e.g.,
				   CYBER 205 ) */
			} else {
				lemin = MIN (ngpmin, ngnmin);
				/* ( A guess; no known machine ) */
				iwarn = TRUE;
			}

		} else if ( (i__1 = ngpmin - ngnmin, abs (i__1)) == 1 && gpmin == gnmin) {
			if (gpmin - MIN (ngpmin, ngnmin) == 3) {
				lemin = MAX (ngpmin, ngnmin) - 1 + lt;
				/* ( Twos-complement machines with gradual underflow; no
				   known machine ) */
			} else {
				lemin = MIN (ngpmin, ngnmin);
				/* ( A guess; no known machine ) */
				iwarn = TRUE;
			}

		} else {
			/* Computing MIN */
			i__1 = MIN (ngpmin, ngnmin), i__1 = MIN (i__1, gpmin);
			lemin = MIN (i__1, gnmin);
			/* ( A guess; no known machine ) */
			iwarn = TRUE;
		}
		/* Comment out this if block if EMIN is ok */
		if (iwarn) {
			first = TRUE;
			std::cout<<"\n WARNING. The value EMIN may be incorrect:- EMIN = "<<lemin<<" after inspection, the value EMIN looks acceptable. please comment out \n the IF block as marked within the code of routine DLAMC2, \n otherwise supply EMIN explicitly."<<std::endl;
			std::cout<<"SystemData.cpp: Line: 210."<<std::endl;
		/****	Melder_warning (L"\n\n WARNING. The value EMIN may be incorrect:- " "EMIN = ", Melder_integer (lemin),
			                L"\nIf, after inspection, the value EMIN looks acceptable"
			                "please comment out \n the IF block as marked within the"
			                "code of routine DLAMC2, \n otherwise supply EMIN" "explicitly.\n");  ****/
		}
		/* ** Assume IEEE arithmetic if we found denormalised numbers above,
		   or if arithmetic seems to round in the IEEE style, determined in
		   routine DLAMC1. A true IEEE machine should have both things true;
		   however, faulty machines may have one or the other. */

		ieee = ieee || lieee1;

		/* Compute RMIN by successive division by BETA. We could compute RMIN
		   as BASE**( EMIN - 1 ), but some machines underflow during this
		   computation. */

		lrmin = 1.;
		i__1 = 1 - lemin;
		for (i = 1; i <= 1 - lemin; ++i) {
			d__1 = lrmin * rbase;
			lrmin = dlamc3_ (&d__1, &zero);
			/* L30: */
		}

		/* Finally, call DLAMC5 to compute EMAX and RMAX. */

		dlamc5_ (&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
	}

	*beta = lbeta;
	*t = lt;
	*rnd = lrnd;
	*eps = leps;
	*emin = lemin;
	*rmin = lrmin;
	*emax = lemax;
	*rmax = lrmax;

	return 0;
}							/* dlamc2_ */

static double dlamc3_ (double *a, double *b)
/* Purpose =======
   dlamc3_ is intended to force A and B to be stored prior to doing the
   addition of A and B , for use in situations where optimizers might hold
   one of these in a register.
   Arguments =========

   A, B (input) DOUBLE PRECISION The values A and B.
   ===================================================================== */
{
	volatile double ret_val;

	ret_val = *a + *b;
	return ret_val;
}	  /* dlamc3_ */


static int dlamc4_ (long *emin, double *start, long *base) {
	/* -- LAPACK auxiliary routine (version 2.0) -- Univ. of Tennessee, Univ.
	   of California Berkeley, NAG Ltd., Courant Institute, Argonne National
	   Lab, and Rice University October 31, 1992

	   Purpose =======

	   DLAMC4 is a service routine for DLAMC2.

	   Arguments =========

	   EMIN (output) EMIN The minimum exponent before (gradual) underflow,
	   computed by

	   setting A = START and dividing by BASE until the previous A can not be
	   recovered.

	   START (input) DOUBLE PRECISION The starting point for determining
	   EMIN.

	   BASE (input) INTEGER The base of the machine.

	   ===================================================================== */
	/* System generated locals */
	long i__1;
	double d__1;

	/* Local variables */
	static double zero, a;
	static long i;
	static double rbase, b1, b2, c1, c2, d1, d2;
	static double one;

	a = *start;
	one = 1.;
	rbase = one / *base;
	zero = 0.;
	*emin = 1;
	d__1 = a * rbase;
	b1 = dlamc3_ (&d__1, &zero);
	c1 = a;
	c2 = a;
	d1 = a;
	d2 = a;
	/* + WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. $ ( D1.EQ.A ).AND.( D2.EQ.A
	   ) )LOOP */
 L10:
	if (c1 == a && c2 == a && d1 == a && d2 == a) {
		-- (*emin);
		a = b1;
		d__1 = a / *base;
		b1 = dlamc3_ (&d__1, &zero);
		d__1 = b1 * *base;
		c1 = dlamc3_ (&d__1, &zero);
		d1 = zero;
		i__1 = *base;
		for (i = 1; i <= *base; ++i) {
			d1 += b1;
			/* L20: */
		}
		d__1 = a * rbase;
		b2 = dlamc3_ (&d__1, &zero);
		d__1 = b2 / rbase;
		c2 = dlamc3_ (&d__1, &zero);
		d2 = zero;
		i__1 = *base;
		for (i = 1; i <= *base; ++i) {
			d2 += b2;
			/* L30: */
		}
		goto L10;
	}
	/* + END WHILE */

	return 0;
}								/* dlamc4_ */


static int dlamc5_ (long *beta, long *p, long *emin, long *ieee, long *emax, double *rmax) {
	/*
	   First compute LEXP and UEXP, two powers of 2 that bound abs(EMIN). We
	   then assume that EMAX + abs(EMIN) will sum approximately to the bound
	   that is closest to abs(EMIN). (EMAX is the exponent of the required
	   number RMAX). */
	/* Table of constant values */
	static double c_b5 = 0.;

	/* System generated locals */
	long i__1;
	double d__1;

	/* Local variables */
	static long lexp;
	static double oldy;
	static long uexp, i;
	static double y, z;
	static long nbits;
	static double recbas;
	static long exbits, expsum, try__;

	lexp = 1;
	exbits = 1;
L10:
	try__ = lexp << 1;
	if (try__ <= - (*emin)) {
		lexp = try__;
		++exbits;
		goto L10;
	}
	if (lexp == - (*emin)) {
		uexp = lexp;
	} else {
		uexp = try__;
		++exbits;
	}

	/* Now -LEXP is less than or equal to EMIN, and -UEXP is greater than or
	   equal to EMIN. EXBITS is the number of bits needed to store the
	   exponent. */

	if (uexp + *emin > -lexp - *emin) 
		expsum = lexp << 1;
	else 
		expsum = uexp << 1;
	
	/* EXPSUM is the exponent range, approximately equal to EMAX - EMIN + 1 . */
	 
	*emax = expsum + *emin - 1;
	nbits = exbits + 1 + *p;

	/* NBITS is the total number of bits needed to store a floating-point number. */
	   
	if (nbits % 2 == 1 && *beta == 2) {
		/* Either there are an odd number of bits used to store a
		   floating-point number, which is unlikely, or some bits are not
		   used in the representation of numbers, which is possible, (e.g.
		   Cray machines) or the mantissa has an implicit bit, (e.g. IEEE
		   machines, Dec Vax machines), which is perhaps the most likely. We
		   have to assume the last alternative. If this is true, then we need
		   to reduce EMAX by one because there must be some way of
		   representing zero in an implicit-bit system. On machines like
		   Cray, we are reducing EMAX by one unnecessarily. */
		-- (*emax);
	}

	if (*ieee) {
		/* Assume we are on an IEEE machine which reserves one exponent for infinity and NaN. */	   
		-- (*emax);
	}

	/* Now create RMAX, the largest machine number, which should be equal to
	   (1.0 - BETA**(-P)) * BETA**EMAX . First compute 1.0 - BETA**(-P),
	   being careful that the result is less than 1.0 . */

	recbas = 1. / *beta;
	z = *beta - 1.;
	y = 0.;
	i__1 = *p;
	for (i = 1; i <= *p; ++i) {
		z *= recbas;
		if (y < 1.) 
			oldy = y;
		y = dlamc3_ (&y, &z);
		/* L20: */
	}
	if (y >= 1.) {
		y = oldy;
	}

	/* Now multiply by BETA**EMAX to get RMAX. */
	i__1 = *emax;
	for (i = 1; i <= *emax; ++i) {
		d__1 = y * (*beta);
		y = dlamc3_ (&d__1, &c_b5);
		/* L30: */
	}

	*rmax = y;
	return 0;
}								/* dlamc5_ */

double getSystemData(const char *cmach){
   	/* Initialized data */
	static long first = TRUE;

	/* System generated locals */
	long i__1;
	double ret_val;

	/* Builtin functions */
	/* Local variables */
	static double base;
	static long beta;
	static double emin, prec, emax;
	static long imin, imax;
	static long lrnd;
	static double rmin, rmax, t, rmach;
	static double smal, sfmin;
	static long it;
	static double rnd, eps;

	if (first) {
		first = FALSE;
		dlamc2_ (&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
		base = (double) beta;
		t = (double) it;
		if (lrnd) {
			rnd = 1.;
			i__1 = 1 - it;
			eps = pow_di (&base, &i__1) / 2;
		} else {
			rnd = 0.;
			i__1 = 1 - it;
			eps = pow_di (&base, &i__1);
		}
		prec = eps * base;
		emin = (double) imin;
		emax = (double) imax;
		sfmin = rmin;
		smal = 1. / rmax;
		if (smal >= sfmin) {
			/* Use smal plus a bit, to avoid the possibility of rounding
			   causing overflow when computing 1/sfmin. */
			sfmin = smal * (eps + 1.);
		}
	}
	
	if (lsame_ (cmach, "E")) {
		rmach = eps;
	} else if (lsame_ (cmach, "S")) {
		rmach = sfmin;
	} else if (lsame_ (cmach, "B")) {
		rmach = base;
	} else if (lsame_ (cmach, "P")) {
		rmach = prec;
	} else if (lsame_ (cmach, "N")) {
		rmach = t;
	} else if (lsame_ (cmach, "R")) {
		rmach = rnd;
	} else if (lsame_ (cmach, "M")) {
		rmach = emin;
	} else if (lsame_ (cmach, "U")) {
		rmach = rmin;
	} else if (lsame_ (cmach, "L")) {
		rmach = emax;
	} else if (lsame_ (cmach, "O")) {
		rmach = rmax;
	}

	ret_val = rmach;
	return ret_val;
}

#undef MAX
#undef MIN
