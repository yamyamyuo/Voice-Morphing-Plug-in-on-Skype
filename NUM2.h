#include "NUM2.h"
#include "NUM.h"
#include "NUMfft_core.h"
#include "Structure.h"
#include "SystemData.h"
#include <iostream>

NUMfft_Table NUMfft_Table_create(long n, double *trigcache, long *splitcache)
{
    NUMfft_Table me = (NUMfft_Table)malloc(sizeof(structNUMfft_Table));
	if(!me)  return NULL;
	my n = n;
	my trigcache = trigcache;
	my splitcache = splitcache;
	return me;
}

void NUMfft_Table_init (NUMfft_Table me, long n) {
	my n = n;
	my trigcache = (double *)malloc(sizeof(double) * 3 * n);    /// NUMvector <double> (0, 3 * n - 1);
	my splitcache = (long *)malloc(sizeof(long) * 32);          /// NUMvector <long> (0, 31);
	NUMrffti (n, my trigcache, my splitcache);
}

void NUMfft_forward (NUMfft_Table me, double *data) {
	if (my n == 1) 
		return;
	drftf1 (my n, &data[1], my trigcache, my trigcache + my n, my splitcache);
}

void NUMfft_backward (NUMfft_Table me, double *data) {
	if (my n == 1) 
		return;
	drftb1 (my n, &data[1], my trigcache, my trigcache + my n, my splitcache);
}

void NUMrealft (double *data, long n, int isign) {
	isign == 1 ? NUMforwardRealFastFourierTransform (data, n) :
	NUMreverseRealFastFourierTransform (data, n);
}

void NUMforwardRealFastFourierTransform (double *data, long n) {
	NUMfft_Table table = NUMfft_Table_create();
	NUMfft_Table_init (table, n);
	NUMfft_forward (table, data);

	if (n > 1) {
		// To be compatible with old behaviour
		double tmp = data[n - 1];
		for (long i = n; i > 2; i--) {
			data[i] = data[i - 1];
		}
		data[2] = tmp;
	}
}

void NUMreverseRealFastFourierTransform (double *data, long n) {
	NUMfft_Table table = NUMfft_Table_create();

	if (n > 1) {
	 // To be compatible with old behaviour
		double tmp = data[1];
		for (long i = 2; i < n; i++) {
			data[i] = data[i + 1];
		}
		data[n] = tmp;
	}

	NUMfft_Table_init (table, n);
	NUMfft_backward (table, data);
}

double NUMminimize_brent (double (*f) (double x, void *closure), double a, 
                          double b, void *closure, double tol, double *fx) 
/*
	The function returns an estimate for the minimum location with accuracy
		3 * SQRT_EPSILON * abs(x) + tol.
	The function always obtains a local minimum which coincides with
	the global one only if a function under investigation being unimodular.
	If a function being examined possesses no local minimum within
	the given range, the function returns 'a' (if f(a) < f(b)), otherwise
	it returns the right range boundary value b.

	Algorithm

	The function makes use of the golden section procedure combined with
	parabolic interpolation.
	At every step, the program operates at three abscissae - x, v, and w.
	x - the last and the best approximation to the minimum location,
		i.e. f(x) <= f(a) or/and f(x) <= f(b)
		(if the function f has a local minimum in (a,b), then both
		conditions are fulfiled after one or two steps).
	v, w are previous approximations to the minimum location. They may
	coincide with a, b, or x (although the algorithm tries to make all
 	u, v, and w distinct). Points x, v, and w are used to construct
	interpolating parabola whose minimum will be treated as a new
	approximation to the minimum location if the former falls within
	[a,b] and reduces the range enveloping minimum more efficient than
	the golden section procedure.
	When f(x) has a second derivative positive at the minimum location
	(not coinciding with a or b) the procedure converges superlinearly
	at a rate order about 1.324
*/						  
{
	double x, v, fv, w, fw;
	const double golden = 1 - NUM_goldenSection;           /// 1 - 0.618...
	const double sqrt_epsilon = sqrt (getSystemData("Epsilon"));  /// Quantization error (relative machine precision)
	long itermax = 60;
    
	if(tol <= 0 || a >= b){
	    std::cout<<"tol = " <<tol<<" a = "<<a<<" b = "<<b<<" .Condition: tol > 0  && a < b."<<std::endl;
		std::cout<<"SoundCompute.cpp: Line 122."<<std::endl;
		exit(0);
	}    /**** Melder_assert (tol > 0 && a < b); ****/

	/* First step - golden section */
	v = a + golden * (b - a);
	fv = (*f) (v, closure);
	x = v;  w = v;
	*fx = fv;  fw = fv;

	for (long iter = 1; iter <= itermax; iter++) {
		double range = b - a;
		double middle_range = (a + b) / 2;
		double tol_act = sqrt_epsilon * fabs (x) + tol / 3;
		double new_step; /* Step at this iteration */

		if (fabs (x - middle_range) + range / 2 <= 2 * tol_act)
			return x;

		/* Obtain the golden section step */
		new_step = golden * (x < middle_range ? b - x : a - x);

		/* Decide if the parabolic interpolation can be tried	*/
		if (fabs (x - w) >= tol_act) {
			/*
				Interpolation step is calculated as p/q;
				division operation is delayed until last moment.
			*/
			double p, q, t;

			t = (x - w) * (*fx - fv);
			q = (x - v) * (*fx - fw);
			p = (x - v) * q - (x - w) * t;
			q = 2 * (q - t);

			if (q > 0) 
				p = -p;
		    else 
				q = -q;

			/*
				If x+p/q falls in [a,b], not too close to a and b,
				and isn't too large, it is accepted.
				If p/q is too large then the golden section procedure can
				reduce [a,b] range.
			*/
			if (fabs (p) < fabs (new_step * q) &&
			        p > q * (a - x + 2 * tol_act) &&
			        p < q * (b - x - 2 * tol_act)) {
				new_step = p / q;
			}
		}

		/* Adjust the step to be not less than tolerance. */
		if (fabs (new_step) < tol_act) {
			new_step = new_step > 0 ? tol_act : - tol_act;
		}
		/* Obtain the next approximation to min	and reduce the enveloping range */
		{
			double t = x + new_step;	/* Tentative point for the min	*/
			double ft = (*f) (t, closure);

			/*
				If t is a better approximation, reduce the range so that
				t would fall within it. If x remains the best, reduce the range
				so that x falls within it.
			*/

			if (ft <= *fx) {
				if (t < x) {
					b = x;
				} else {
					a = x;
				}

				v = w;  w = x;  x = t;
				fv = fw;  fw = *fx;  *fx = ft;
			} else {
				if (t < x) {
					a = t;
				} else {
					b = t;
				}

				if (ft <= fw || w == x) {
					v = w; w = t;
					fv = fw; fw = ft;
				} else if (ft <= fv || v == x || v == w) {
					v = t;
					fv = ft;
				}
			}
		}
	}
	std::cout<<"NUMminimize_brent: maximum number of iterations ("<< (int) itermax << ") exceeded.  Position:  SoundCompute.cpp: Line 228"<<std::endl;
	return x;
}
