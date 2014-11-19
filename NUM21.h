#ifndef _NUM2_h_
#define _NUM2_h_

#include "Structure.h"


NUMfft_Table NUMfft_Table_create(long n = 0, double *trigcache = 0, long *splitcache = 0);
void NUMfft_Table_init (NUMfft_Table table, long n);

void NUMfft_forward (NUMfft_Table table, double *data);
/*
	Function:
		Calculates the Fourier Transform of a set of n real-valued data points.
		Replaces this data in array data [1...n] by the positive frequency half
		of its complex Fourier Transform, with a minus sign in the exponent.
	Preconditions:
		data != NULL;
		table must have been initialised with NUMfft_Table_init_f/d
	Postconditions:
		data[1] contains real valued first component (Direct Current)
		data[2..n-1] even index : real part; odd index: imaginary part of DFT.
		data[n] contains real valued last component (Nyquist frequency)

	Output parameters:

	data  r(1) = the sum from i=1 to i=n of r(i)

        If l =(int) (n+1)/2

          then for k = 2,...,l

             r(2*k-2) = the sum from i = 1 to i = n of

                  r(i)*cos((k-1)*(i-1)*2*pi/n)

             r(2*k-1) = the sum from i = 1 to i = n of

                 -r(i)*sin((k-1)*(i-1)*2*pi/n)

        if n is even

             r(n) = the sum from i = 1 to i = n of

                  (-1)**(i-1)*r(i)

		i.e., the ordering of the output array will be for n even
			r(1),(r(2),i(2)),(r(3),i(3)),...,(r(l-1),i(l-1)),r(l).
		Or ...., (r(l),i(l)) for n uneven.

 *****  note
             this transform is unnormalized since a call of NUMfft_forward
             followed by a call of NUMfft_backward will multiply the input
             sequence by n.
*/

void NUMfft_backward (NUMfft_Table table, double *data);
/*
	Function:
		Calculates the inverse transform of a complex array if it is the transform of real data.
		(Result in this case must be multiplied by 1/n.)
	Preconditions:
		n is an integer power of 2.
		data != NULL;
		data [1] contains real valued first component (Direct Current)
		data [2..n-1] even index : real part; odd index: imaginary part of DFT.
		data [n] contains real valued last component (Nyquist frequency)

		table must have been initialised with NUMfft_Table_init_f/d

	Output parameters

	data       for n even and for i = 1,...,n

             r(i) = r(1)+(-1)**(i-1)*r(n)

                  plus the sum from k=2 to k=n/2 of

                   2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

                  -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)

        for n odd and for i = 1,...,n

             r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of

                  2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

                 -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)

 *****  note
             this transform is unnormalized since a call of NUMfft_forward
             followed by a call of NUMfft_backward will multiply the input
             sequence by n.
*/

void NUMrealft (double *data, long n, int isign);

void NUMforwardRealFastFourierTransform (double *data, long n);
void NUMreverseRealFastFourierTransform (double *data, long n);

double NUMminimize_brent (double (*f) (double x, void *closure), double a, double b,
	                      void *closure, double tol, double *fx);
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

#endif
