#include "NUM.h"
#include "NUM2.h"
#include "Structure.h"
#include "time.h"

//NUM.cpp
#define NUM_interpolate_simple_cases \
	    if (nx < 1) return NUMundefined; \
	    if (x > nx) return y[nx]; \
	    if (x < 1) return y[1]; \
	    if (x == midleft) return y[midleft]; \
	    if (maxDepth > midright - 1) maxDepth = midright - 1; \
	    if (maxDepth > nx - midleft) maxDepth = nx - midleft; \
	    if (maxDepth <= NUM_VALUE_INTERPOLATE_NEAREST) return y[(long) floor (x + 0.5)]; \
	    if (maxDepth == NUM_VALUE_INTERPOLATE_LINEAR) return y[midleft] + (x - midleft) * (y[midright] - y[midleft]); \
	    if (maxDepth == NUM_VALUE_INTERPOLATE_CUBIC) { \
			  double yl = y[midleft], yr = y[midright]; \
			  double dyl = 0.5 * (yr - y[midleft - 1]), dyr = 0.5 * (y[midright + 1] - yl); \
		      double fil = x - midleft, fir = midright - x; \
			return yl * fir + yr * fil - fil * fir * (0.5 * (dyr - dyl) + (fil - 0.5) * (dyl + dyr - 2 * (yr - yl))); \
	}
 /***** 1 < x < nx && x not integer: interpolate.   line 5 *****/
 
//NUM.cpp 	369 / 403 ? alternative ro error   has choosen 369
double NUM_interpolate_sinc (double y [], long nx, double x, long maxDepth) {
	long ix, midleft = floor(x), midright = midleft + 1, left, right;
	double result = 0.0, a, halfsina, aa, daa, cosaa, sinaa, cosdaa, sindaa;
	
	NUM_interpolate_simple_cases
	
	left = midright - maxDepth, right = midleft + maxDepth;
	a = NUMpi * (x - midleft);
	halfsina = 0.5 * sin (a);
	aa = a / (x - left + 1); cosaa = cos (aa); sinaa = sin (aa);
	daa = NUMpi / (x - left + 1); cosdaa = cos (daa); sindaa = sin (daa);
	
	for (ix = midleft; ix >= left; ix --) {
		double d = halfsina / a * (1.0 + cosaa), help;
		result += y [ix] * d;
		a += NUMpi;
		help = cosaa * cosdaa - sinaa * sindaa;
		sinaa = cosaa * sindaa + sinaa * cosdaa;
		cosaa = help;
		halfsina = - halfsina;
	}

	a = NUMpi * (midright - x);
	halfsina = 0.5 * sin (a);
	aa = a / (right - x + 1); cosaa = cos (aa); sinaa = sin (aa);
	daa = NUMpi / (right - x + 1); cosdaa = cos (daa); sindaa = sin (daa);
	
	for (ix = midright; ix <= right; ix ++) {
		double d = halfsina / a * (1.0 + cosaa), help;
		result += y [ix] * d;
		a += NUMpi;
		help = cosaa * cosdaa - sinaa * sindaa;
		sinaa = cosaa * sindaa + sinaa * cosdaa;
		cosaa = help;
		halfsina = - halfsina;
	}

	return result;
}

struct improve_params {
	int depth;
	double *y;
	long ixmax;
	int isMaximum;
};
 

static double improve_evaluate (double x, void *closure) {
	struct improve_params *me = (struct improve_params *) closure;
	double y = NUM_interpolate_sinc (my y, my ixmax, x, my depth);
	return my isMaximum ? - y : y;
}

double NUMimproveExtremum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real, int isMaximum) {
	struct improve_params params;
	double result;
	if (ixmid <= 1) { *ixmid_real = 1; return y[1]; }
	if (ixmid >= nx) { *ixmid_real = nx; return y[nx]; } 
	if (interpolation <= NUM_PEAK_INTERPOLATE_NONE) { *ixmid_real = ixmid; return y[ixmid]; } 
	if (interpolation == NUM_PEAK_INTERPOLATE_PARABOLIC) {
		double dy = 0.5 * (y[ixmid + 1] - y[ixmid - 1]);
		double d2y = 2 * y[ixmid] - y[ixmid - 1] - y[ixmid + 1];
		*ixmid_real = ixmid + dy / d2y;
		return y[ixmid] + 0.5 * dy * dy / d2y; 
	}
	/* Sinc interpolation. */
	params.y = y;
	params.depth = interpolation == NUM_PEAK_INTERPOLATE_SINC70 ? 70 : 700;
	params.ixmax = nx;
	params.isMaximum = isMaximum;
	/*return isMaximum ?
		- NUM_minimize (ixmid - 1, ixmid, ixmid + 1, improve_evaluate, & params, 1e-10, 1e-11, ixmid_real) :
		  NUM_minimize (ixmid - 1, ixmid, ixmid + 1, improve_evaluate, & params, 1e-10, 1e-11, ixmid_real);*/
	*ixmid_real = NUMminimize_brent (improve_evaluate, ixmid - 1, ixmid + 1, & params, 1e-10, & result);
	return isMaximum ? - result : result;
}

double NUMimproveMaximum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real)
	{ return NUMimproveExtremum (y, nx, ixmid, interpolation, ixmid_real, 1); }

double NUMimproveMinimum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real)
	{ return NUMimproveExtremum (y, nx, ixmid, interpolation, ixmid_real, 0); }


void NUMrandomRestart (unsigned long seed) {
	/*
		Based on Knuth, p. 187,602.

		Knuth had:
			int s = seed;
		This is incorrect (even if an int is 32 bit), since a negative value causes a loop in the test
			if (s != 0)
				s >>= 1;
		because the >> operator on negative integers adds a sign bit to the left.
	 */
	long t, j;
	double u [2 * LONG_LAG - 1], ul [2 * LONG_LAG - 1];
	double ulp = 1.0 / (1L << 30) / (1L << 22), ss;
	ss = 2.0 * ulp * (seed + 2);   /* QUESTION: does this work if seed exceeds 2^32 - 3? See Knuth p. 187. */
	for (j = 0; j < LONG_LAG; j ++) {
		u [j] = ss;
		ul [j] = 0.0;
		ss += ss;
		if (ss >= 1.0) ss -= 1.0 - 2 * ulp;
	}
	for (; j < 2 * LONG_LAG - 1; j ++)
		u [j] = ul [j] = 0.0;
	u [1] += ulp;
	ul [1] = ulp;
	t = STREAM_SEPARATION - 1;
	while (t > 0) {
		for (j = LONG_LAG - 1; j > 0; j --) {
			ul [j + j] = ul [j];
			u [j + j] = u [j];
		}
		for (j = 2 * LONG_LAG - 2; j > LAG_DIFF; j -= 2) {
			ul [2 * LONG_LAG - 1 - j] = 0.0;
			u [2 * LONG_LAG - 1 - j] = u [j] - ul [j];
		}
		for (j = 2 * LONG_LAG - 2; j >= LONG_LAG; j --) if (ul [j] != 0) {
			ul [j - LAG_DIFF] = ulp - ul [j - LAG_DIFF];
			u [j - LAG_DIFF] += u [j];
			if (u [j - LAG_DIFF] >= 1.0) u [j - LAG_DIFF] -= 1.0;
			ul [j - LONG_LAG] = ulp - ul [j - LONG_LAG];
			u [j - LONG_LAG] += u [j];
			if (u [j - LONG_LAG] >= 1.0) u [j - LONG_LAG] -= 1.0;
		}
		if ((seed & 1) != 0) {
			for (j = LONG_LAG; j > 0; j --) {
				ul [j] = ul [j - 1];
				u [j] = u [j - 1];
			}
			ul [0] = ul [LONG_LAG];
			u [0] = u [LONG_LAG];
			if (ul [LONG_LAG] != 0) {
				ul [SHORT_LAG] = ulp - ul [SHORT_LAG];
				u [SHORT_LAG] += u [LONG_LAG];
				if (u [SHORT_LAG] >= 1.0) u [SHORT_LAG] -= 1.0;
			}
		}
		if (seed != 0) {
			seed >>= 1;
		} else {
			t --;
		}
	}
	for (j = 0; j < SHORT_LAG; j ++)
		randomArray [j + LAG_DIFF] = u [j];
	for (; j < LONG_LAG; j ++)
		randomArray [j - SHORT_LAG] = u [j];
	randomArrayPointer1 = 0;
	randomArrayPointer2 = LAG_DIFF;
	iquality = 0;
	randomInited = 1;
}

double NUMrandomFraction (void) {
	/*
		Knuth uses a long random array of length QUALITY to copy values from randomArray.
		We save 8 kilobytes by using randomArray as a cyclic array (10% speed loss).
	*/
	long p1, p2;
	double newValue;
	if (! randomInited) NUMrandomRestart (time (NULL));
	p1 = randomArrayPointer1, p2 = randomArrayPointer2;
	if (p1 >= LONG_LAG) p1 = 0;
	if (p2 >= LONG_LAG) p2 = 0;
	newValue = randomArray [p1] + randomArray [p2];
	if (newValue >= 1.0) newValue -= 1.0;
	randomArray [p1] = newValue;
	p1 ++;
	p2 ++;
	if (++ iquality == LONG_LAG) {
		for (; iquality < QUALITY; iquality ++) {
			double newValue2;
			/*
				Possible future minor speed improvement:
					the cyclic array is walked down instead of up.
					The following tests will then be for 0.
			*/
			if (p1 >= LONG_LAG) p1 = 0;
			if (p2 >= LONG_LAG) p2 = 0;
			newValue2 = randomArray [p1] + randomArray [p2];
			if (newValue2 >= 1.0) newValue2 -= 1.0;
			randomArray [p1] = newValue2;
			p1 ++;
			p2 ++;
		}
		iquality = 0;
	}
	randomArrayPointer1 = p1;
	randomArrayPointer2 = p2;
	return newValue;
}

double NUMrandomUniform (double lowest, double highest) {
	return lowest + (highest - lowest) * NUMrandomFraction ();
}

long NUMrandomInteger (long lowest, long highest) {
	return lowest + (long) ((highest - lowest + 1) * NUMrandomFraction ());
}

double NUMrandomGauss (double mean, double standardDeviation) {
	//	Knuth, p. 122.
	static int secondAvailable = 0;
	static double y;
	double s, x;
	if (secondAvailable) {
		secondAvailable = FALSE;
		return mean + standardDeviation * y;
	} else {
		repeat {
			x = 2.0 * NUMrandomFraction () - 1.0;   /* Inside the square [-1; 1] x [-1; 1]. */
			y = 2.0 * NUMrandomFraction () - 1.0;
			s = x * x + y * y;
		} until (s < 1.0);   /* Inside the unit circle. */
		if (s == 0.0) {
			x = y = 0.0;
		} else {
			double factor = sqrt (-2.0 * log (s) / s);
			x *= factor, y *= factor;
		}
		secondAvailable = TRUE;
		return mean + standardDeviation * x;
	}
}
