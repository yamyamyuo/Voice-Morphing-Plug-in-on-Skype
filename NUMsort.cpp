#include "NUM.h"
#include <string.h>

double NUMhertzToMel (double hertz) {
	return hertz < 0 ? NUMundefined : 550.0 * log (1.0 + hertz / 550.0);
}

double NUMmelToHertz (double mel) {
	return mel < 0 ? NUMundefined : 550.0 * (exp (mel / 550.0) - 1);
}

double NUMhertzToErb (double hertz) {
	return hertz < 0 ? NUMundefined : 11.17 * log ((hertz + 312.0) / (hertz + 14680.0)) + 43.0;
}

double NUMsemitonesToHertz (double semitones) {
	return semitones == NUMundefined ? NUMundefined : 100.0 * exp (semitones * (NUMln2 / 12.0));
}

double NUMhertzToSemitones (double hertz) {
	return hertz <= 0.0 ? NUMundefined : 12.0 * log (hertz / 100.0) / NUMln2;
}

double NUMerbToHertz (double erb) {
	double dum = exp ((erb - 43.0) / 11.17);
	return erb < 0 ? NUMundefined : (14680.0 * dum - 312.0) / (1.0 - dum);
}


#define MACRO_NUMsort(TYPE)  { \
	long l, r, j, i; \
	TYPE k; \
	if (n < 2) return;   /* Already sorted. */ \
	/* This n<2 step is absent from Press et al.'s implementation, */ \
	/* which will therefore not terminate on if(--ir==1). */ \
	/* Knuth's initial assumption is now fulfilled: n >= 2. */ \
if (n == 2) { \
	if (a [1] > a [2]) { TYPE min = a [2]; a [2] = a [1]; a [1] = min; } \
	return; \
} \
if (n <= 12) { \
	for (i = 1; i < n; i ++) { \
		TYPE min = a [i]; \
		long pmin = i; \
		for (j = i + 1; j <= n; j ++) if (a [j] < min) { \
			min = a [j]; \
			pmin = j; \
		} \
		a [pmin] = a [i]; \
		a [i] = min; \
	} \
	return; \
} \
	l = (n >> 1) + 1; \
	r = n; \
	for (;;) { \
		if (l > 1) { \
			l --; \
			k = a [l]; \
		} else /* l == 1 */ { \
			k = a [r]; \
			a [r] = a [1]; \
			r --; \
			if (r == 1) { a [1] = k; return; } \
		} \
		j = l; \
		for (;;) { \
			i = j; \
			j = j << 1; \
			if (j > r) break; \
			if (j < r && a [j] < a [j + 1]) j ++; \
			if (k >= a [j]) break; \
			a [i] = a [j]; \
		} \
		a [i] = k; \
	} \
}

void NUMsort_d (long n, double a [])
	MACRO_NUMsort (double)

void NUMsort_i (long n, int a [])
	MACRO_NUMsort (int)

void NUMsort_l (long n, long a [])
	MACRO_NUMsort (long)

void NUMsort_str (long n, wchar_t *a []) {
	long l, r, j, i;
	wchar_t *k;
	if (n < 2) return;
	l = (n >> 1) + 1;
	r = n;
	for (;;) {
		if (l > 1) {
			l --;
			k = a [l];
		} else { 
			k = a [r];
			a [r] = a [1];
			r --;
			if (r == 1) { a [1] = k; return; }
		}
		j = l;
		for (;;) {
			i = j;
			j = j << 1;
			if (j > r) break;
			if (j < r && wcscmp (a [j], a [j + 1]) < 0) j ++;
			if (wcscmp (k, a [j]) >= 0) break;
			a [i] = a [j];
		}
		a [i] = k;
	}
}

void NUMsort_p (long n, void *a [], int (*compare) (const void *, const void *)) {
	long l, r, j, i;
	void *k;
	if (n < 2) return;
	l = (n >> 1) + 1;
	r = n;
	for (;;) {
		if (l > 1) {
			l --;
			k = a [l];
		} else { 
			k = a [r];
			a [r] = a [1];
			r --;
			if (r == 1) { a [1] = k; return; }
		}
		j = l;
		for (;;) {
			i = j;
			j = j << 1;
			if (j > r) break;
			if (j < r && compare (a [j], a [j + 1]) < 0) j ++;
			if (compare (k, a [j]) >= 0) break;
			a [i] = a [j];
		}
		a [i] = k;
	}
}

double NUMquantile (long n, double a [], double factor) {
	double place = factor * n + 0.5;
	long left = floor (place);
	if (n < 1) return 0.0;
	if (n == 1) return a [1];
	if (left < 1) left = 1;
	if (left >= n) left = n - 1;
	if (a [left + 1] == a [left]) return a [left];
	return a [left] + (place - left) * (a [left + 1] - a [left]);
}

/* End of file NUMsort.cpp */
