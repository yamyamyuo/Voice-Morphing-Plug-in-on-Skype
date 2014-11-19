#ifndef _NUM_h_
#define _NUM_h_

#include "math.h"

//      print e(1)
#define NUMe  2.7182818284590452353602874713526624977572
//      print 1/l(2)
#define NUMlog2e  1.4426950408889634073599246810018921374266
//      print l(10)/l(2)
#define NUMlog2_10  3.3219280948873623478703194294893901758648
//      print 1/l(10)
#define NUMlog10e  0.4342944819032518276511289189166050822944
//      print l(2)/l(10)
#define NUMlog10_2  0.3010299956639811952137388947244930267682
//      print l(2)
#define NUMln2  0.6931471805599453094172321214581765680755
//      print l(10)
#define NUMln10  2.3025850929940456840179914546843642076011
//      print a(1)*8
#define NUM2pi  6.2831853071795864769252867665590057683943
//      print a(1)*4
#define NUMpi  3.1415926535897932384626433832795028841972
//      print a(1)*2
#define NUMpi_2  1.5707963267948966192313216916397514420986
//      print a(1)
#define NUMpi_4  0.7853981633974483096156608458198757210493
//      print 0.25/a(1)
#define NUM1_pi  0.3183098861837906715377675267450287240689
//      print 0.5/a(1)
#define NUM2_pi  0.6366197723675813430755350534900574481378
//      print sqrt(a(1)*4)
#define NUMsqrtpi  1.7724538509055160272981674833411451827975
//      print sqrt(a(1)*8)
#define NUMsqrt2pi  2.5066282746310005024157652848110452530070
//      print 1/sqrt(a(1)*8)
#define NUM1_sqrt2pi  0.3989422804014326779399460599343818684759
//      print 1/sqrt(a(1))
#define NUM2_sqrtpi  1.1283791670955125738961589031215451716881
//      print l(a(1)*4)
#define NUMlnpi  1.1447298858494001741434273513530587116473
//      print sqrt(2)
#define NUMsqrt2  1.4142135623730950488016887242096980785697
//      print sqrt(0.5)
#define NUMsqrt1_2  0.7071067811865475244008443621048490392848
//      print sqrt(3)
#define NUMsqrt3  1.7320508075688772935274463415058723669428
//      print sqrt(5)
#define NUMsqrt5  2.2360679774997896964091736687312762354406
//      print sqrt(6)
#define NUMsqrt6  2.4494897427831780981972840747058913919659
//      print sqrt(7)
#define NUMsqrt7  2.6457513110645905905016157536392604257102
//      print sqrt(8)
#define NUMsqrt8  2.8284271247461900976033774484193961571393
//      print sqrt(10)
#define NUMsqrt10  3.1622776601683793319988935444327185337196
//      print sqrt(5)/2-0.5
#define NUM_goldenSection  0.6180339887498948482045868343656381177203
// The Euler-Mascheroni constant cannot be computed by bc.
// Instead we use the 40 digits computed by Johann von Soldner in 1809.
#define NUM_euler  0.5772156649015328606065120900824024310422
#define NUMundefined  HUGE_VAL
#define NUMdefined(x)  ((x) != NUMundefined)
#define NUMlog2(x)  (log (x) * NUMlog2e)


/***** ÔÝ´æÓë NUMsort.cpp ******/
double NUMhertzToMel (double hertz);
double NUMmelToHertz (double mel);
double NUMhertzToErb (double hertz);
double NUMhertzToMel (double hertz);
double NUMsemitonesToHertz (double semitones);
double NUMhertzToSemitones (double hertz);
double NUMerbToHertz (double erb);


/**** Sorting (NUMsort.cpp)****/
void NUMsort_d (long n, double ra []);   /* Heap sort. */

void NUMsort_i (long n, int ra []);

void NUMsort_l (long n, long ra []);

void NUMsort_str (long n, wchar_t *a []);

double NUMquantile (long n, double a [], double factor);

void NUMsort_p (long n, void *a [], int (*compare) (const void *, const void *));


/**** Interpolation and optimization (NUMInOp.cpp)****/

#define NUM_VALUE_INTERPOLATE_NEAREST  0
#define NUM_VALUE_INTERPOLATE_LINEAR  1
#define NUM_VALUE_INTERPOLATE_CUBIC  2

#define NUM_PEAK_INTERPOLATE_NONE  0
#define NUM_PEAK_INTERPOLATE_PARABOLIC  1
#define NUM_PEAK_INTERPOLATE_CUBIC  2
#define NUM_PEAK_INTERPOLATE_SINC70  3
#define NUM_PEAK_INTERPOLATE_SINC700  4
#define NUM_VALUE_INTERPOLATE_SINC70  70
#define NUM_VALUE_INTERPOLATE_SINC700  700

static double improve_evaluate (double x, void *closure);

double NUM_interpolate_sinc (double y [], long nx, double x, long interpolationDepth);

double NUMimproveExtremum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real, int isMaximum);

double NUMimproveMaximum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real);

double NUMimproveMinimum (double *y, long nx, long ixmid, int interpolation, double *ixmid_real);

/******************random part(NUMinop.cpp)***************/
#define LONG_LAG  100
#define SHORT_LAG  37
#define STREAM_SEPARATION  70
#define QUALITY  1009
#define LAG_DIFF  (LONG_LAG - SHORT_LAG)

#define repeat  do
#define until(cond)  while (!(cond))


static double randomArray [LONG_LAG];
static int randomInited = 0;
static long randomArrayPointer1, randomArrayPointer2, iquality;

void NUMrandomRestart (unsigned long seed);

double NUMrandomFraction (void);

double NUMrandomUniform (double lowest, double highest);

long NUMrandomInteger (long lowest, long highest);

double NUMrandomGauss (double mean, double standardDeviation);

#endif
