//Pitch_to_PointProcess.h

#include "Structure.h"

//Pitch_to_PointProcess.cpp
static double findExtremum_3 (double *channel1_base, double *channel2_base, long d, long n, int includeMaxima, int includeMinima);

static double Sound_findExtremum (Sound me, double tmin, double tmax, int includeMaxima, int includeMinima);

static double Sound_findMaximumCorrelation (Sound me, double t1, double windowLength, double tmin2, double tmax2, double *tout, double *peak);

long PointProcess_getLowIndex (PointProcess me, double t);

void PointProcess_addPoint (PointProcess me, double t);

void PointProcess_init (I, double tmin, double tmax, long initialMaxnt);

PointProcess PointProcess_create (double tmin, double tmax, long initialMaxnt);

int Pitch_getVoicedIntervalAfter (Pitch me, double after, double *tleft, double *tright);

PointProcess Sound_Pitch_to_PointProcess_cc (Sound sound, Pitch pitch);

long PointProcess_getNearestIndex (PointProcess me, double t);
