//Pitch.h
#include "Structure.h"

//Pitch.cpp

#define Pitch_LEVEL_FREQUENCY  1
#define Pitch_LEVEL_STRENGTH  2

#define Pitch_STRENGTH_UNIT_min  0
#define Pitch_STRENGTH_UNIT_AUTOCORRELATION  0
#define Pitch_STRENGTH_UNIT_NOISE_HARMONICS_RATIO  1
#define Pitch_STRENGTH_UNIT_HARMONICS_NOISE_DB  2
#define Pitch_STRENGTH_UNIT_max  2

#define Pitch_NEAREST  0
#define Pitch_LINEAR  1

#define kPitch_unit_HERTZ 0
#define kPitch_unit_HERTZ_LOGARITHMIC 1
#define kPitch_unit_MEL 2
#define kPitch_unit_LOG_HERTZ 3
#define kPitch_unit_SEMITONES_1 4
#define kPitch_unit_SEMITONES_100 5
#define kPitch_unit_SEMITONES_200 6
#define kPitch_unit_SEMITONES_440 7
#define kPitch_unit_ERB 8
#define ERB(f)  NUMhertzToErb (f)

void Pitch_Frame_init (Pitch_Frame frame, int nCandidate);

Pitch Pitch_create(double tmin, double tmax, long nt, double dt, double t1, double ceiling, int maxCandidates);

void Pitch_setCeiling (Pitch pitch, double ceiling);

int Pitch_getMaxnCandidates (Pitch me);

void Pitch_pathFinder (Pitch me, double silenceThreshold, double voicingThreshold, double octaveCost,
	                   double octaveJumpCost, double voicedUnvoicedCost, double ceiling, int pullFormants);

double Pitch_getValueAtTime (Pitch me, double time, int unit, int interpolate);

bool Pitch_isVoiced_i(Pitch me, long iframe);

double Sampled_getValueAtSample (Pitch me, long isamp, long ilevel, int unit);

double Sampled_getValueAtX (Pitch me, double x, long ilevel, int unit, int interpolate);

double v_getValueAtSample (Pitch me, long iframe, long ilevel, int unit);

double v_convertStandardToSpecialUnit (double value, long ilevel, int unit);

void Pitch_scaleDuration(Pitch pitch, double multiplier);

void Pitch_scalePitch(Pitch pitch, double multiplier);

bool intersectRangeWithDomain(Pitch me, double *x1, double *x2);

long Sampled_getWindowSamples (Pitch me, double xmin, double xmax, long *ixmin, long *ixmax);

double getQuantitle(Pitch me, double xmin, double xmax, double quantile, long ilevel, int unit);

double Pitch_getQuantile (Pitch me, double tmin, double tmax, double quantile, int unit);

Pitch Pitch_scaleTime_old (Pitch me, double scaleFactor);

double Pitch_getMinimum (Pitch me, double tmin, double tmax, int unit, int interpolate);

void Pitch_getMinimumAndTime (Pitch me, double tmin, double tmax, int unit, int interpolate,
	                          double *return_minimum, double *return_timeOfMinimum);

void Sampled_getMinimumAndX (Pitch me, double xmin, double xmax, long ilevel, int unit, int interpolate,
	double *return_minimum, double *return_xOfMinimum);

bool IntersectRangeWithDomain (Pitch me, double *x1, double *x2);
