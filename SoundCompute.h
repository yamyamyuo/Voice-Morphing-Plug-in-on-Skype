#ifndef _SOUNDCOMPUTE_H_
#define _SOUNDCOMPUTE_H_

#include "Structure.h"
//SoundCompute.cpp
MorphingFactor Factor_create(double fm = 1, double pm = 1, double prm = 1, double dm = 1);

LocationTuple Location_create(double x = 0, double y = 0);

MorphingFactor location_to_factors(LocationTuple tuple);

Sound Sound_create(double xmin, double xmax, long nx, double dx, double x1, 
		   long numOfChannels = 1);

Sound Sound_createSimple(double duration, double samplingFrequency, long numOfChannels = 1);

Sound Sound_upsample (Sound me);

Sound Sound_resample(Sound me, double samplingFrequency, long precision);

void Sound_getMaximumAndX(Sound me, double xmin, double xmax, long channel, int interpolation,
	       double *return_maximum, double *return_xOfMaximum);	

void Sound_getMinimumAndX(Sound me, double xmin, double xmax, long channel, int interpolation,
	       double *return_minimum, double *return_xOfMinimum);
						  
void Sound_getMaximumAndXAndChannel (Sound me, double xmin, double xmax, int interpolation,
	       double *return_maximum, double *return_xOfMaximum, long *return_channelOfMaximum);
	
void Sound_getMinimumAndXAndChannel(Sound me, double xmin, double xmax, int interpolation, 
           double *return_minimum, double *return_xOfMinimum, long *return_channelOfMinimum);

double Sound_getMaximum(Sound me, double xmin, double xmax, int interpolation);									
double Sound_getMinimum(Sound me, double xmin, double xmax, int interpolation);

double Sound_getAbsoluteExtremum(Sound me, double xmin, double xmax, int interpolation);
		  
long Sampled_getWindowSamples (Sound me, double xmin, double xmax, long *ixmin, long *ixmax);
double Sound_getValueAtX (Sound me, double x, long ilevel, int interpolation);

VoicedUnvoiced GeVicedUnvoiced(Pitch pitch, PointProcess points, 
							double periodsPerWindow = 3.0, double minimumPitch = 75);

void Sound_substructMean(Sound me);

void Sound_overrideSamplingFrequency (Sound me, double rate);

Sound Sound_changespeaker(Sound me, double pitchMin, double pitchMax, 
	                       double formantMultiplier , // > 0
	                       double pitchMultiplier, // > 0
						   double pitchRaneMultiplier, // any number
						   double durationMultiplier  // > 0
						   );

Sound Sound_and_pitch_changespeaker(Sound me, Pitch him, 
	                       double formantMultiplier , // > 0
	                       double pitchMultiplier,     // > 0
						   double pitchRaneMultiplier, // any number
						   double durationMultiplier  // > 0
						   );

/* Outphased */
Sound Sound_changeGender_old (Sound me, double fmin, double fmax, double formantRatio,
	double new_pitch, double pitchRangeFactor, double durationFactor);

Sound Sound_and_Pitch_changeGender_old (Sound me, Pitch him, double formantRatio,
	double new_pitch, double pitchRangeFactor, double durationFactor);

             //long AnyTier_timeToLowIndex (I, double time);

double RealTier_getValueAtTime (RealTier me, double t);

double RealTier_getArea (RealTier me, double tmin, double tmax);
/* Inside points: linear intrapolation. */
/* Outside points: constant extrapolation. */
/* No points: NUMundefined. */


//TD-PSOLA procedure

Sound Sound_Point_Pitch_Duration_to_Sound (Sound me, PointProcess pulses,
	PitchTier pitch, DurationTier duration, double maxT,double formantMultiplier);

//Manipulation.cpp
void copyRise (Sound me, double tmin, double tmax, Sound thee, double tmaxTarget);

void copyFall (Sound me, double tmin, double tmax, Sound thee, double tminTarget);

void copyBell (Sound me, double tmid, double leftWidth, double rightWidth, 
			  Sound thee, double tmidTarget);

void copyBell2 (Sound me, PointProcess source, long isource, double leftWidth, 
	            double rightWidth,Sound thee, double tmidTarget, double maxT);
void copyFall2(Sound me, PointProcess source, long isource, double leftWidth, double rightWidth,
	            Sound thee, double tmidTarget, double maxT);
void copyRise2 (Sound me, PointProcess source, long isource, double leftWidth, double rightWidth,
	            Sound thee, double tmidTarget, double maxT);
#endif
