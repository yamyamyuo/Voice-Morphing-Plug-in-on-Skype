#ifndef _REALTIER_H
#define _REALTIER_H

#include "Structure.h"

void PitchTier_modifyExcursionRange (PitchTier me, double tmin, double tmax, double multiplier, double fref_Hz);

RealPoint RealPoint_create (double time, double value);

void RealTier_init (RealTier me, double tmin, double tmax);

DurationTier DurationTier_create (double tmin, double tmax);

DurationTier PointProcess_upto_DurationTier (PointProcess me);

PitchTier PitchTier_create (double tmin, double tmax);

long SortedSet_getposition(SortedSetOfDouble me, SimpleDouble data); 

void SortedSet_insertItem(SortedSetOfDouble me, SimpleDouble data, long position); 

void SortedSet_addItem(SortedSetOfDouble me, SimpleDouble data);

void SortedSetOfDouble_init(SortedSetOfDouble me, long initialCapacity);

SortedSetOfDouble SortedSetOfDouble_create (void);

void PitchTier_shiftFrequencies (PitchTier me, double tmin, double tmax, double shift, int units);

void PitchTier_multiplyFrequencies (PitchTier me, double tmin, double tmax, double factor);

void RealTier_addPoint(RealTier me, double t, double value);

double RealTier_getValueAtTime(RealTier me, double t);

PitchTier Pitch_to_PitchTier (Pitch me);

Pitch Pitch_PuitchTier_to_Pitch(Pitch me, PitchTier tier);

long AnyTier_timeToLowIndex (I, double time);

long AnyTier_timeToHighIndex (I, double time);

double RealTier_getMinimumValue (RealTier me);

void PitchTier_modifyRange_old (PitchTier me, double tmin, double tmax, double factor, double fmid);

#endif
