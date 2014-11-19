#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_
/*
 * Sructure.h
 *  Define some useful structure for sound 
 *  Includeing contend:
 *  Schema for describing the structure and Intermediate structure for handling some result 
 *  Referance: Praat 1992-2011 by Paul Boersma
 *  Version: 1.0
 *  Author: Wenjie Feng
 *  Time: 2012.09.28
*/

#include <string>

typedef void *Any; /* Prevent compile-time type checking. */

/*
	Use the macros 'I' and 'thou' for objects in the formal parameter lists
	(if the explicit type cannot be used).
	Use the macros 'iam' and 'thouart' as the first declaration in a function definition.
	After this, the object 'me' or 'thee' has the right class (for the compiler),
	so that you can use the macros 'my' and 'thy' to refer to members.
	Example: int Person_getAge (I) { iam (Person); return my age; }
*/


#define I Any void_me
#define thou  Any void_thee
#define iam(klas)  klas me = (klas)void_me
#define thouart(klas)  klas thee = (klas)void_thee
#define my  me->
#define thy  thee->
#define his  him->
#define FALSE 0
#define TRUE 1

#define Vector_CHANNEL_AVERAGE  0
#define Vector_CHANNEL_1  1
#define Vector_CHANNEL_2  2
#define Vector_VALUE_INTERPOLATION_NEAREST  0
#define Vector_VALUE_INTERPOLATION_LINEAR  1
#define Vector_VALUE_INTERPOLATION_CUBIC  2
#define Vector_VALUE_INTERPOLATION_SINC70  3
#define Vector_VALUE_INTERPOLATION_SINC700  4

struct structSound {
    double xmin;     // start time(seconds).
	double xmax;     // end time(secondes).
	long nx;         // Number of samples.
	double dx;       // Sampling period (second).
	double x1;       // Time of the first sample (seconds).
	long ny;         // channel number
	double **z;      // z[ny][nx]; // Amplitude
};

typedef struct structSound* Sound;//Ö¸Õë

struct structPitch_Candidate {
    double frequency;
	double strength;
};

typedef struct structPitch_Candidate* Pitch_Candidate;

struct structPitch_Frame {
    double intensity;
	long nCandidates;
	Pitch_Candidate candidate;       // candidate[nCandidate].
};

typedef struct structPitch_Frame* Pitch_Frame;

struct structPitch {
    double xmin;           //start time(seconds)
	double xmax;           //tmax > tmin end time(seconds
	long nx;               //number of time slices   number of frames
	double dx;			   //time step(seconds)   analysis step
	double x1;             //centre of first time slice (seconds).
	double ceiling;        //candidates with a higher frequency are unvoiced.
	int maxnCandidates;    //maximum number of candidates per time slice.
	Pitch_Frame frame;     //frame[1..nx]
	/*
	  frame[1..nx].nCandidates	// Number of candidates in each time slice, including the unvoiced candidate.
	  frame[1..nx].candidate[1..nCandidates].frequency
		   The frequency of each candidate (Hz), 0 means aperiodic or silent.
		   candidate[1].frequency is the frequency of the currently best candidate.
	  frame[1..nx].candidate[1..nCandidates].strength
		   The strength of each candidate, a real number between 0 and 1:
		   0 means not periodic at all, 1 means perfectly periodic;
		   if the frequency of the candidate is 0, its strength is a real number
		   that represents the maximum periodicity that
		   can still be considered to be due to noise (e.g., 0.4).
		   candidate[1].strength is the strength of the currently best candidate.
	  frame[1..nx].intensity
		   The relative intensity of each frame, a real number between 0 and 1. 
	*/
};

typedef struct structPitch* Pitch;

struct structNUMfft_Table
{
  long n;
  double *trigcache;
  long *splitcache;
};

typedef struct structNUMfft_Table *NUMfft_Table;

struct structPointProcess {
    double xmin;    // start time (s)
	double xmax;    // end time (s)
	long maxnt;     // default 10  maximum number of pitch point
	long nt;        // number of pitch point
	double *t;      // t[maxnt]    vector of the pitch point
};

typedef struct structPointProcess *PointProcess;

struct structVoicedUnvoiced{
   bool *vuv;
   long nt;
};

typedef struct structVoicedUnvoiced *VoicedUnvoiced;

struct structSimpleDouble {
    double number;
};

typedef struct structSimpleDouble *SimpleDouble, *AnyPoint;

struct structRealPoint : structSimpleDouble{
	double value;
};

typedef struct structRealPoint *RealPoint;

struct structSortedSetOfDouble {
   long _capacity, size;
 //  bool _dontOwnItems;

   SimpleDouble *item;
   /*
      Attributes:
	     _capacity >= size    // private; grows as you add items.
	     size			      // the current number of items.
	     item [1..size]		  // the items.
   */
};

typedef struct structSortedSetOfDouble *SortedSetOfDouble;

struct structRealTier{
   double xmin;
   double xmax;
   SortedSetOfDouble points; 
};

typedef struct structRealTier *RealTier;

struct structPitchTier : structRealTier{ };

typedef struct structPitchTier *PitchTier;

struct structDurationTier : structRealTier{ };

typedef struct structDurationTier *DurationTier;


struct  structLocationTuple{
    double x;
	double y;
};

typedef struct structLocationTuple *LocationTuple;

struct structMorphingFactor{
     double formantMultiplier;
	 double pitchMultiplier;
	 double pitchRangeMultiplier;
     double durationMultiplier;
};

typedef struct structMorphingFactor *MorphingFactor;

#endif
